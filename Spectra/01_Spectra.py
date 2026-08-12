import Spectra
import pandas as pd
import unicodedata
import pickle
import scanpy as sc
import os


def clean_text(s):
    """
    Normalize string, remove unicode diacritics and non‑ascii characters.
    Handle pd.NA / np.nan input, return empty string for missing values.

    Parameters
    ----------
    s : str or pd.NA
        Raw input string to be sanitized

    Returns
    -------
    str
        Cleaned ascii‑only string
    """
    if pd.isna(s):
        return ""
    s = str(s).strip()
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")


def build_spectra_gene_dict(xlsx_path):
    """
    Construct Spectra‑compatible nested gene‑set dictionary from excel table.
    Excel table must contain columns: ["cell_type", "process", "gene"].
    Output dict structure: {cell_type: {process_name: [geneA, geneB, ...]}}.

    Parameters
    ----------
    xlsx_path : str
        File path to input excel (.xlsx) file

    Returns
    -------
    dict
        Nested dictionary formatted for Spectra gene_set_dictionary argument
    """
    df = pd.read_excel(xlsx_path)
    df["cell_type"] = df["cell_type"].apply(clean_text)
    df["process"] = df["process"].apply(clean_text)
    df["gene"] = df["gene"].apply(clean_text)

    df = df.dropna(subset=["cell_type", "process", "gene"])
    df = df.drop_duplicates(subset=["cell_type", "process", "gene"])

    gene_dict = {}
    for ct in df["cell_type"].unique():
        sub_df = df[df["cell_type"] == ct]
        proc_dict = {}
        for proc in sub_df["process"].unique():
            raw_genes = sub_df.loc[sub_df["process"] == proc, "gene"].tolist()
            clean_genes = [g.strip() for g in raw_genes if len(g.strip()) > 0]
            proc_dict[proc] = clean_genes
        gene_dict[ct] = proc_dict
    return gene_dict


def run_spectra_metacell(
    work_dir,
    gene_dict_pkl_name="gene_dict.pkl",
    adata_h5ad_name="seacell_metacell_logavg_filtered.h5ad",
    cell_type_key="MostAbundantCellType",
    lam=0.1,
    delta=0.001,
    kappa=None,
    rho=0.001,
    num_epochs=2000,
    n_top_vals=50,
    overlap_threshold=0.2,
    min_gs_num=3,
    use_weights=True,
    use_cell_types=True,
    use_highly_variable=False,
    clean_gs=True,
    label_factors=True
):
    """
    Wrapper for Spectra.est_spectra on SEACells metacell anndata.
    Perform gene‑dict repair, sanity checks, model fitting and export all outputs.

    Parameters
    ----------
    work_dir : str
        Working directory for input file reading and output writing
    gene_dict_pkl_name : str
        Filename of pre‑built gene‑set dictionary pickle inside work_dir
    adata_h5ad_name : str
        Filename of metacell anndata (.h5ad) inside work_dir
    cell_type_key : str
        Column name in adata.obs storing cell‑type annotations
    lam : float
        Spectra lambda: balance graph loss vs expression loss
    delta : float
        Spectra delta: lower bound for gene scaling factors
    kappa : float or None
        Background rate of 1s in adjacency graph; None → estimate from data
    rho : float or None
        Background rate of 0s in adjacency graph; None → estimate from data
    num_epochs : int
        Training epochs passed to model.train()
    n_top_vals : int
        Number of top marker genes per factor to extract
    overlap_threshold : float
        Szymkiewicz–Simpson overlap cutoff for factor labeling
    min_gs_num : int
        Minimum expressed gene count per gene‑set, active when clean_gs=True
    use_weights : bool
        Whether to use graph‑derived edge weights during training
    use_cell_types : bool
        Fit cell‑type‑specific factors using cell‑type annotations
    use_highly_variable : bool
        Restrict model to highly variable genes; set False for metacell data
    clean_gs : bool
        Enable Spectra internal gene‑set dictionary cleaning
    label_factors : bool
        Assign factor labels based on overlap with input gene‑sets

    Returns
    -------
    tuple
        (model, adata)
        model: trained SPECTRA_Model object
        adata: anndata object with Spectra results stored in .obsm / .uns
    """

    pkl_path = os.path.join(work_dir, gene_dict_pkl_name)
    pkl_fixed_path = os.path.join(work_dir, "gene_dict_fixed.pkl")
    h5ad_path = os.path.join(work_dir, adata_h5ad_name)
    out_prefix = os.path.join(work_dir, "spectra_run")

    with open(pkl_path, "rb") as f:
        gene_dict = pickle.load(f)

    adata = sc.read_h5ad(h5ad_path)

    obs_cts = set(adata.obs[cell_type_key].unique())
    dict_cts = set(gene_dict.keys())
    missing_keys = obs_cts - dict_cts
    for ck in missing_keys:
        gene_dict[ck] = dict()
    print(f"Added empty placeholder for cell_types: {missing_keys}")

    with open(pkl_fixed_path, "wb") as f:
        pickle.dump(gene_dict, f)
    print(f"Fixed gene_dict saved to: {pkl_fixed_path}")

    print("gene_dict top keys:", list(gene_dict.keys()))
    print(f"adata shape: {adata.shape}")
    print("adata obs columns:", adata.obs.columns.tolist())

    dict_cts_noglobal = set(gene_dict.keys()) - {"global"}
    print("\n==== cell type match check ====")
    print(f"dict specific cell keys: {dict_cts_noglobal}")
    print(f"adata obs cell types: {obs_cts}")
    print("dict only cell types:", dict_cts_noglobal - obs_cts)
    print("adata only cell types:", obs_cts - dict_cts_noglobal)

    all_genes_in_dict = set()
    for ctype, gs_dict in gene_dict.items():
        for gs_name, glist in gs_dict.items():
            all_genes_in_dict.update(glist)
    adata_genes = set(adata.var_names)
    overlap = all_genes_in_dict & adata_genes
    no_overlap = all_genes_in_dict - adata_genes
    print("\n==== gene overlap check ====")
    print(f"Total genes in gene_dict: {len(all_genes_in_dict)}")
    print(f"Overlap with adata.var_names: {len(overlap)}")
    print(f"Missing genes(first 20): {sorted(list(no_overlap))[:20]}")

    valid_gs = []
    invalid_gs = []
    for ctype, proc_dict in gene_dict.items():
        for proc_name, raw_glist in proc_dict.items():
            inter = [g for g in raw_glist if g in adata_genes]
            if len(inter) >= min_gs_num:
                valid_gs.append((ctype, proc_name, len(inter)))
            else:
                invalid_gs.append((ctype, proc_name, len(inter)))
    print(f"\nValid gene_set(>={min_gs_num} overlapping genes): {len(valid_gs)}")
    print(f"Invalid gene_set(<{min_gs_num} overlapping genes): {len(invalid_gs)}")
    print("invalid examples:", invalid_gs[:10])

    model = Spectra.est_spectra(
        adata=adata,
        gene_set_dictionary=gene_dict,
        use_highly_variable=use_highly_variable,
        cell_type_key=cell_type_key,
        use_weights=use_weights,
        lam=lam,
        delta=delta,
        kappa=kappa,
        rho=rho,
        use_cell_types=use_cell_types,
        n_top_vals=n_top_vals,
        label_factors=label_factors,
        overlap_threshold=overlap_threshold,
        clean_gs=clean_gs,
        min_gs_num=min_gs_num,
        num_epochs=num_epochs
    )

    factors = adata.uns['SPECTRA_factors']
    markers = adata.uns['SPECTRA_markers']
    cell_scores = adata.obsm['SPECTRA_cell_scores']

    print(f"\nfactors shape: {factors.shape}")
    print(f"markers shape: {markers.shape}")
    print(f"cell_scores shape: {cell_scores.shape}")

    pd.DataFrame(
        factors,
        index=[f"Factor_{i+1}" for i in range(factors.shape[0])],
        columns=adata.var_names
    ).to_csv(f"{out_prefix}_factors.csv", encoding="utf-8-sig")

    pd.DataFrame(
        cell_scores,
        index=adata.obs_names,
        columns=[f"Factor_{i+1}" for i in range(cell_scores.shape[1])]
    ).to_csv(f"{out_prefix}_cell_scores.csv", encoding="utf-8-sig")

    marker_df = pd.DataFrame(markers)
    marker_df.index = [f"Factor_{i+1}" for i in range(marker_df.shape[0])]
    marker_df.to_csv(f"{out_prefix}_markers.csv", encoding="utf-8-sig")

    adata.write_h5ad(os.path.join(work_dir, "seacell_metacell_spectra.h5ad"))
    print("\nAll output saved, finished.")

    return model, adata


# -----------------------------------------------------------------------------
# Example usage (comment / uncomment as needed)
# -----------------------------------------------------------------------------
# if __name__ == "__main__":
#     # Step1: build gene dict from excel (local desktop)
#     # xlsx_file = r"E:/AID cohort/code/NMF/Spectra/mydict.xlsx"
#     # gd = build_spectra_gene_dict(xlsx_file)
#     # with open(r"E:/AID cohort/code/NMF/Spectra/gene_dict.pkl", "wb") as fout:
#     #     pickle.dump(gd, fout)
#
#     # Step2: run spectra on cluster
#     # wd = "/mnt/data/work/masy/scAltas/2026scRNAcohort"
#     # spectra_model, adata_out = run_spectra_metacell(work_dir=wd)
