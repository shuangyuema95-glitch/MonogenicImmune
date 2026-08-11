#!/usr/bin/env python
"""
Aggregate single‑cell data to SEACell metacells, reproduce Nature Medicine aggregation logic.
- Filter metacells: retain only lowPurity == False
- adata.X: arithmetic mean of log‑normalized single‑cell expression (for Spectra input)
- adata.raw.X: sum of raw UMI counts per metacell (for DE / DESeq2 / decoupler)
Input:
    seacell_cell_mapping_purityFilt.pkl: cell‑to‑SEACell mapping table
    seacell_metacell_stats.csv: metacell‑level QC metadata
    716_pbmc1.h5ad: original single‑cell AnnData
    retainGenes.csv: target gene list to subset
Output:
    seacell_metacell_logavg_filtered.h5ad: metacell‑level AnnData
    check.txt: runtime log and validation report
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
from tqdm import tqdm


def clean_gene(g: str) -> str:
    """Strip whitespace for gene symbol matching."""
    return str(g).strip()


def write_log(log_str: str, log_path: str):
    """Append message to log file and stdout."""
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(log_str + "\n")
    print(log_str)


def aggregate_metacells(
    base_dir: str,
    pkl_filt_name: str = "seacell_cell_mapping_purityFilt.pkl",
    meta_csv_name: str = "seacell_metacell_stats.csv",
    h5ad_name: str = "716_pbmc1.h5ad",
    gene_list_name: str = "retainGenes.csv",
    out_h5ad_name: str = "seacell_metacell_logavg_filtered.h5ad",
    log_name: str = "check.txt"
):
    """
    Perform SEACell metacell aggregation with strict QC filtering and validation.

    Parameters
    ----------
    base_dir : str
        Root working directory containing SEACELLS subfolder and input h5ad.
    pkl_filt_name : str
        Filename of cell‑to‑SEACell mapping pickle inside SEACELLS folder.
    meta_csv_name : str
        Filename of metacell statistics csv inside SEACELLS folder.
    h5ad_name : str
        Original single‑cell h5ad filename under base_dir.
    gene_list_name : str
        One‑column csv of target gene symbols under base_dir.
    out_h5ad_name : str
        Output metacell‑level h5ad filename inside SEACELLS folder.
    log_name : str
        Runtime log filename inside SEACELLS folder.
    """
    seacell_dir = os.path.join(base_dir, "SEACELLS")
    fp_pkl_filt = os.path.join(seacell_dir, pkl_filt_name)
    fp_meta_csv = os.path.join(seacell_dir, meta_csv_name)
    fp_h5ad = os.path.join(base_dir, h5ad_name)
    fp_gene = os.path.join(base_dir, gene_list_name)
    out_h5ad = os.path.join(seacell_dir, out_h5ad_name)
    check_log = os.path.join(seacell_dir, log_name)

    open(check_log, "w", encoding="utf-8").close()

    df_map = pd.read_pickle(fp_pkl_filt)
    df_meta = pd.read_csv(fp_meta_csv)
    df_meta = df_meta.set_index("SEACell")

    write_log("========== Metacell aggregation & validation log ==========", check_log)
    write_log(f"SEACell mapping pkl shape: {df_map.shape}", check_log)
    write_log(f"metacell stats csv shape: {df_meta.shape}", check_log)
    write_log(f"df_meta columns: {df_meta.columns.tolist()}", check_log)

    assert "MostAbundantCellType" in df_meta.columns, "Missing MostAbundantCellType in metadata csv"
    assert "lowPurity" in df_meta.columns, "Missing lowPurity in metadata csv"

    # filter valid metacells: lowPurity == False
    df_meta_filtered = df_meta[df_meta["lowPurity"] == False].copy()
    valid_seacell = set(df_meta_filtered.index)
    write_log(f"\nTotal raw metacells: {df_meta.shape[0]}", check_log)
    write_log(f"Retained (lowPurity=False): {len(valid_seacell)}", check_log)
    write_log(f"Discarded (lowPurity=True): {df_meta.shape[0] - len(valid_seacell)}", check_log)

    # subset cell‑mapping table to valid metacells only
    df_map = df_map[df_map["SEACell"].isin(valid_seacell)].copy()
    seacell_ids = sorted(df_map["SEACell"].unique())
    write_log(f"Metacells ready for aggregation: {len(seacell_ids)}", check_log)
    assert len(seacell_ids) == len(valid_seacell), "Mismatch: metacell IDs between pkl and csv"

    obs_df = df_meta_filtered.reindex(seacell_ids).copy()
    write_log(f"obs_df columns: {obs_df.columns.tolist()}", check_log)
    assert "MostAbundantCellType" in obs_df.columns, "MostAbundantCellType lost in obs_df"

    # load target gene list
    gene_df = pd.read_csv(fp_gene, header=None)
    keep_raw = gene_df.iloc[:, 0].dropna().tolist()
    keep_clean = [clean_gene(g) for g in keep_raw]
    keep_set = set(keep_clean)
    write_log(f"\nretainGenes raw count: {len(keep_raw)}, cleaned unique: {len(keep_set)}", check_log)

    # load single‑cell anndata and raw layer
    adata = sc.read_h5ad(fp_h5ad)
    raw_adata = adata.raw.to_adata()
    raw_adata.var_names = adata.var_names.copy()

    adata_gene_clean = [clean_gene(g) for g in adata.var_names]
    adata.var_names = adata_gene_clean
    raw_adata.var_names = adata_gene_clean
    adata_gene_set = set(adata_gene_clean)

    common_genes = sorted(list(keep_set & adata_gene_set))
    adata = adata[:, common_genes].copy()
    raw_adata = raw_adata[:, common_genes].copy()
    write_log(f"Intersection genes with adata: {len(common_genes)}", check_log)

    only_keep = keep_set - adata_gene_set
    only_adata = adata_gene_set - keep_set
    write_log(f"\n[Gene matching summary]", check_log)
    write_log(f"Genes in retainGenes but missing in adata: {len(only_keep)}", check_log)
    if len(only_keep) > 0:
        write_log(f"First 20 missing genes: {list(only_keep)[:20]}", check_log)
    write_log(f"Genes present in adata but not in retainGenes: {len(only_adata)}", check_log)

    # filter barcodes to those present in filtered mapping table
    barcodes_keep = df_map.index.tolist()
    adata = adata[barcodes_keep, :].copy()
    raw_adata = raw_adata[barcodes_keep, :].copy()

    X_list = []
    X_raw_list = []
    write_log("\nStart aggregating per‑metacell expression ...", check_log)
    for sid in tqdm(seacell_ids):
        cb_sub = df_map.loc[df_map["SEACell"] == sid].index.tolist()
        sub_X = adata[cb_sub, :].X.toarray()
        sub_X_raw = raw_adata[cb_sub, :].X.toarray()
        avg_log = np.mean(sub_X, axis=0)
        sum_raw = np.sum(sub_X_raw, axis=0)
        X_list.append(avg_log)
        X_raw_list.append(sum_raw)

    X_metacell = np.vstack(X_list)
    Xraw_metacell = np.vstack(X_raw_list)

    var_df = adata.var.drop(columns=[c for c in adata.var.columns if "_index" in str(c)], errors="ignore")
    meta_adata = sc.AnnData(X=X_metacell, obs=obs_df, var=var_df)
    meta_raw = sc.AnnData(X=Xraw_metacell, var=var_df.copy())
    meta_adata.raw = meta_raw

    meta_adata.write_h5ad(out_h5ad)
    write_log(f"\nOutput h5ad: {out_h5ad}", check_log)
    write_log(f"Final metacell shape n_obs={meta_adata.n_obs}, n_vars={meta_adata.n_vars}", check_log)

    write_log("\n==================== Post‑run validation ====================", check_log)
    write_log("obs columns: " + str(meta_adata.obs.columns.tolist()), check_log)
    write_log("\nMostAbundantCellType counts (required for Spectra cell_type_key):", check_log)
    write_log(str(meta_adata.obs["MostAbundantCellType"].value_counts()), check_log)
    write_log("\nlowPurity distribution (must be all False):", check_log)
    write_log(str(meta_adata.obs["lowPurity"].value_counts()), check_log)

    write_log(f"\nraw layer exists: {meta_adata.raw is not None}, raw shape {meta_adata.raw.shape}", check_log)
    write_log("\nFirst SEACell first‑10 log‑mean values:", check_log)
    write_log(str(meta_adata.X[0, :10]), check_log)
    write_log("Corresponding raw summed counts:", check_log)
    write_log(str(meta_adata.raw.X[0, :10]), check_log)

    write_log("\nFirst 10 SEACell IDs:", check_log)
    write_log(str(meta_adata.obs.index[:10].tolist()), check_log)
    write_log("\n========== Aggregation & validation complete ==========\n\n", check_log)


if __name__ == "__main__":
    # Modify base_dir for your environment
    BASE_DIR = "/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50"
    aggregate_metacells(base_dir=BASE_DIR)