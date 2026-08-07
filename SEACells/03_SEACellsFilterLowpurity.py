import os
import glob
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# Full gene-pathway mapping dictionary
gene_list = {
    "NFKB pathway": ["TNFAIP3", "IRAK2", "NOD2", "RELA"],
    "Uncategorized": ["ADA2", "ELF4", "CSF3R", "STAT4"],
    "Osteoclast function": ["ACP5", "OGFRL1"],
    "Ca2+_flux-PLC": ["PLCG2", "PLCG1"],
    "Endolysosomal_nucleic_acid_sensing": ["UNC93B1", "TLR8", "TLR7", "TLR1", "PLD4", "TRAF3"],
    "Arachidonic acid metabolism": ["TBXAS1"],
    "Inflammasome_IL-1β": ["IL1R1", "NLRC4", "NLRP3", "PSTPIP1", "CDC42", "LPIN2"],
    "Immune_metabolic": ["LACC1", "SLC7A7"],
    "PP": ["polygenic"],
    "HC": ["wide-type"],
    "Inborn_errors_of_cell_death": ["RIPK1", "TNFRSF1A", "OTULIN", "RNF31"],
    "Negative_regulation_of_IFN-I": ["ISG15", "USP18", "SOCS1"],
    "Cytoskeleton_and_small_GTPase": ["GNAI2", "KRAS"],
    "TBK1_IRF3": ["STING", "TREX1", "IFIH1", "COPA"],
    "Protein_homeostasis": ["UBA1", "PSMD12"]
}
gene_to_condition0 = {}
for pathway, genenames in gene_list.items():
    for g in genenames:
        gene_to_condition0[g] = pathway

def assign_condition(row):
    """
    Assign condition label to each cell row
    Priority: sample label match > gene pathway match
    Params:
        row: single row from adata.obs DataFrame, contains 'samples' and 'gene' columns
    Return:
        str: condition tag (SLE/JIA/KD/HC/pathway name) or None
    """
    sam = str(row["samples"])
    if "SLE" in sam:
        return "SLE"
    elif "JIA" in sam:
        return "JIA"
    elif "KD" in sam:
        return "KD"
    elif "HC" in sam:
        return "HC"
    else:
        g = row["gene"]
        return gene_to_condition0.get(g, None)

def concat_seacell_mapping(seacell_dir, main_adata):
    """
    Merge all *.seacell.pkl soft assignment files and merge adata obs metadata
    Params:
        seacell_dir: str, folder path storing SEACell output pkl files
        main_adata: sc.AnnData, full annotated single-cell object with obs metadata
    Return:
        pd.DataFrame: cell-level soft SEACell mapping with full obs columns
    """
    pkl_paths = glob.glob(os.path.join(seacell_dir, "*.seacell.pkl"))
    df_list = []
    for p in tqdm(pkl_paths):
        df = pd.read_pickle(p)
        df_list.append(df)
    df_all = pd.concat(df_list)
    obs_full = main_adata.obs.copy()
    df_all = df_all.merge(obs_full, left_index=True, right_index=True, how="left")
    return df_all

def compute_metacell_stats(df_cell, annotation_col="Level1", sample_col="samples", disease_col="status", cond_col="condition", purity_threshold=0.75):
    """
    Calculate aggregated metacell-level QC statistics from cell-level mapping table
    Params:
        df_cell: pd.DataFrame, full cell-level SEACell soft assignment table
        annotation_col: str, cell type column name in obs (Level1)
        sample_col: str, sample identifier column name
        disease_col: str, disease status column name
        cond_col: str, condition grouping column (SLE/JIA/KD/HC/pathway)
        purity_threshold: float, cutoff for lowPurity flag (metacell dominant cell fraction < threshold = lowPurity=True)
    Return:
        metacell_df: pd.DataFrame, one row per SEACell with all aggregated QC metrics
        ct_prop_df: pd.DataFrame, cell type fraction matrix for each SEACell
    """
    df_tmp = df_cell.copy()
    df_tmp[annotation_col] = df_tmp[annotation_col].astype("category")

    ct_counts = df_tmp.groupby(["SEACell", annotation_col]).size().unstack(fill_value=0)
    ct_total = ct_counts.sum(axis=1)
    ct_frac = ct_counts.div(ct_total, axis=0)

    max_frac = ct_frac.max(axis=1)
    max_ct = ct_frac.idxmax(axis=1)

    metacell_df = pd.DataFrame({
        "SEACell": max_frac.index,
        "MostAbundantCellType": max_ct.values,
        "MostAbundantCellType_proportion": max_frac.values,
        "n_single_cell": ct_total.values
    }).set_index("SEACell")

    sample_map = df_tmp.groupby("SEACell")[sample_col].first()
    disease_map = df_tmp.groupby("SEACell")[disease_col].first()
    cond_map = df_tmp.groupby("SEACell")[cond_col].first()
    compactness_mean = df_tmp.groupby("SEACell")["compactness"].mean()
    separation_mean = df_tmp.groupby("SEACell")["separation"].mean()

    metacell_df[sample_col] = sample_map
    metacell_df[disease_col] = disease_map
    metacell_df[cond_col] = cond_map
    metacell_df["compactness_mean"] = compactness_mean
    metacell_df["separation_mean"] = separation_mean
    metacell_df["lowPurity"] = metacell_df["MostAbundantCellType_proportion"] < purity_threshold

    ct_prop_df = ct_frac.copy()
    return metacell_df, ct_prop_df

def plot_seacell_full_qc(df_cell, metacell_df, ct_prop_df, out_pdf, annotation_col="Level1", disease_col="status", cond_col="condition", purity_threshold=0.75):
    """
    Generate full pre-filter & post-filter QC PDF figures, export each subplot data to separate CSV for R plotting
    Params:
        df_cell: pd.DataFrame, full cell-level SEACell mapping table
        metacell_df: pd.DataFrame, aggregated metacell QC stats
        ct_prop_df: pd.DataFrame, cell type fraction per SEACell
        out_pdf: str, output path for pre-filter summary PDF
        annotation_col: str, dominant cell type column name
        disease_col: str, disease status grouping column
        cond_col: str, condition grouping column
        purity_threshold: float, purity cutoff for filtering metacells
    Return:
        metacell_filt: pd.DataFrame, filtered high-purity metacell subset
        ct_prop_filt: pd.DataFrame, cell type fraction matrix for filtered metacells
    """
    plt.rcParams["font.size"] = 8
    out_csv_dir = "./SEACELLS/plot_csv"
    os.makedirs(out_csv_dir, exist_ok=True)

    # Pre-filter figure: 5 rows, 2 columns, total 10 subplots
    fig = plt.figure(figsize=(18, 26))

    # ax1: Fraction of SEACell_non_trivial_assig
    ax1 = plt.subplot(5, 2, 1)
    assign_counts = df_cell.value_counts("SEACell_non_trivial_assig", normalize=True).sort_index()
    assign_counts.plot.bar(ax=ax1)
    ax1.set_title("Fraction of SEACell_non_trivial_assig")
    ax1.set_ylabel("Fraction cells")
    assign_counts.reset_index().to_csv(os.path.join(out_csv_dir, "ax1_non_trivial_frac.csv"), index=False)

    # ax2: Histogram of cell count per SEACell
    ax2 = plt.subplot(5, 2, 2)
    seacell_cell_count = df_cell["SEACell"].value_counts().reset_index()
    seacell_cell_count.columns = ["SEACell", "cell_count"]
    sns.histplot(data=seacell_cell_count, x="cell_count", bins=75, ax=ax2)
    ax2.set_xlabel("Number of Cells per SEACell")
    ax2.set_ylabel("Number of SEACells")
    seacell_cell_count.to_csv(os.path.join(out_csv_dir, "ax2_seacell_cell_count.csv"), index=False)

    # ax3: Violin plot n_single_cell per SEACell grouped by dominant cell type (pre-filter)
    ax3 = plt.subplot(5, 2, 3)
    sns.violinplot(data=metacell_df, x="MostAbundantCellType", y="n_single_cell", ax=ax3)
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha="right")
    ax3.set_title("n_single_cell per SEACell (pre‑filter)")
    metacell_df[["MostAbundantCellType", "n_single_cell"]].to_csv(os.path.join(out_csv_dir, "ax3_violin_ct_pre.csv"), index=False)

    # ax4: Mean cell-type composition stacked bar per dominant cell type (pre-filter)
    ax4 = plt.subplot(5, 2, 4)
    ct_prop_dist = ct_prop_df.groupby(metacell_df["MostAbundantCellType"], observed=False).mean()
    ct_prop_dist.plot(kind="bar", stacked=True, ax=ax4)
    ax4.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    ax4.set_title("Mean cell‑type composition per SEACell (pre‑filter)")
    ct_prop_dist.reset_index().to_csv(os.path.join(out_csv_dir, "ax4_ct_mean_prop.csv"), index=False)

    # ax5: Absolute SEACell count stacked bar grouped by cell type + lowPurity
    ax5 = plt.subplot(5, 2, 5)
    grouped_ct = metacell_df.groupby(["MostAbundantCellType", "lowPurity"], observed=False).size().unstack(fill_value=0)
    grouped_ct.plot(kind="bar", stacked=True, ax=ax5, legend=False)
    ax5.set_ylabel("Absolute SEACell number")
    grouped_ct.reset_index().to_csv(os.path.join(out_csv_dir, "ax5_ct_type_abs.csv"), index=False)

    # ax6: Normalized SEACell frequency stacked bar grouped by cell type + lowPurity
    ax6 = plt.subplot(5, 2, 6)
    norm_ct = grouped_ct.div(grouped_ct.sum(axis=1), axis=0)
    norm_ct.plot(kind="bar", stacked=True, ax=ax6, legend=False)
    ax6.set_ylabel("SEACell frequency")
    norm_ct.reset_index().to_csv(os.path.join(out_csv_dir, "ax6_ct_type_norm.csv"), index=False)

    # ax7: Absolute SEACell count stacked bar grouped by disease status + lowPurity
    ax7 = plt.subplot(5, 2, 7)
    grouped_disease = metacell_df.groupby([disease_col, "lowPurity"], observed=False).size().unstack(fill_value=0)
    grouped_disease.plot(kind="bar", stacked=True, ax=ax7, legend=False)
    ax7.set_ylabel("Absolute SEACell number")
    grouped_disease.reset_index().to_csv(os.path.join(out_csv_dir, "ax7_disease_abs.csv"), index=False)

    # ax8: Normalized SEACell frequency stacked bar grouped by disease status + lowPurity
    ax8 = plt.subplot(5, 2, 8)
    norm_disease = grouped_disease.div(grouped_disease.sum(axis=1), axis=0)
    norm_disease.plot(kind="bar", stacked=True, ax=ax8, legend=False)
    ax8.set_ylabel("SEACell frequency")
    norm_disease.reset_index().to_csv(os.path.join(out_csv_dir, "ax8_disease_norm.csv"), index=False)

    # ax9: Absolute SEACell count stacked bar grouped by condition + lowPurity
    ax9 = plt.subplot(5, 2, 9)
    grouped_cond = metacell_df.groupby([cond_col, "lowPurity"], observed=False).size().unstack(fill_value=0)
    grouped_cond.plot(kind="bar", stacked=True, ax=ax9, legend=False)
    ax9.set_ylabel("Absolute SEACell number")
    ax9.set_xlabel("Condition")
    grouped_cond.reset_index().to_csv(os.path.join(out_csv_dir, "ax9_condition_abs.csv"), index=False)

    # ax10: Normalized SEACell frequency stacked bar grouped by condition + lowPurity
    ax10 = plt.subplot(5, 2, 10)
    norm_cond = grouped_cond.div(grouped_cond.sum(axis=1), axis=0)
    norm_cond.plot(kind="bar", stacked=True, ax=ax10, legend=False)
    ax10.set_ylabel("SEACell frequency")
    ax10.set_xlabel("Condition")
    norm_cond.reset_index().to_csv(os.path.join(out_csv_dir, "ax10_condition_norm.csv"), index=False)

    plt.tight_layout()
    plt.savefig(out_pdf, dpi=150)
    plt.close()

    # Post-filter figure: 3 rows, 1 column, total 3 subplots
    fig2 = plt.figure(figsize=(14, 12))
    metacell_filt = metacell_df[~metacell_df["lowPurity"]].copy()
    ct_prop_filt = ct_prop_df.loc[metacell_filt.index, :]
    metacell_filt.to_csv(os.path.join(out_csv_dir, "postfilter_metacell_stats.csv"), index=False)

    # ax_a: Violin plot n_single_cell grouped by condition (post-filter)
    ax_a = plt.subplot(3, 1, 1)
    sns.violinplot(data=metacell_filt, x=cond_col, y="n_single_cell", ax=ax_a)
    ax_a.set_xticklabels(ax_a.get_xticklabels(), rotation=45, ha="right")
    ax_a.set_title("n_single_cell per SEACell (post‑filter, group by Condition)")
    metacell_filt[[cond_col, "n_single_cell"]].to_csv(os.path.join(out_csv_dir, "post_axa_violin_cond.csv"), index=False)

    # ax_b: Violin plot n_single_cell grouped by dominant cell type (post-filter)
    ax_b = plt.subplot(3, 1, 2)
    sns.violinplot(data=metacell_filt, x="MostAbundantCellType", y="n_single_cell", ax=ax_b)
    ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=45, ha="right")
    ax_b.set_title("n_single_cell per SEACell (post‑filter purity >= {})".format(purity_threshold))
    metacell_filt[["MostAbundantCellType", "n_single_cell"]].to_csv(os.path.join(out_csv_dir, "post_axb_violin_ct.csv"), index=False)

    # ax_c: Mean cell-type composition stacked bar (post-filter)
    ax_c = plt.subplot(3, 1, 3)
    ct_prop_dist_filt = ct_prop_filt.groupby(metacell_filt["MostAbundantCellType"], observed=False).mean()
    ct_prop_dist_filt.plot(kind="bar", stacked=True, ax=ax_c)
    ax_c.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    ax_c.set_title("Mean cell‑type composition per SEACell (post‑filter)")
    ct_prop_dist_filt.reset_index().to_csv(os.path.join(out_csv_dir, "post_axc_ct_mean_prop.csv"), index=False)

    plt.tight_layout()
    out_pdf2 = out_pdf.replace(".pdf","_postfilter.pdf")
    plt.savefig(out_pdf2, dpi=150)
    plt.close()
    return metacell_filt, ct_prop_filt

def main():
    """
    Main execution pipeline: load adata, assign condition, concat SEACell mapping, compute stats, export CSV, plot QC figures
    Fixed global parameters defined inside function
    """
    seacell_output_dir = "./SEACELLS"
    h5ad_path = "716_pbmc1.h5ad"
    out_cell_mapping_pkl = "./SEACELLS/seacell_cell_mapping.pkl"
    out_cell_csv = "./SEACELLS/cell_level_seacell.csv"
    out_meta_stats_csv = "./SEACELLS/seacell_metacell_stats.csv"
    out_qc_pdf = "./SEACELLS/seacell_qc_summary.pdf"
    out_filtered_cell_mapping_pkl = "./SEACELLS/seacell_cell_mapping_purityFilt.pkl"

    purity_threshold = 0.75
    annotation_col = "Level1"
    sample_col = "samples"
    disease_col = "status"
    cond_col = "condition"

    # Load anndata & generate condition metadata
    adata_full = sc.read_h5ad(h5ad_path)
    obs = adata_full.obs.copy()
    obs["condition"] = obs.apply(assign_condition, axis=1)
    adata_full.obs["condition"] = obs["condition"]
    adata_full = adata_full[adata_full.obs["condition"].notna()].copy()

    # Merge all SEACell soft assignment files
    df_mapping = concat_seacell_mapping(seacell_output_dir, adata_full)
    print(f"Total mapped cells: {df_mapping.shape[0]}")
    total_seacell = df_mapping["SEACell"].nunique()
    print(f"Total generated SEACells: {total_seacell}")

    # Export cell-level mapping files
    df_mapping.to_pickle(out_cell_mapping_pkl)
    df_mapping.to_csv(out_cell_csv, index=False)
    print(f"Cell‑level raw mapping saved: {out_cell_mapping_pkl}, {out_cell_csv}")

    # Calculate metacell aggregated statistics
    metacell_df, ct_prop_df = compute_metacell_stats(
        df_mapping,
        annotation_col=annotation_col,
        sample_col=sample_col,
        disease_col=disease_col,
        cond_col=cond_col,
        purity_threshold=purity_threshold
    )

    # Export full metacell stats table
    metacell_df.to_csv(out_meta_stats_csv, index=False)
    print(f"Metacell stats csv saved: {out_meta_stats_csv}")

    # Draw QC figures & export per-subplot CSV data
    metacell_filt, ct_prop_filt = plot_seacell_full_qc(
        df_mapping,
        metacell_df,
        ct_prop_df,
        out_qc_pdf,
        annotation_col=annotation_col,
        disease_col=disease_col,
        cond_col=cond_col,
        purity_threshold=purity_threshold
    )

    # Print filtering summary
    print("\n==== Metacell purity summary ====")
    print(f"Total metacells: {metacell_df.shape[0]}")
    print(f"Pass purity (>= {purity_threshold}): {metacell_filt.shape[0]}")
    print(f"Filtered out low‑purity: {(metacell_df['lowPurity']).sum()}")

    # Export filtered cell-SEACell mapping
    seacells_keep = set(metacell_filt.index)
    df_mapping_filt = df_mapping[df_mapping["SEACell"].isin(seacells_keep)].copy()
    df_mapping_filt.to_pickle(out_filtered_cell_mapping_pkl)
    print(f"Purity‑filtered cell‑SEACell mapping saved: {out_filtered_cell_mapping_pkl}")

if __name__ == "__main__":
    main()