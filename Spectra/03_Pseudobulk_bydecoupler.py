# -*- coding: utf-8 -*-
"""
Pseudobulk aggregation by sample + Level2 (Nature Medicine ULM pipeline)

dc.pp.pseudobulk: use adata.X (log-normalized), compute group mean
post-aggregation: per-gene z-score (zero mean, unit variance)
subset to retainGenes before aggregation, filter groups with cell_count >= 5
output h5ad for downstream decoupler ULM analysis
"""

import os
import numpy as np
import scanpy as sc
import pandas as pd
import decoupler as dc


def load_gene_list(fp):
    """Load one-gene-per-line text file with explicit UTF-8 encoding."""
    with open(fp, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def main():
    work_dir = "/mnt/data/work/masy/scAltas/2026scRNAcohort"
    raw_h5ad = os.path.join(work_dir, "716_pbmc1.h5ad")
    gene_list_fp = os.path.join(work_dir, "retainGenes.txt")
    pb_h5ad = os.path.join(work_dir, "pseudobulk_sample_Level2_logmean.h5ad")
    min_cell_thresh = 5
    expected_genes = 21002

    print("=" * 70)
    print("[STEP1] load retain gene list")
    keep_genes = load_gene_list(gene_list_fp)
    print(f"number of retain genes: {len(keep_genes)}")
    assert len(set(keep_genes)) == expected_genes, \
        f"retainGenes count mismatch: {len(set(keep_genes))} != {expected_genes}"

    print("=" * 70)
    print("[STEP2] load log-normalized adata, subset to retain genes")
    adata = sc.read_h5ad(raw_h5ad)
    adata = adata[:, adata.var_names.isin(keep_genes)].copy()

    # ensure grouping cols are str; dc.pp.pseudobulk builds index as
    # f"{sample}_{group}" and crashes on numpy.int32 / int64 sample IDs.
    adata.obs["samples"] = adata.obs["samples"].astype(str)
    adata.obs["Level2"] = adata.obs["Level2"].astype(str)

    n_cells, n_genes = adata.shape
    n_sample = adata.obs["samples"].nunique()
    n_level2 = adata.obs["Level2"].nunique()
    print(f"raw adata shape after gene subset: cells={n_cells}, genes={n_genes}")
    print(f"unique samples: {n_sample}")
    print(f"unique Level2 celltypes: {n_level2}")

    print("=" * 70)
    print("[STEP3] dc.pp.pseudobulk: log data -> group mean")
    pbulk_adata = dc.pp.pseudobulk(
        adata=adata,
        sample_col="samples",
        groups_col="Level2",
        mode="mean",
    )
    pb_n_obs, pb_n_vars = pbulk_adata.shape
    print(f"pseudobulk shape before filtering: n_obs={pb_n_obs}, n_genes={pb_n_vars}")
    print(f"pseudobulk obs columns: {list(pbulk_adata.obs.columns)}")
    print(f"pseudobulk unique samples: {pbulk_adata.obs['samples'].nunique()}")
    print(f"pseudobulk unique Level2: {pbulk_adata.obs['Level2'].nunique()}")

    # compute cell_count from original adata, keyed by obs_names string
    # ("sample_group") which is exactly what decoupler builds internally.
    # This avoids pandas MultiIndex .loc type-mismatch issues.
    print("=" * 70)
    print("[STEP3.5] compute cell_count from original adata (keyed by obs_names)")
    cell_count_series = adata.obs.groupby(["samples", "Level2"]).size()
    cell_count_dict = {f"{s}_{g}": int(c) for (s, g), c in cell_count_series.items()}
    pbulk_adata.obs["cell_count"] = [
        cell_count_dict.get(idx, 0) for idx in pbulk_adata.obs_names
    ]
    n_missing = int((pbulk_adata.obs["cell_count"] == 0).sum())
    print(f"cell_count range: [{pbulk_adata.obs['cell_count'].min()}, {pbulk_adata.obs['cell_count'].max()}]")
    print(f"groups with missing count (set to 0): {n_missing}")

    # filter groups with cell_count >= min_cell_thresh
    pbulk_adata = pbulk_adata[pbulk_adata.obs["cell_count"] >= min_cell_thresh, :].copy()
    pb_filt_obs, pb_filt_vars = pbulk_adata.shape
    print(f"filter threshold cell_count >= {min_cell_thresh}")
    print(f"pseudobulk shape after filtering: n_obs={pb_filt_obs}, n_genes={pb_filt_vars}")
    print(f"filtered pseudobulk unique samples: {pbulk_adata.obs['samples'].nunique()}")
    print(f"filtered pseudobulk unique Level2: {pbulk_adata.obs['Level2'].nunique()}")

    print("=" * 70)
    print("[STEP4] per-gene z-score (zero mean, unit variance across pseudobulks)")
    # NM paper: "log transformed and scaled the expression values to zero mean
    # and unit variance, to reduce the impact of highly expressed genes."
    # Input is already log-normalized, so apply z-score directly.
    # Scale by gene (axis=0): each gene normalized across all pseudobulk samples.
    X = pbulk_adata.X.toarray() if hasattr(pbulk_adata.X, "toarray") else np.asarray(pbulk_adata.X)
    gene_mean = X.mean(axis=0)
    gene_std = X.std(axis=0, ddof=0)
    gene_std_safe = np.where(gene_std == 0, 1.0, gene_std)
    X_scaled = (X - gene_mean) / gene_std_safe
    pbulk_adata.X = X_scaled
    n_zero_std = int((gene_std == 0).sum())
    print(f"z-score done. zero-std genes (set to 0): {n_zero_std}")
    print(f"scaled matrix range: [{X_scaled.min():.3f}, {X_scaled.max():.3f}]")

    print("=" * 70)
    print("[STEP5] save main h5ad output")
    pbulk_adata.write_h5ad(pb_h5ad)
    print(f"pseudobulk h5ad saved: {pb_h5ad}")

    print("=" * 70)
    print("[STEP6] final summary of saved pseudobulk")
    final_n_obs, final_n_vars = pbulk_adata.shape
    final_samples = sorted(pbulk_adata.obs["samples"].unique())
    final_level2 = sorted(pbulk_adata.obs["Level2"].unique())
    print(f"pseudobulk total genes (n_vars): {final_n_vars}")
    print(f"pseudobulk total profiles (n_obs):   {final_n_obs}")
    print(f"number of unique samples:           {len(final_samples)}")
    print(f"number of unique Level2 celltypes:  {len(final_level2)}")
    print("--- sample list ---")
    print(", ".join(final_samples))
    print("--- Level2 list ---")
    print(", ".join(final_level2))
    print("--- per-sample Level2 counts ---")
    sample_l2_counts = pbulk_adata.obs.groupby("samples")["Level2"].nunique()
    for s in final_samples:
        print(f"  {s}: {sample_l2_counts[s]} Level2 subtypes")
    print("log-mean pseudobulk aggregation + z-score finished")


if __name__ == "__main__":
    main()

