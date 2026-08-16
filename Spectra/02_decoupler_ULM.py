# -*- coding: utf-8 -*-
"""
Spectra factor level ULM analysis using decoupler.mt.ulm
data: factor × gene matrix (NO transpose)
net: long format dataframe with source, target, weight
fix: remove duplicate source‑target pairs, handle wide matrix output, manual BH‑FDR
"""
import os
import pickle
import pandas as pd
import decoupler as dc
from statsmodels.stats.multitest import multipletests


def main():
    work_dir = "/mnt/data/work/masy/scAltas/2026scRNAcohort"

    pkl_path = os.path.join(work_dir, "gene_dict.pkl")
    factor_csv = os.path.join(work_dir, "spectra_run_factors.csv")
    net_out = os.path.join(work_dir, "decoupler_net_weight.csv")
    binary_out = os.path.join(work_dir, "binary_gene_pathway_matrix.csv")
    ulm_out = os.path.join(work_dir, "mt_ulm_factor_pathway_result.csv")

    print("=" * 70)
    print("[STEP1] Build long format net from gene_dict")
    with open(pkl_path, "rb") as f:
        gene_dict = pickle.load(f)

    net_rows = []
    for cell_type, pathway_dict in gene_dict.items():
        for pathway_name, gene_list in pathway_dict.items():
            for gene in gene_list:
                net_rows.append({"source": pathway_name, "target": gene, "weight": 1.0})
    df_net = pd.DataFrame(net_rows)

    # drop duplicate source target pairs for decoupler validation
    print(f"before deduplicate net shape: {df_net.shape}")
    df_net = df_net.drop_duplicates(subset=["source", "target"], keep="first").copy()
    print(f"after deduplicate net shape: {df_net.shape}")

    df_net.to_csv(net_out, index=False)
    print(f"net file saved: {net_out}")

    # binary gene pathway matrix for archive, not used in computation
    sources = sorted(df_net["source"].unique())
    targets = sorted(df_net["target"].unique())
    binary_mat = pd.DataFrame(0, index=targets, columns=sources, dtype=int)
    for _, row in df_net.iterrows():
        binary_mat.loc[row["target"], row["source"]] = row["weight"]
    binary_mat.to_csv(binary_out)
    print(f"binary matrix shape(gene × pathway): {binary_mat.shape}")

    print("=" * 70)
    print("[STEP2] Load spectra factor matrix, NO TRANSPOSE")
    df_data = pd.read_csv(factor_csv, index_col=0)
    print(f"data matrix (Factor × Gene): {df_data.shape}")
    print(f"first 5 gene columns: {df_data.columns[:5].tolist()}")

    valid_genes = set(df_data.columns)
    df_net_filter = df_net[df_net["target"].isin(valid_genes)].copy()
    print(f"filtered net rows (match genes in data): {df_net_filter.shape[0]}")

    print("=" * 70)
    print("[STEP3] run dc.mt.ulm official api")
    ret = dc.mt.ulm(
        data=df_data,
        net=df_net_filter,
        tmin=3,
        tval=False,
        verbose=True
    )
    scores_wide, pvals_wide = ret
    print(f"scores_wide shape(factor × pathway): {scores_wide.shape}")
    print(f"pvals_wide shape(factor × pathway): {pvals_wide.shape}")

    # convert wide matrix to long format
    scores_long = scores_wide.stack().reset_index()
    scores_long.columns = ["factor_id", "gene_set", "estimate"]

    pvals_long = pvals_wide.stack().reset_index()
    pvals_long.columns = ["factor_id", "gene_set", "p_value"]

    res_df = pd.merge(scores_long, pvals_long, on=["factor_id", "gene_set"])

    # manual BH FDR correction
    _, fdr_corr, _, _ = multipletests(res_df["p_value"].values, method="fdr_bh")
    res_df["fdr"] = fdr_corr

    print("=" * 70)
    print("[STEP4] save result and summary statistics")
    res_df.to_csv(ulm_out, index=False)
    print(f"result saved: {ulm_out}")
    print(f"total factor‑pathway pairs: {res_df.shape[0]}")
    print(f"unique pathways in output: {res_df['gene_set'].nunique()}")
    sig_mask = (res_df["fdr"] < 0.05) & (res_df["estimate"] > 0)
    print(f"significant positive pairs (fdr<0.05): {sig_mask.sum()}")
    print("mt.ulm finished")


if __name__ == "__main__":
    main()

