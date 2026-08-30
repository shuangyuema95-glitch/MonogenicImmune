# -*- coding: utf-8 -*-
"""
Second-round ULM: pathway activity quantification
data: pseudobulk expression (pseudobulk x gene, z-scored)  [response]
net: annotated factor-gene weights (factor x gene, continuous) [explanatory]
tval=True -> return t-statistic
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
from statsmodels.stats.multitest import multipletests

WORK_DIR = "/mnt/data/work/masy/scAltas/2026scRNAcohort"
PSEUDOBULK_H5AD = os.path.join(WORK_DIR, "pseudobulk_sample_Level2_logmean.h5ad")
FACTOR_CSV = os.path.join(WORK_DIR, "lam0001_dict242_218factors_labeled.csv")
ULM_OUT = os.path.join(WORK_DIR, "lam0001_dict242_second_ulm_result.csv")
TSTAT_WIDE_OUT = os.path.join(WORK_DIR, "lam0001_dict242_second_ulm_tstat_wide.csv")

def main():
    print("=" * 70)
    print("[STEP1] Load pseudobulk expression (data / response)")
    pbulk = sc.read_h5ad(PSEUDOBULK_H5AD)
    print(f"  pseudobulk shape: {pbulk.shape}")
    print(f"  obs columns: {list(pbulk.obs.columns)}")
    print(f"  first 3 obs_names: {pbulk.obs_names[:3].tolist()}")

    X = pbulk.X.toarray() if hasattr(pbulk.X, "toarray") else np.asarray(pbulk.X)
    df_data = pd.DataFrame(X, index=pbulk.obs_names, columns=pbulk.var_names)
    print(f"  df_data shape (pseudobulk x gene): {df_data.shape}")

    print("=" * 70)
    print("[STEP2] Load 218 factors, build long-format net (explanatory)")
    df_factor_raw = pd.read_csv(FACTOR_CSV)
    print(f"  raw csv shape: {df_factor_raw.shape}")
    print(f"  raw columns first 5: {df_factor_raw.columns[:5].tolist()}")

    # auto-detect label column: string column containing 'Factor_'
    label_col = None
    for col in df_factor_raw.columns:
        if df_factor_raw[col].dtype == object:
            sample = df_factor_raw[col].dropna().astype(str).head(30)
            if sample.str.contains('Factor_', na=False).any():
                label_col = col
                break
    if label_col is None:
        label_col = df_factor_raw.columns[0]
        print(f"  [WARN] no Factor_ column found, using first column: {label_col}")

    print(f"  label column detected: '{label_col}'")
    df_factor = df_factor_raw.set_index(label_col)
    # keep only numeric columns (drops any stray row-number columns like 'Unnamed: 0')
    df_factor = df_factor.select_dtypes(include=[np.number])
    print(f"  factor matrix shape (factor x gene): {df_factor.shape}")  # expect (218, 21002)
    print(f"  first 3 factor labels: {df_factor.index[:3].tolist()}")
    print(f"  first 3 genes: {df_factor.columns[:3].tolist()}")
    print(f"  any NA: {df_factor.isna().any().any()}")

    # wide (factor x gene) -> long (source, target, weight)
    df_net = df_factor.stack().reset_index()
    df_net.columns = ["source", "target", "weight"]
    print(f"  net long shape: {df_net.shape}")
    print(f"  unique factors (source): {df_net['source'].nunique()}")
    print(f"  unique genes (target): {df_net['target'].nunique()}")
    print(f"  weight range: [{df_net['weight'].min():.4f}, {df_net['weight'].max():.4f}]")

    # filter net genes to match data columns
    valid_genes = set(df_data.columns)
    df_net = df_net[df_net["target"].isin(valid_genes)].copy()
    print(f"  net after gene filter: {df_net.shape}")

    print("=" * 70)
    print("[STEP3] Run second-round ULM (tval=True -> t-statistic)")
    scores_wide, pvals_wide = dc.mt.ulm(
        data=df_data,
        net=df_net,
        tmin=3,
        tval=True,
        verbose=True
    )
    print(f"  scores_wide (pseudobulk x factor): {scores_wide.shape}")
    print(f"  pvals_wide (pseudobulk x factor): {pvals_wide.shape}")

    print("=" * 70)
    print("[STEP4] Long format + BH-FDR + metadata")
    scores_long = scores_wide.stack().reset_index()
    scores_long.columns = ["pseudobulk", "factor", "t_stat"]
    pvals_long = pvals_wide.stack().reset_index()
    pvals_long.columns = ["pseudobulk", "factor", "p_value"]
    res_df = pd.merge(scores_long, pvals_long, on=["pseudobulk", "factor"])

    _, res_df["fdr"], _, _ = multipletests(res_df["p_value"].values, method="fdr_bh")

    # split pseudobulk name into sample + Level2
    # obs_names like '108PLD4_B01_Naive_BACH2' -> first underscore splits sample vs rest
    # but sample may contain underscores; use obs metadata instead
    res_df = res_df.merge(
        pbulk.obs[["samples", "Level2"]].reset_index().rename(columns={"index": "pseudobulk"}),
        on="pseudobulk", how="left"
    )

    print(f"  total pairs: {len(res_df)}")
    print(f"  unique pseudobulks: {res_df['pseudobulk'].nunique()}")
    print(f"  unique factors: {res_df['factor'].nunique()}")
    print(f"  significant (fdr<0.05): {(res_df['fdr'] < 0.05).sum()}")

    print("=" * 70)
    print("[STEP5] Save outputs")
    res_df.to_csv(ULM_OUT, index=False)
    print(f"  long result saved: {ULM_OUT}")
    scores_wide.to_csv(TSTAT_WIDE_OUT)
    print(f"  t-stat wide matrix saved: {TSTAT_WIDE_OUT}")
    print("second-round ULM finished")

if __name__ == "__main__":
    main()