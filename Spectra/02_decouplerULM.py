import os
import pickle
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

work_dir = "/mnt/data/work/masy/scAltas/2026scRNAcohort"
pkl_path = os.path.join(work_dir, "gene_dict.pkl")
factor_csv = os.path.join(work_dir, "spectra_run_factors.csv")
binary_out = os.path.join(work_dir, "binary_gene_pathway_matrix.csv")
out_csv = os.path.join(work_dir, "ulm_factor_pathway_result.csv")

print("="*60)
print("[STEP 1] Load gene_dict from pkl")
with open(pkl_path, "rb") as f:
    gene_dict = pickle.load(f)

# build binary matrix directly in memory
rows = []
for celltype, pw_dict in gene_dict.items():
    for pw_name, gene_list in pw_dict.items():
        for g in gene_list:
            rows.append({"source": pw_name, "target": g})
df_long = pd.DataFrame(rows)
sources = sorted(df_long["source"].unique())
targets = sorted(df_long["target"].unique())
binary = pd.DataFrame(0, index=targets, columns=sources, dtype=int)
for _, row in df_long.iterrows():
    binary.loc[row["target"], row["source"]] = 1

print(f"binary matrix shape(gene × pathway): {binary.shape}")
binary.to_csv(binary_out)
print(f"save binary matrix to {binary_out}")

print("="*60)
print("[STEP 2] Load spectra factor matrix and align genes")
df_factors = pd.read_csv(factor_csv, index_col=0)
Y = df_factors.T.copy()
common_genes = Y.index.intersection(binary.index)
print(f"common genes count: {len(common_genes)}")
Y = Y.loc[common_genes]
X = binary.loc[common_genes]
print(f"Y(gene × factor): {Y.shape}")
print(f"X(gene × pathway): {X.shape}")

print("="*60)
print("[STEP 3] Run OLS‑ULM regression, replicate Spectra paper")
pathway_list = X.columns.tolist()
factor_list = Y.columns.tolist()
res_records = []
for pw in pathway_list:
    x_vec = X[pw].values
    n_gene = int(x_vec.sum())
    if n_gene < 3:
        print(f"skip {pw}, gene count {n_gene} <3")
        continue
    Xm = sm.add_constant(x_vec)
    for fac in factor_list:
        y_vec = Y[fac].values
        fit = sm.OLS(y_vec, Xm).fit()
        res_records.append({
            "gene_set": pw,
            "factor_id": fac,
            "estimate": fit.params[1],
            "p_value": fit.pvalues[1]
        })

res_df = pd.DataFrame(res_records)
_, fdr, _, _ = multipletests(res_df["p_value"].values, method="fdr_bh")
res_df["fdr"] = fdr
res_df.to_csv(out_csv, index=False)

print("="*60)
print("[STEP 4] Finished")
print(f"result saved: {out_csv}")
print(f"total pairs: {res_df.shape[0]}")
mask = (res_df["fdr"] < 0.05) & (res_df["estimate"] > 0)
print(f"significant fdr<0.05 & estimate>0: {mask.sum()}")
