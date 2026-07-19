import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd
import matplotlib.pyplot as plt
import sccoda
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization
import sccoda.datasets as scd
import arviz as az
import anndata as ad
import re
import numpy as np

input="/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/composition/Level2"
os.chdir(input)

############Find reference cell
def get_min_variance_ref_cell(df_count):
    count_mat = df_count.copy()
    group_col = "samples"
    cell_cols = [col for col in count_mat.columns if col != group_col]
    var_dict = {}
    for cell in cell_cols:
        group_mean = count_mat.groupby(group_col)[cell].mean()
        var_dict[cell] = np.var(group_mean)
    var_df = pd.DataFrame(list(var_dict.items()), columns=["cell_type", "group_variance"])
    var_df = var_df.sort_values("group_variance", ascending=True).reset_index(drop=True)
    ref_cell = var_df.iloc[0]["cell_type"]
    print(f"Auto selected reference cell (min group variance): {ref_cell}")
    print("Variance ranking (top5 least variable):")
    print(var_df.head(5))
    return ref_cell

#---result is: DC07_cDC2_ISG15

############scCODA running
def run_scCODA_batch_no_cov(
        csv_path,
        control_group="HC",
        fdr=0.4,
        sampler="nuts",
        output_csv=""
):
    cell_counts = pd.read_csv(csv_path)
    cell_counts.columns = cell_counts.columns.str.strip()
    print("Cleaned columns:", cell_counts.columns.tolist())

    ref_cell = get_min_variance_ref_cell(cell_counts)

    # 关键修复：covariate_columns 传入["samples"]，该列会自动存入data_all.obs
    data_all = dat.from_pandas(cell_counts, covariate_columns=["samples"])
    # 不需要去除后缀，直接复制samples作为Condition
    data_all.obs["Condition"] = data_all.obs["samples"]

    groups = [g for g in data_all.obs["Condition"].unique() if g != control_group]
    all_res = []

    for group in groups:
        print(f"\n=== Running: {group} vs {control_group} | Sampler: {sampler} ===")
        sub_data = data_all[data_all.obs["Condition"].isin([control_group, group])]

        model = mod.CompositionalAnalysis(
            sub_data,
            formula="Condition",
            reference_cell_type=ref_cell
        )

        if sampler.lower() == "nuts":
            sim = model.sample_nuts()
        else:
            sim = model.sample_hmc()

        sim.set_fdr(est_fdr=fdr)
        eff = sim.effect_df.reset_index()
        cred = sim.credible_effects().reset_index(name="Significant")
        df = pd.merge(eff, cred, on=["Covariate", "Cell Type"])
        df["Group"] = group
        all_res.append(df)

    final_df = pd.concat(all_res, ignore_index=True)
    final_df = final_df.rename(columns={"log2-fold change": "log2FC"})
    final_df.to_csv(output_csv, index=False)
    print(f"\nDone! Output file: {output_csv}")
    return final_df

res_condition = run_scCODA_batch_no_cov(
    csv_path="condition_count.csv",
    output_csv="scCODA_condition_out.csv"
)

