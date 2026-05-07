import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd
import pickle as pkl
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

input="/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/composition"
os.chdir(input)
pcaNumber = 50

adata_merged = sc.read_h5ad(f"{pcaNumber}_final_merged_object.h5ad")
#adata_merged.obs['cell_type']=adata_merged.obs['leiden_harmony']
meta = pd.read_csv("E:\\AID cohort\\code\\metainfo.txt", sep="\t")

samples_series = adata_merged.obs['samples']
genes = samples_series.str.replace(r'^\d+', '', regex=True)
adata_merged.obs['gene'] = genes

type_gene_dict = {"NF-κB": ["TNFAIP3", "RELA", "NOD2"],
    "Uncategorized": ["ADA2", "ELF4", "CSF3R"],
    "Osteoclast function": ["OGFRL1", "ACP5"],
    "Ca2+ flux-PLC axis": ["PLCG1", "PLCG2"],
    "Endolysosomal nucleic acid sensing": ["TLR7", "TLR8", "TRAF3", "PLD4"],
    "Arachidonic acid metabolism": ["TBXAS1"],
    "Inflammasome/IL1β axis": ["IL1R1", "NLRC4", "PSTPIP1", "CDC42", "LPIN2"],
    "Immune metabolic regulation axis": ["LACC1", "SLC7A7"],
    "Inbore errors of cell death": ["RIPK1", "OTULIN", "RNF31"],
    "Negative regulation of type I interferon signaling": ["ISG15", "USP18", "SOCS1"],
    "Cytoskeleton and small GTPase signaling": ["GNAI2", "KRAS"],
    "TBK1-IRF3 axis": ["STING", "TREX1", "DDX58", "IFIH1"],
    "Protein homeostasis axis": ["UBA1", "PSMD12"]}

type_gene_list = [(t, g) for t, genes_list in type_gene_dict.items() for g in genes_list]
type_gene_df = pd.DataFrame(type_gene_list, columns=['type', 'gene'])
adata_merged.obs = adata_merged.obs.merge(type_gene_df, on='gene', how='left')

meta['dataset'] = meta['dataset'].str.strip()
meta['status'] = meta['status'].str.strip()
meta['sex'] = meta['sex'].str.strip()
meta['age'] = meta['age'].str.strip()
obs = adata_merged.obs.copy()
obs = obs.merge(meta[['dataset','status','sex','age']],left_on='samples', right_on='dataset',how='left')
adata_merged.obs = obs
print(adata_merged.obs[['samples','gene','status','type','sex','age']].head(20))



def generate_cell_count_matrices(adata, cell_type_col="cell_type"):
    """
    Generate 4 cell count matrices and save to CSV.
    Only 'samples' includes covariates: age_clean, sex, status, age_group.
    """
    obs = adata.obs.copy()

    # Count function
    def get_counts(group_col):
        return pd.crosstab(obs[group_col], obs[cell_type_col])

    # 1. Counts by sample
    df_samples = get_counts("samples")

    # Clean age
    obs["age_val"] = obs["age"].astype(str).str.replace(r"[a-zA-Z]", "", regex=True)
    obs["age_val"] = pd.to_numeric(obs["age_val"], errors="coerce")

    # Covariates (unique per sample)
    cov = obs[["samples", "age_val", "sex", "status"]].drop_duplicates("samples")

    # Safe merge (NO INDEX JOIN ERROR)
    df_samples = df_samples.reset_index().merge(cov, on="samples", how="left").set_index("samples")

    # Age binning
    bins = [0, 16, 30, 40, 50, 60, 70, np.inf]
    labels = ["1-16", "17-30", "32-40", "40-50", "51-60", "61-70", ">70"]
    df_samples["age_group"] = pd.cut(df_samples["age_val"], bins=bins, labels=labels, right=False)

    # 2. Other counts
    df_gene = get_counts("gene")
    df_status = get_counts("status")
    df_type = get_counts("type")

    # Save
    df_samples.to_csv("cell_counts_by_samples.csv")
    df_gene.to_csv("cell_counts_by_gene.csv")
    df_status.to_csv("cell_counts_by_status.csv")
    df_type.to_csv("cell_counts_by_type.csv")
generate_cell_count_matrices(adata_merged)




def run_scCODA_batch_and_export(
    csv_path="haber_counts.csv",
    control_group="Control",
    reference_cell_type="Goblet",
    fdr=0.4,
    sampler="hmc",
    output_csv="scCODA_results.csv"
):
    """Full pipeline: read CSV cell counts → run scCODA for all groups vs control → export plot table. """
    cell_counts = pd.read_csv(csv_path)
    data_all = dat.from_pandas(cell_counts, covariate_columns=["Mouse"])
    data_all.obs["Condition"] = data_all.obs["Mouse"].str.replace(r"_[0-9]+", "", regex=True)

    groups = [g for g in data_all.obs["Condition"].unique() if g != control_group]
    all_res = []

    for group in groups:
        print(f"\n=== Running: {group} vs {control_group} | Sampler: {sampler} ===")
        sub_data = data_all[data_all.obs["Condition"].isin([control_group, group])]

        model = mod.CompositionalAnalysis(
            sub_data,
            formula="Condition",
            reference_cell_type=reference_cell_type
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

# Run HMC (fast for plotting)
result_table_hmc = run_scCODA_batch_and_export(
    csv_path="haber_counts.csv",
    control_group="Control",
    reference_cell_type="Goblet",
    fdr=0.4,
    sampler="hmc",
    output_csv="scCODA_HMC_for_R.csv"
)

# Run NUTS (for publication)
result_table_nuts = run_scCODA_batch_and_export(
    csv_path="haber_counts.csv",
    control_group="Control",
    reference_cell_type="Goblet",
    fdr=0.4,
    sampler="nuts",
    output_csv="scCODA_NUTS_for_R.csv"
)