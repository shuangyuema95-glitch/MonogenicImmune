import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

print("="*60)
print("[Step1] Start reading h5ad data")
print("="*60)
adata = sc.read_h5ad("716_pbmc1.h5ad")
print("[Step1] Extract raw count matrix via adata.raw.to_adata()")
raw_adata = adata.raw.to_adata()
# Fix raw var_names bug
raw_adata.var_names = adata.var_names.copy()
del adata
print("[Step1] Original adata removed to save memory")

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
# =========================================================

print("\n"+"="*60)
print("[Step2] Assign condition metadata")
print("="*60)
obs = raw_adata.obs.copy()
obs["condition"] = obs.apply(assign_condition, axis=1)
raw_adata.obs["condition"] = obs["condition"]
raw_adata = raw_adata[raw_adata.obs["condition"].notna()].copy()
print("[Step2] condition metadata constructed successfully")

print("\n"+"="*60)
print("[Step3] Run pseudobulk aggregation: only group by samples")
print("="*60)
pdata = dc.pp.pseudobulk(adata=raw_adata, sample_col="samples", groups_col=None, mode="sum")
print(f"[Step3] Pseudobulk matrix shape (samples x genes): {pdata.shape}")

check_df = pdata.obs[["samples","condition"]].drop_duplicates().sort_values("samples")
print("[Step3] Preview sample-condition mapping:")
print(check_df.head(40))
check_df.to_csv("CHECK_global_sample_condition_mapping.csv", index=False)
print("[Step3] Count per condition:")
print(pdata.obs["condition"].value_counts())

del raw_adata
print("[Step3] raw_adata deleted")

print("\n"+"="*60)
print("[Step4] Initialize PyDESeq2 model | design = ~ condition")
print("="*60)
inference = DefaultInference(n_cpus=8)

dds = DeseqDataSet(
    adata=pdata,
    design="~ condition",
    refit_cooks=True,
    inference=inference
)
dds.deseq2()
print("[Step4] DESeq2 dispersion estimation finished")

print("\n"+"="*60)
print("[Step5] Iterate all disease groups vs HC")
print("[Threshold] padj<0.05 & abs(log2FoldChange) >1")
print("="*60)
padj_cut = 0.01
lfc_cut  = 2

all_degs = []
all_conds = sorted(pdata.obs["condition"].unique())
test_conds = [c for c in all_conds if c != "HC"]

for cond in test_conds:
    target_count = (pdata.obs["condition"] == cond).sum()
    hc_count = (pdata.obs["condition"] == "HC").sum()
    if target_count < 2 or hc_count < 2:
        print(f"\n[Skip] {cond}_VS_HC: target={target_count}, HC={hc_count}, require >=2")
        continue

    print(f"\n[Compare] {cond} VS HC")
    stat_res = DeseqStats(dds, contrast=["condition", cond, "HC"], inference=inference)
    stat_res.summary()
    res_df = stat_res.results_df.copy()
    res_df["gene"] = res_df.index
    res_df["comparison"] = f"{cond}_VS_HC"
    sig_df = res_df[(res_df["padj"] < padj_cut) & (np.abs(res_df["log2FoldChange"]) > lfc_cut)].copy()
    all_degs.append(sig_df)
    print(f"[Compare] Significant DE genes count: {sig_df.shape[0]}")

print("\n"+"="*60)
print("[Step6] Merge all DE results and save csv")
print("="*60)
out_df = pd.concat(all_degs, axis=0).reset_index(drop=True)
out_df = out_df[["gene", "comparison", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
out_df.to_csv("DEG00515_pseudobulk_condition_VS_HC.csv", index=False)
print("[Step6] Output saved: DEG_global_pseudobulk_condition_VS_HC.csv")
print("[All Global Task Finished!]")
