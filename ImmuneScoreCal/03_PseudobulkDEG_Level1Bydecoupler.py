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

print("\n"+"="*60)
print("[Step2] Assign condition metadata")
print("="*60)
obs = raw_adata.obs.copy()
obs["condition"] = obs.apply(assign_condition, axis=1)
raw_adata.obs["condition"] = obs["condition"]
raw_adata = raw_adata[raw_adata.obs["condition"].notna()].copy()
print("[Step2] condition metadata constructed successfully")

print("\n"+"="*60)
print("[Step3] Run pseudobulk aggregation: samples + Level1")
print("="*60)
pdata = dc.pp.pseudobulk(adata=raw_adata, sample_col="samples", groups_col="Level1", mode="sum")
print(f"[Step3] Pseudobulk matrix shape: {pdata.shape}")

check_df = pdata.obs[["samples","Level1","condition"]].drop_duplicates().sort_values(["samples","Level1"])
check_df.to_csv("CHECK_Level1_sample_L1_condition_mapping.csv", index=False)

del raw_adata
print("[Step3] raw_adata deleted")

print("\n"+"="*60)
print("[Step4] Iterate each Level1 cell type for DE analysis")
print("[Threshold] padj<0.01 & abs(log2FoldChange) >2")
print("="*60)
padj_cut = 0.001
lfc_cut  = 2.5
inference = DefaultInference(n_cpus=8)
all_results = []
celltype_list = sorted(pdata.obs["Level1"].unique())

for ct in celltype_list:
    print(f"\n{'='*50}")
    print(f"[CellType] Process: {ct}")
    sub_pb = pdata[pdata.obs["Level1"] == ct].copy()
    cond_counts = sub_pb.obs["condition"].value_counts()
    print(f"Condition count inside {ct}:")
    print(cond_counts)
    
    if "HC" not in cond_counts.index:
        print(f"[SKIP] {ct} : missing HC control group")
        continue
    
    test_conds = [c for c in cond_counts.index if c != "HC"]
    try:
        dds = DeseqDataSet(
            adata=sub_pb,
            design="~ condition",
            refit_cooks=True,
            inference=inference
        )
        dds.deseq2()
    except Exception as err:
        print(f"[ERROR] {ct} model failed: {err}")
        continue

    for cond in test_conds:
        print(f"[Compare] {ct} | {cond} VS HC")
        stat_res = DeseqStats(dds, contrast=["condition", cond, "HC"], inference=inference)
        stat_res.summary()
        res_df = stat_res.results_df.copy()
        res_df["gene"] = res_df.index
        res_df["celltype"] = ct
        res_df["comparison"] = f"{cond}_VS_HC"
        sig_df = res_df[(res_df["padj"] < padj_cut) & (np.abs(res_df["log2FoldChange"]) > lfc_cut)].copy()
        all_results.append(sig_df)
        print(f"[Compare] DE gene count: {sig_df.shape[0]}")

print("\n"+"="*60)
print("[Step5] Merge all celltype DE results and save csv")
print("="*60)
final_df = pd.concat(all_results, axis=0).reset_index(drop=True)
final_df = final_df[["gene", "celltype", "comparison", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
final_df.to_csv("DEG_Level1_pseudobulk_condition_VS_HC.csv", index=False)
print("[Step5] Output saved: DEG_Level1_pseudobulk_condition_VS_HC.csv")
print("[All Level1 Task Finished!]")
