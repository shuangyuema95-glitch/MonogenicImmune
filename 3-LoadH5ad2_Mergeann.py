import os
import scanpy as sc
import pandas as pd
import anndata
import time
import pickle

input = "E:\\AID cohort\\code\\result\\"
samples = ["137ACP5", "84ADA2", "159ADA2", "163COPA", "88GNAI2", "78IL1R1", "79IL1R1_2", "156ISG15", "157LACC1", "90OGFRL1", "36PLCG1", "37PLCG2", "108PLD4", "109PLD4", "110PLD4", "155PSMD12", "138RIPK1", "139RIPK1", "140RIPK1", "141RIPK1", "152RNF31", "149SOCS1", "162STING", "69TBXAS1", "70TBXAS1", "89TLR7", "82TNFAIP3", "83TNFAIP3", "154TRAF3", "151USP18", "134ELF4", "115KRAS", "81NLRC4", "91TLR8", "148UBA1", "1IRAK2", "3IRAK2", "165TREX1", "168TREX1", "169TREX1", "175IFIH1", "176IFIH1", "177IFIH1", "178ISG15", "179STAT4", "180SLC7A7", "181NOD2", "12HC", "17HC", "22HC", "75HC", "76HC", "77HC", "158HC", "160HC", "161HC", "9HC", "HC_GSE139324_GSM4138162", "HC_GSE139324_GSM4138163", "HC_GSE139324_GSM4138164", "HC_GSE139324_GSM4138165", "HC_GSE139324_GSM4138166", "HC_GSE139324_GSM4138167", "HC_10XGenomics10K", "HC_GSE156989_GSM4749756", "HC_GSE156989_GSM4749762", "HC_GSE156989_GSM4749768", "HC_GSE199445_GSM5973143", "HC_GSE199445_GSM5973144", "HC_GSE199445_GSM5973145", "HC_GSE199445_GSM5973146", "HC_GSE168732_GSM5160432", "HC_GSE168732_GSM5160434", "HC_GSE168732_GSM5160435", "SLE_GSE156989_GSM4749774", "SLE_GSE156989_GSM4749779", "SLE_GSE156989_GSM4749784", "SLE_GSE142016_GSM4217718", "SLE_GSE142016_GSM4217719", "SLE_GSE142016_GSM4217720", "SLE_GSE263931_GSM8207595", "SLE_GSE263931_GSM8207597", "SLE_GSE263931_GSM8207599", "SLE_GSE263931_GSM8207601", "SLE_GSE263931_GSM8207603", "SLE_GSE263931_GSM8207605", "SLE_GSE263931_GSM8207607", "SLE_GSE224198_GSM7017326", "SLE_GSE224198_GSM7017329", "SLE_GSE224198_GSM7017331", "SLE_GSE224198_GSM7017334", "JIA_GSE205095_GSM6205132", "JIA_GSE205095_GSM6205136", "JIA_GSE205095_GSM6205138", "JIA_GSE205095_GSM6205140", "JIA_GSE205095_GSM6205142", "JIA_GSE205095_GSM6205144", "KD_GSE168732_GSM5160417", "KD_GSE168732_GSM5160420", "KD_GSE168732_GSM5160422", "KD_GSE168732_GSM5160424", "KD_GSE168732_GSM5160427", "KD_GSE168732_GSM5160430", "STING_GSE226598_GSM7079986", "STING_GSE226598_GSM7079989", "STING_GSE226598_GSM7079991", "PLCG2_Figshare", "OTULIN_GSE199445_GSM5973147", "OTULIN_GSE199445_GSM5973148", "OTULIN_GSE199445_GSM5973149"]

def batch_filter_by_global_doublet(input, samples, top_k=3000):
    start_time = time.time()

    filtered_adatas = []

    print("--- loading h5ad files ---")
    for sample in samples:
        file_path = os.path.join(input, f"{sample}Redl.h5ad")

        if os.path.exists(file_path):
            ad = sc.read_h5ad(file_path)
            ad.obs["samples"] = sample
            filtered_adatas.append(ad)
        else:
            print(f"missing: {file_path}")


    # HVGs selection by counting all the number of all genes in each dataset
    print("--- extracting HVGs per dataset ---")
    hvgs_list = []
    for ad in filtered_adatas:
        # Use the boolean mask for selecting "highly_variable" genes
        hvg_genes = ad.var_names[ad.var['highly_variable']].tolist()  # Only True ones
        hvgs_list.append(set(hvg_genes))

    # Save HVGs list
    hvgs_save_path = os.path.join(input, "hvgs_list.pkl")
    with open(hvgs_save_path, "wb") as f:
        pickle.dump(hvgs_list, f)

    print("--- building HVG frequency map ---")
    gene_count = {}
    for hvg_set in hvgs_list:
        for g in hvg_set:
            gene_count[g] = gene_count.get(g, 0) + 1

    sorted_genes = sorted(gene_count.items(), key=lambda x: x[1], reverse=True)
    # Cutoff rule: frequency at top_k boundary
    if len(sorted_genes) < top_k:
        cutoff = sorted_genes[-1][1]
    else:
        cutoff = sorted_genes[top_k - 1][1]

    individual_crossintersect = {g for g, c in gene_count.items() if c >= cutoff}

    print(f"HVG final: {len(individual_crossintersect)} genes")
    print(f"Cutoff frequency: {cutoff}")
    cross_save_path = os.path.join(input, "top3000rank_interhvgs.csv")
    pd.Series(sorted(list(individual_crossintersect))).to_csv(cross_save_path, index=False, header=False)

    # Merging datasets
    print("--- merging datasets ---")

    combined_adata = anndata.concat(filtered_adatas, index_unique="-", join="outer")
    combined_adata.write(os.path.join(input, f"{len(samples)}.h5ad"), compression="gzip")

    end_time = time.time()
    print(f"done in {(end_time - start_time) / 60:.2f} min")
    print("Data saved successfully!")

# Call the function
batch_filter_by_global_doublet(input, samples)

# Reload pkl file to check file format
with open('hvgs_list.pkl', 'rb') as f:
    hvgs_list = pickle.load(f)