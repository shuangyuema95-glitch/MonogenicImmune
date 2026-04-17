import os
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import time
import shutil
import scipy.sparse as sp
from scipy.io import mmwrite
import anndata as ad
import gzip

input = "E:\\AID cohort\\code\\result\\"
samples = ["137ACP5", "84ADA2", "159ADA2", "163COPA", "88GNAI2", "78IL1R1", "79IL1R1_2", "156ISG15", "157LACC1", "90OGFRL1", "36PLCG1", "37PLCG2", "108PLD4", "109PLD4", "110PLD4", "155PSMD12", "138RIPK1", "139RIPK1", "140RIPK1", "141RIPK1", "152RNF31", "149SOCS1", "162STING", "69TBXAS1", "70TBXAS1", "89TLR7", "82TNFAIP3", "83TNFAIP3", "154TRAF3", "151USP18", "134ELF4", "115KRAS", "81NLRC4", "91TLR8", "148UBA1", "1IRAK2", "3IRAK2", "165TREX1", "168TREX1", "169TREX1", "175IFIH1", "176IFIH1", "177IFIH1", "178ISG15", "179STAT4", "180SLC7A7", "181NOD2", "12HC", "17HC", "22HC", "75HC", "76HC", "77HC", "158HC", "160HC", "161HC", "9HC", "HC_GSE139324_GSM4138162", "HC_GSE139324_GSM4138163", "HC_GSE139324_GSM4138164", "HC_GSE139324_GSM4138165", "HC_GSE139324_GSM4138166", "HC_GSE139324_GSM4138167", "HC_10XGenomics10K", "HC_GSE156989_GSM4749756", "HC_GSE156989_GSM4749762", "HC_GSE156989_GSM4749768", "HC_GSE199445_GSM5973143", "HC_GSE199445_GSM5973144", "HC_GSE199445_GSM5973145", "HC_GSE199445_GSM5973146", "HC_GSE168732_GSM5160432", "HC_GSE168732_GSM5160434", "HC_GSE168732_GSM5160435", "SLE_GSE156989_GSM4749774", "SLE_GSE156989_GSM4749779", "SLE_GSE156989_GSM4749784", "SLE_GSE142016_GSM4217718", "SLE_GSE142016_GSM4217719", "SLE_GSE142016_GSM4217720", "SLE_GSE263931_GSM8207595", "SLE_GSE263931_GSM8207597", "SLE_GSE263931_GSM8207599", "SLE_GSE263931_GSM8207601", "SLE_GSE263931_GSM8207603", "SLE_GSE263931_GSM8207605", "SLE_GSE263931_GSM8207607", "SLE_GSE224198_GSM7017326", "SLE_GSE224198_GSM7017329", "SLE_GSE224198_GSM7017331", "SLE_GSE224198_GSM7017334", "JIA_GSE205095_GSM6205132", "JIA_GSE205095_GSM6205136", "JIA_GSE205095_GSM6205138", "JIA_GSE205095_GSM6205140", "JIA_GSE205095_GSM6205142", "JIA_GSE205095_GSM6205144", "KD_GSE168732_GSM5160417", "KD_GSE168732_GSM5160420", "KD_GSE168732_GSM5160422", "KD_GSE168732_GSM5160424", "KD_GSE168732_GSM5160427", "KD_GSE168732_GSM5160430", "STING_GSE226598_GSM7079986", "STING_GSE226598_GSM7079989", "STING_GSE226598_GSM7079991", "PLCG2_Figshare", "OTULIN_GSE199445_GSM5973147", "OTULIN_GSE199445_GSM5973148", "OTULIN_GSE199445_GSM5973149"]
outdir= input
def batch_filter_by_global_doublet(input, outdir, samples):
    start_time = time.time()

    adatas_list = []
    all_scores = []

    # Load data
    print("--- loading h5ad files ---")
    for sample in samples:
        file_path = os.path.join(input, f"{sample}.h5ad")
        if os.path.exists(file_path):
            ad = sc.read_h5ad(file_path)
            ad.obs["samples"] = sample
            adatas_list.append(ad)
            all_scores.extend(ad.obs["doublet_score"].tolist())
        else:
            print(f"missing: {file_path}")

    # Global doublet filtering
    global_threshold = np.nanquantile(all_scores, 0.95)
    print(f"doublet threshold: {global_threshold:.4f}")

    filtered_adatas = []
    stats_list = []

    for adata in adatas_list:

        sample = adata.obs["samples"].iloc[0]

        n_before = adata.n_obs

        # Filtering the doublets
        adata_filtered = adata[adata.obs["doublet_score"] <= global_threshold].copy()

        n_after = adata_filtered.n_obs
        n_removed = n_before - n_after
        removed_pct = (n_removed / n_before * 100) if n_before > 0 else 0

        # save filtered h5ad
        out_path = os.path.join(outdir, f"{sample}Redl.h5ad")
        adata_filtered.write(out_path, compression='gzip')

        # save the processed data to 10X matrix (normalized X)
        # Use adata_filtered directly
        adata_filtered.var_names_make_unique()
        X = adata_filtered.X

        if not sp.issparse(X):
            X = sp.csr_matrix(X)

        X = X.T

        sample_outdir = os.path.join(outdir, f"filter{sample}")
        os.makedirs(sample_outdir, exist_ok=True)

        obs_path = os.path.join(sample_outdir, "cell_metadata.xlsx")
        adata_filtered.obs.to_excel(obs_path, index=True)

        var_path = os.path.join(sample_outdir, "gene_metadata.xlsx")
        adata_filtered.var.to_excel(var_path, index=True)

        mtx_path = os.path.join(sample_outdir, "matrix.mtx")
        mmwrite(mtx_path, X)

        with open(mtx_path, 'rb') as f_in, gzip.open(mtx_path + ".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(mtx_path)

        barcodes_path = os.path.join(sample_outdir, "barcodes.tsv.gz")
        with gzip.open(barcodes_path, 'wt') as f:
            f.write("\n".join(adata_filtered.obs_names) + "\n")

        features_path = os.path.join(sample_outdir, "features.tsv.gz")
        with gzip.open(features_path, 'wt') as f:
            for gene in adata_filtered.var_names:
                f.write(f"{gene}\t{gene}\n")

        # save raw counts matrix
        X = adata_filtered.layers["counts"]

        if not sp.issparse(X):
            X = sp.csr_matrix(X)

        X = X.T

        sample_outdir = os.path.join(outdir, f"Counts{sample}")
        os.makedirs(sample_outdir, exist_ok=True)

        mtx_path = os.path.join(sample_outdir, "matrix.mtx")
        mmwrite(mtx_path, X)

        with open(mtx_path, 'rb') as f_in, gzip.open(mtx_path + ".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(mtx_path)

        barcodes_path = os.path.join(sample_outdir, "barcodes.tsv.gz")
        with gzip.open(barcodes_path, 'wt') as f:
            f.write("\n".join(adata_filtered.obs_names) + "\n")

        features_path = os.path.join(sample_outdir, "features.tsv.gz")
        with gzip.open(features_path, 'wt') as f:
            for gene in adata_filtered.var_names:
                f.write(f"{gene}\t{gene}\n")


        stats_list.append({
            "Sample": sample,
            "Cells_Before": n_before,
            "Doublets_Removed": n_removed,
            "Cells_After": n_after,
            "Removal_Rate_Percent": round(removed_pct, 2),
            "Global_Threshold_Used": round(global_threshold, 4)
        })
        print(f"\n--- Removing doublet and saving matrix for: {sample} were done---")

        filtered_adatas.append(adata_filtered)

    # save stats
    stats_df = pd.DataFrame(stats_list)
    stats_excel_path = os.path.join(outdir, "Doublet_Filtering_Detailed_Stats.xlsx")
    stats_df.to_excel(stats_excel_path, index=False)

    combined_adata = anndata.concat(filtered_adatas, index_unique="-", join="outer")
    combined_adata.write(os.path.join(input, f"{len(samples)}.h5ad"), compression="gzip")

    print(f"Cells counting done in: {stats_excel_path}\n")
    print(f"The threshold of doublets was: {global_threshold}")
    end_time = time.time()
    print(f"\n===== Remove doublets for all samples done =====\nTotal runtime: {round((end_time - start_time) / 60, 2)} minutes")

batch_filter_by_global_doublet(input, outdir, samples)