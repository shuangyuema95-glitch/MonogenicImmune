import os
import shutil
import pandas as pd
import scanpy as sc
import numpy as np
import scrublet as scr
import scipy.sparse as sp
from scipy.io import mmwrite
import anndata as ad
import openpyxl
import time
import re
import gzip

from Cython.Compiler.Parsing import p_c_arg_list

# Set directory for input files
fold2 = "E:/AID cohort/code"
os.chdir(fold2)
outdir="E:/AID cohort/code/result"
samples=["137ACP5","84ADA2","159ADA2","163COPA","88GNAI2","78IL1R1","79IL1R1_2","156ISG15","157LACC1","90OGFRL1","36PLCG1","37PLCG2","108PLD4","109PLD4","110PLD4","155PSMD12","138RIPK1","139RIPK1","140RIPK1","141RIPK1","152RNF31","149SOCS1","162STING","69TBXAS1","70TBXAS1","89TLR7","82TNFAIP3","83TNFAIP3","154TRAF3","151USP18","134ELF4","115KRAS","81NLRC4","91TLR8","148UBA1","1IRAK2","3IRAK2","165TREX1","168TREX1","169TREX1","175IFIH1","176IFIH1","177IFIH1","178ISG15","179STAT4","180SLC7A7","181NOD2","12HC","17HC","22HC","75HC","76HC","77HC","158HC","160HC","161HC","9HC","HC_GSE139324_GSM4138162","HC_GSE139324_GSM4138163","HC_GSE139324_GSM4138164","HC_GSE139324_GSM4138165","HC_GSE139324_GSM4138166","HC_GSE139324_GSM4138167","HC_10XGenomics10K","HC_GSE156989_GSM4749756","HC_GSE156989_GSM4749762","HC_GSE156989_GSM4749768","HC_GSE199445_GSM5973143","HC_GSE199445_GSM5973144","HC_GSE199445_GSM5973145","HC_GSE199445_GSM5973146","HC_GSE168732_GSM5160432","HC_GSE168732_GSM5160434","HC_GSE168732_GSM5160435","HC_CELLxGENE","SLE_GSE156989_GSM4749774","SLE_GSE156989_GSM4749779","SLE_GSE156989_GSM4749784","SLE_GSE142016_GSM4217718","SLE_GSE142016_GSM4217719","SLE_GSE142016_GSM4217720","SLE_GSE263931_GSM8207595","SLE_GSE263931_GSM8207597","SLE_GSE263931_GSM8207599","SLE_GSE263931_GSM8207601","SLE_GSE263931_GSM8207603","SLE_GSE263931_GSM8207605","SLE_GSE263931_GSM8207607","SLE_GSE224198_GSM7017326","SLE_GSE224198_GSM7017329","SLE_GSE224198_GSM7017331","SLE_GSE224198_GSM7017334","JIA_GSE205095_GSM6205132","JIA_GSE205095_GSM6205136","JIA_GSE205095_GSM6205138","JIA_GSE205095_GSM6205140","JIA_GSE205095_GSM6205142","JIA_GSE205095_GSM6205144","KD_GSE168732_GSM5160417","KD_GSE168732_GSM5160420","KD_GSE168732_GSM5160422","KD_GSE168732_GSM5160424","KD_GSE168732_GSM5160427","KD_GSE168732_GSM5160430","STING_GSE226598_GSM7079986","STING_GSE226598_GSM7079989","STING_GSE226598_GSM7079991","PLCG2_Figshare","OTULIN_GSE199445_GSM5973147","OTULIN_GSE199445_GSM5973148","OTULIN_GSE199445_GSM5973149"]
metainfo=pd.read_csv("metainfo.txt", sep="\t",header=0)
# Define QC thresholds like in R
QCdata = pd.DataFrame({
    'metric': ["3' V3","3' V4","3' V2","5' V1","5' V2","Drop-seq"],
    'mt_threshold': [25,30,20,20,20,20],
    'rb_threshold': [60,60,65,60,60,60],
    'feature_threshold': [6000,8000,4000,6000,6000,4000],
    'count_threshold': [50000,60000,25000,50000,40000,25000]
})

# Quality control procedures
def prefilter_raw(samples,fold2, outdir):
    """
    Perform QC filtering, HVG selection, PCA, and doublet detection for multiple 10X scRNA-seq samples.
    """
    start_time = time.time()
    res_dict = {}

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for sample in samples:
        print(f"\n--- Processing sample: {sample} ---")

        # Load counts
        counts = sc.read_10x_mtx(os.path.join(fold2, sample), var_names='gene_symbols', cache=True)
        beforeQC_n = counts.n_obs
        beforeQC_g = counts.n_vars
        sc.pp.filter_genes(counts, min_cells=3)
        counts.obs['samples'] = sample

        # Basic Seurat-style QC metrics
        ribo_genes = [g for g in counts.var_names if re.match(r"^RP[SL]", g)]
        print(len(ribo_genes), "Ribosomal genes found")

        if len(ribo_genes) > 0:
            counts.obs['percent_ribo'] = (counts[:, ribo_genes].X.sum(axis=1) / counts.X.sum(axis=1)).A1 * 100
        else:
            counts.obs['percent_ribo'] = 0.0

        counts.obs['percent_mt'] = np.sum(counts[:, counts.var_names.str.startswith('MT-')].X, axis=1).A1 / np.sum(counts.X, axis=1).A1 * 100
        counts.obs['log10GenesPerUMI'] = np.log10(np.sum(counts.X > 0, axis=1).A1) / np.log10(np.sum(counts.X, axis=1).A1)#
        counts.obs['log10GenesPerUMI'] = (np.log10((counts.X>0).sum(axis=1).A1)/np.log10(counts.X.sum(axis=1).A1))
        counts.obs['n_counts'] = counts.X.sum(axis=1).A1
        counts.obs['n_genes'] = (counts.X > 0).sum(axis=1).A1

        # Match sample chemistry type
        chem_type = metainfo.loc[metainfo['dataset'] == sample, 'chemistry'].values[0]
        thres_row = QCdata.loc[QCdata['metric'] == chem_type]
        if thres_row.empty:
            raise ValueError(f"No QC threshold found for sample {sample} with chemistry {chem_type}")
        thres = thres_row.iloc[0]

        # Filter cells
        counts = counts[
            (np.array(counts.obs_vector('n_counts')) < thres['count_threshold']) &
            (np.array(counts.obs_vector('n_genes')) > 200) &
            (np.array(counts.obs_vector('n_genes')) < thres['feature_threshold']) &
            (counts.obs['percent_mt'] < thres['mt_threshold']) &
            (counts.obs['log10GenesPerUMI'] > 0.7)]

        counts.raw = counts.copy()  # Keep raw counts as input for doublets prediction

        exclude_genes = [g for g in counts.var_names if (g.startswith((
            'TRAV', 'TRAJ', 'TRAC',
            'TRBV', 'TRBJ', 'TRBC',
            'TRGV', 'TRGJ', 'TRGC',
            'TRDV', 'TRDJ', 'TRDC',
            'IGK', 'IGL', 'IGH',
            'RPS', 'RPL', 'MT-', 'HSP')))]

        # Exclude confounding genes
        keep_genes = [g for g in counts.var_names if g not in exclude_genes]
        counts = counts[:, keep_genes]
        counts.layers["counts"] = counts.X.copy()

        sc.pp.highly_variable_genes(counts, n_top_genes=4000, flavor='seurat_v3')  # Seurat-style HVG selection
        sc.pp.normalize_total(counts, target_sum=1e4)
        sc.pp.log1p(counts)

        # ================= Scrublet (Doublet detection) =================
        print("Running Scrublet...")

        # Use raw counts for scrublet
        counts_matrix = counts.raw.X

        # If sparse matrix, convert to CSR format (Scrublet works better with CSR)
        if not sp.isspmatrix_csr(counts_matrix):
            counts_matrix = counts_matrix.tocsr()

        # Scrublet initiation and doublet detection
        scrub = scr.Scrublet(counts_matrix)  # Initialize Scrublet
        doublet_scores, predicted_doublets = scrub.scrub_doublets()  # Predict doublets

        # Save doublet scores and predictions in obs
        counts.obs['doublet_score'] = doublet_scores
        counts.obs['predicted_doublet'] = predicted_doublets

        #afterQC_n = np.sum(counts.obs['predicted_doublet'] == False)  # Number of singlets
        afterQC_n = counts.n_obs
        afterQC_g = counts.n_vars  # Number of genes

        # QC Summary
        summary_df = pd.DataFrame({
            'sample': [sample],
            'BeforeControl_cells': [beforeQC_n],
            'BeforeControl_genes': [beforeQC_g],
            'AfterControl_cells': [afterQC_n],
            'AfterControl_genes': [afterQC_g],
            'median_nFeature': [np.median(counts.obs_vector('n_genes'))],
            'median_nCount': [np.median(counts.obs_vector('n_counts'))],
            'median_mt': [np.median(counts.obs['percent_mt'])],
            'median_ribo': [np.median(counts.obs['percent_ribo'])]
        })

        # Save the processed data to .h5ad format (compressed)
        out_path = os.path.join(outdir, f"{sample}.h5ad")
        counts.write(out_path, compression='gzip')

        # save the processed data to 10X matrix
        counts.var_names_make_unique()
        X = counts.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X)

        X = X.T
        sample_outdir = os.path.join(outdir, f"filter{sample}")
        os.makedirs(sample_outdir, exist_ok=True)

        obs_path = os.path.join(sample_outdir, "cell_metadata.xlsx")
        counts.obs.to_excel(obs_path, index=True)

        var_path = os.path.join(sample_outdir, "gene_metadata.xlsx")
        counts.var.to_excel(var_path, index=True)

        mtx_path = os.path.join(sample_outdir, "matrix.mtx")
        mmwrite(mtx_path, X)
        with open(mtx_path, 'rb') as f_in, gzip.open(mtx_path + ".gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(mtx_path)

        barcodes_path = os.path.join(sample_outdir, "barcodes.tsv.gz")
        with gzip.open(barcodes_path, 'wt') as f:
            f.write("\n".join(counts.obs_names) + "\n")

        features_path = os.path.join(sample_outdir, "features.tsv.gz")
        with gzip.open(features_path, 'wt') as f:
            for gene in counts.var_names:
                f.write(f"{gene}\t{gene}\n")

        # save the raw counts of filter data to 10X matrix
        counts.var_names_make_unique()
        #mask = counts.raw.var_names.isin(keep_genes)
        #X = counts.raw.X[:, mask]
        X = counts.layers['counts']
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
            f.write("\n".join(counts.obs_names) + "\n")

        features_path = os.path.join(sample_outdir, "features.tsv.gz")
        with gzip.open(features_path, 'wt') as f:
            for gene in counts.var_names:
                f.write(f"{gene}\t{gene}\n")

        # Add to the results dictionary
        res_dict[sample] = {'summary': summary_df}

    # Save summary as a dataframe and output to Excel
    summary_df_all = pd.concat([res_dict[sample]['summary'] for sample in samples], ignore_index=True)
    summary_df_all = summary_df_all.merge(metainfo[['dataset', 'chemistry']],left_on='sample',right_on='dataset',how='left')
    summary_df_all = summary_df_all.drop(columns='dataset')
    summary_df_all.to_excel(os.path.join(outdir, "QC_summary.xlsx"), index=False)

    # End time
    end_time = time.time()
    print(f"\n===== QC done =====\nTotal runtime: {round((end_time - start_time) / 60, 2)} minutes")
prefilter_raw(samples, fold2, outdir)

# reload all h5ad files to calculate 99th percentile of scroublet score across all datasets
import os
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import time
input = "E:\\AID cohort\\code\\result\\"
samples=["137ACP5","84ADA2","159ADA2","163COPA","88GNAI2","78IL1R1","79IL1R1_2","156ISG15","157LACC1","90OGFRL1","36PLCG1","37PLCG2","108PLD4","109PLD4","110PLD4","155PSMD12","138RIPK1","139RIPK1","140RIPK1","141RIPK1","152RNF31","149SOCS1","162STING","69TBXAS1","70TBXAS1","89TLR7","82TNFAIP3","83TNFAIP3","154TRAF3","151USP18","134ELF4","115KRAS","81NLRC4","91TLR8","148UBA1","1IRAK2","3IRAK2","165TREX1","168TREX1","169TREX1","175IFIH1","176IFIH1","177IFIH1","178ISG15","179STAT4","180SLC7A7","181NOD2","12HC","17HC","22HC","75HC","76HC","77HC","158HC","160HC","161HC","9HC","HC_GSE139324_GSM4138162","HC_GSE139324_GSM4138163","HC_GSE139324_GSM4138164","HC_GSE139324_GSM4138165","HC_GSE139324_GSM4138166","HC_GSE139324_GSM4138167","HC_10XGenomics10K","HC_GSE156989_GSM4749756","HC_GSE156989_GSM4749762","HC_GSE156989_GSM4749768","HC_GSE199445_GSM5973143","HC_GSE199445_GSM5973144","HC_GSE199445_GSM5973145","HC_GSE199445_GSM5973146","HC_GSE168732_GSM5160432","HC_GSE168732_GSM5160434","HC_GSE168732_GSM5160435","HC_CELLxGENE","SLE_GSE156989_GSM4749774","SLE_GSE156989_GSM4749779","SLE_GSE156989_GSM4749784","SLE_GSE142016_GSM4217718","SLE_GSE142016_GSM4217719","SLE_GSE142016_GSM4217720","SLE_GSE263931_GSM8207595","SLE_GSE263931_GSM8207597","SLE_GSE263931_GSM8207599","SLE_GSE263931_GSM8207601","SLE_GSE263931_GSM8207603","SLE_GSE263931_GSM8207605","SLE_GSE263931_GSM8207607","SLE_GSE224198_GSM7017326","SLE_GSE224198_GSM7017329","SLE_GSE224198_GSM7017331","SLE_GSE224198_GSM7017334","JIA_GSE205095_GSM6205132","JIA_GSE205095_GSM6205136","JIA_GSE205095_GSM6205138","JIA_GSE205095_GSM6205140","JIA_GSE205095_GSM6205142","JIA_GSE205095_GSM6205144","KD_GSE168732_GSM5160417","KD_GSE168732_GSM5160420","KD_GSE168732_GSM5160422","KD_GSE168732_GSM5160424","KD_GSE168732_GSM5160427","KD_GSE168732_GSM5160430","STING_GSE226598_GSM7079986","STING_GSE226598_GSM7079989","STING_GSE226598_GSM7079991","PLCG2_Figshare","OTULIN_GSE199445_GSM5973147","OTULIN_GSE199445_GSM5973148","OTULIN_GSE199445_GSM5973149"]

def batch_filter_by_global_doublet(input, samples):
    start_time = time.time()

    adatas_list = []
    all_scores = []

    # --- load all primary h5ad files in one list
    print("--- begining load h5ad files ---")
    for sample in samples:
        file_path = os.path.join(input, f"{sample}.h5ad")
        if os.path.exists(file_path):
            tmp_adata = sc.read_h5ad(file_path)
            tmp_adata.obs['samples'] = sample
            adatas_list.append(tmp_adata)
            all_scores.extend(tmp_adata.obs['doublet_score'].tolist())
        else:
            print(f"warning: there no files exits in {file_path}")

    # --- calculating 95th percentile of doublet scores across all datasets ---
    global_threshold = np.nanquantile(all_scores, 0.95)
    print(f"The threshold of doublet score was {global_threshold:.4f}")

    filtered_adatas = [] # saving h5ad after removeing doublets
    stats_list = []  # statistic data summary

    for adata in adatas_list:
        s_name = adata.obs['samples'].iloc[0]
        n_before = adata.n_obs
        adata_filtered = adata[adata.obs['doublet_score'] <= global_threshold].copy()
        n_after = adata_filtered.n_obs
        n_removed = n_before - n_after
        removed_pct = (n_removed / n_before * 100) if n_before > 0 else 0


        stats_list.append({
            'Sample': s_name,
            'Cells_Before': n_before,
            'Doublets_Removed': n_removed,
            'Cells_After': n_after,
            'Removal_Rate_Percent': round(removed_pct, 2),
            'Global_Threshold_Used': round(global_threshold, 4)
        })

        filtered_adatas.append(adata_filtered)

    stats_df = pd.DataFrame(stats_list)

    stats_excel_path = os.path.join(input, "Doublet_Filtering_Detailed_Stats.xlsx")
    stats_df.to_excel(stats_excel_path, index=False)
    print(f"cells counting done in：{stats_excel_path}")

    # --- combining all h5ad for next integration
    print("--- Integrating staring ---")
    combined_adata = ad.concat(filtered_adatas, index_unique="-", join='outer')


    #sc.pp.highly_variable_genes(
    #    combined_adata,
    #    n_top_genes=4000,
    #    flavor='seurat_v3',
    #    layer="counts"
    #)

    # save the combined adata list
    final_output = os.path.join(input, "combined_all_samples_post_doublet_QC.h5ad")
    combined_adata.write(final_output, compression='gzip')

    print(f"Combining all datasets, and the total cell number was: {combined_adata.n_obs}")

    end_time = time.time()
    print(f"\n===== QC done =====\nTotal runtime: {round((end_time - start_time) / 60, 2)} minutes")
    return combined_adata
adata_final = batch_filter_by_global_doublet(input,samples)
