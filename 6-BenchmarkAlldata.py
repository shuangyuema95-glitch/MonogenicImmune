import os
import scanpy as sc
import scvi
import scib_metrics
from scib_metrics.benchmark import Benchmarker
import anndata
import pandas as pd
import bbknn
import harmonypy
import scanorama
import torch
from scalex import SCALEX
import numpy as np
import psutil
import time

np.Inf = np.inf
process = psutil.Process()
torch.cuda.empty_cache()
print("The number of GPU: ", torch.cuda.device_count())

input_path = "E:/AID cohort/code/result/"
os.chdir(input_path)

pcaNumber=50
adata_merged = sc.read_h5ad(f"{pcaNumber}_final_merged_object.h5ad")

def format_benchmark_results(df_harmony, df_bbknn, df_scanorama, df_scVI, df_scalex):
    """
    Standardize and clean benchmark results into a unified format.
    Parameters
    df_harmony : pd.DataFrame, Raw harmony benchmark result
    df_bbknn : pd.DataFrame, Raw bbknn benchmark result
    df_scanorama : pd.DataFrame, Raw scanorama benchmark result
    df_scVI : pd.DataFrame, Raw scVI benchmark result
    df_scalex : pd.DataFrame, Raw scalex benchmark result
    Returns, pd.DataFrame, Cleaned, sorted, and formatted benchmark results
    """
    # Clean each result
    df_harmony_clean = df_harmony.iloc[[0]].copy()
    df_bbknn_clean = df_bbknn.iloc[[0]].copy()
    df_scanorama_clean = df_scanorama.iloc[[0]].copy()
    df_scvi_clean = df_scVI.iloc[[0]].copy()
    df_scalex_clean = df_scalex.iloc[[0]].copy()

    # Assign method names
    df_harmony_clean["Embedding"] = "Harmony"
    df_bbknn_clean["Embedding"] = "BBKNN"
    df_scanorama_clean["Embedding"] = "Scanorama"
    df_scvi_clean["Embedding"] = "scVI"
    df_scalex_clean["Embedding"] = "scalex"

    # Concatenate
    all_clean = pd.concat(
        [df_harmony_clean, df_bbknn_clean, df_scanorama_clean, df_scvi_clean, df_scalex_clean],
        axis=0, ignore_index=True)

    # Reorder columns
    cols = ["Embedding"] + [col for col in all_clean.columns if col != "Embedding"]
    all_clean = all_clean[cols]
    # Round numeric columns
    numeric_cols = all_clean.select_dtypes(include=['number']).columns
    all_clean[numeric_cols] = all_clean[numeric_cols].round(9)
    return all_clean

def run_full_benchmark(adata_merged, batch_key, pcaNumber, n_jobs=8):
    """
    Run full-dataset benchmark for all integration methods.
    Parameters
    adata_merged : anndata.AnnData
    batch_key : str
    pcaNumber : int
    n_jobs : int
    Returns pd.DataFrame
    """

    pipeline_start_time = time.time()
    def run_single(adata, bk, lk, ek):
        bm = Benchmarker(adata, batch_key=bk, label_key=lk, embedding_obsm_keys=[ek], n_jobs=n_jobs)
        bm.benchmark()
        return bm.get_results()

    df_harmony = run_single(adata_merged, batch_key, "leiden_harmony", "X_pca_harmony")
    df_bbknn = run_single(adata_merged, batch_key, "leiden_bbknn", "X_umap_BBKNN")
    df_scanorama = run_single(adata_merged, batch_key, "leiden_scanorama", "X_scanorama")
    df_scVI = run_single(adata_merged, batch_key, "leiden_scVI", "X_scVI")
    df_scalex = run_single(adata_merged, batch_key, "scalex_leiden", "X_scalex")

    all_clean = format_benchmark_results(df_harmony, df_bbknn, df_scanorama, df_scVI, df_scalex)
    all_clean.to_csv(f"{pcaNumber}_Full_Benchmark_Final.csv", index=False)

    pipeline_end_time = time.time()
    total_elapsed_sec = pipeline_end_time - pipeline_start_time
    total_elapsed_min = total_elapsed_sec / 60
    print("\n========================================")
    print(f"Total benchmark all data time: {total_elapsed_sec:.2f} s / {total_elapsed_min:.2f} min")
    print("========================================")


run_full_benchmark(adata_merged,batch_key="samples",pcaNumber=50)