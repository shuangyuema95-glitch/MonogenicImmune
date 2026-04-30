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

def run_subset_benchmark(adata_merged, batch_key, pcaNumber, subsetcell=2000, n_repeats=10, n_jobs=8):
    """
    Run repeated subset benchmark for integration stability.
    Parameters
    adata_merged : anndata after batch correction
    batch_key : str
    pcaNumber : int
    subsetcell : int
    n_repeats : int
    n_jobs : int
    Returns pd.DataFrame
    """
    pipeline_start_time = time.time()
    np.random.seed(42)
    all_cells = adata_merged.obs_names.values
    fixed_cells = [np.random.choice(all_cells, subsetcell, replace=False) for _ in range(n_repeats)]

    def bench(adata, bk, lk, ek):
        dfs = []
        for ids in fixed_cells:
            sub = adata[adata.obs_names.isin(ids)]
            bm = Benchmarker(sub, bk, lk, [ek], n_jobs=n_jobs)
            bm.benchmark()
            dfs.append(bm.get_results().iloc[[0]])
        return pd.concat(dfs).mean(0).to_frame().T

    df_harmony = bench(adata_merged, batch_key, "leiden_harmony", "X_pca_harmony")
    df_bbknn = bench(adata_merged, batch_key, "leiden_bbknn", "X_umap_BBKNN")
    df_scanorama = bench(adata_merged, batch_key, "leiden_scanorama", "X_scanorama")
    df_scVI = bench(adata_merged, batch_key, "leiden_scVI", "X_scVI")
    df_scalex = bench(adata_merged, batch_key, "scalex_leiden", "X_scalex")

    all_clean = format_benchmark_results(df_harmony, df_bbknn, df_scanorama, df_scVI, df_scalex)
    all_clean.to_csv(f"{pcaNumber}_Subset{subsetcell}_Benchmark_Final.csv", index=False)

    pipeline_end_time = time.time()
    total_elapsed_sec = pipeline_end_time - pipeline_start_time
    total_elapsed_min = total_elapsed_sec / 60
    print("\n========================================")
    print(f"Total benchmark subset time: {total_elapsed_sec:.2f} s / {total_elapsed_min:.2f} min")
    print("========================================")

run_subset_benchmark(adata_merged, batch_key="samples", pcaNumber=50)