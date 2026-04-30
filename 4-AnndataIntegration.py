import os
import scanpy as sc
import scvi
import scib_metrics
from scib_metrics.benchmark import Benchmarker
import time
import anndata
import pandas as pd
import bbknn
import harmonypy
import scanorama
import scanpy.external as sce
import torch
import pickle
from scalex import SCALEX
import numpy as np
import psutil

np.Inf = np.inf
process = psutil.Process()
torch.cuda.empty_cache()
print("The number of GPU: ", torch.cuda.device_count())

input_path = "E:/AID cohort/code/result/"
os.chdir(input_path)

hvgtop3000 = pd.read_csv("top3000rank_interhvgs.csv", header=None)
s_genes = pd.read_csv("s_gene.csv", header=0)
g2m_genes = pd.read_csv("g2m_gene.csv", header=0)
hvgs = hvgtop3000[0].tolist()
s_genes = s_genes['gene'].tolist()
g2m_genes = g2m_genes['gene'].tolist()


def get_memory_gb():
    """Get the current RSS memory usage of the process in gigabytes (GB).
        Returns:float: Memory usage in GB."""
    return process.memory_info().rss / (1024 ** 3)

def get_gpu_memory_gb():
    """Get the peak GPU memory allocated by PyTorch in gigabytes (GB).
       Returns:float: GPU memory usage in GB. Returns 0 if CUDA is not available."""
    if torch.cuda.is_available():
        return torch.cuda.max_memory_allocated() / (1024 ** 3)
    else:
        return 0


def export_all_for_R(adata, filename):
    """Extract all UMAP embeddings and metadata from AnnData and save to CSV for R.
    Parameters:adata (AnnData): Input AnnData object containing obs and obsm.
    filename (str): Output CSV file path to save the combined data."""
    df = adata.obs.copy()
    df['cell_id'] = df.index
    for key in list(adata.obsm.keys()):
        if key.startswith('X_umap'):
            emb = adata.obsm[key][:, :2]
            df[f'{key}_1'] = emb[:, 0]
            df[f'{key}_2'] = emb[:, 1]

    df.to_csv(filename, index=False)
    print(f"All coordinates of UMAP and all leiden clusters has been output in: {filename}")


def run_one_integration_pipeline(raw_adata_path="2.h5ad",pcaNumber=50, batch_key="samples",
    n_top_genes=4000,hvg_strategy="individual",use_hvg=True,pklfile="adata_list_for_scalex.pkl"):
    """
     1. Perform data normalization, scaling, PCA, and integration using multiple methods.
     2. Evaluate batch effects and integration quality using various metrics.
     3. export umap coordinates and metadata to CSV at different stages
     """

    adata_merged = sc.read_h5ad(raw_adata_path)
    resource_logs = []
    pipeline_start_time = time.time()

    #####HVG Selection and Re-normalization
    if use_hvg:
        if hvg_strategy == "merged":
            print(f"[HVG] merged strategy (n_top_genes={n_top_genes})")
            sc.pp.highly_variable_genes(adata_merged, n_top_genes=n_top_genes, flavor="seurat_v3", batch_key=batch_key)
        elif hvg_strategy == "individual":
            print("[HVG] individual strategy")
            adata_merged.var["highly_variable"] = adata_merged.var_names.isin(hvgs)

    print(len(adata_merged.var_names[adata_merged.var["highly_variable"]]))

    #####Re-normalization after merging
    sc.pp.normalize_total(adata_merged, target_sum=1e4, layer='counts')  # normalize using raw counts
    sc.pp.log1p(adata_merged)  # log-transform normalized counts
    adata_merged.layers['log_all'] = adata_merged.X.copy()  # store log-normalized full matrix
    adata_merged.raw = adata_merged  # store primary normalization data to a novel adata attribute
    # print(adata_merged.layers['log_all'].max()) #same dims with primary adata_merge.X
    # print(adata_merged.raw.X.max()) #same dims with primary adata_merge.X

    sc.tl.score_genes_cell_cycle(adata_merged, s_genes=s_genes, g2m_genes=g2m_genes, layer='log_all')
    sc.pp.regress_out(adata_merged, keys=['S_score', 'G2M_score', 'percent_mt', 'percent_ribo', 'n_counts'])
    sc.pp.scale(adata_merged, max_value=10)
    adata_merged.layers['scaled'] = adata_merged.X.copy()  # store scaled HVG matrix
    # print(adata_merged.layers['scaled'].max())
    # print(adata_merged.layers['counts'].max())

    #####PCA
    sc.tl.pca(adata_merged, n_comps=pcaNumber, use_highly_variable=True)
    sc.pp.neighbors(adata_merged, use_rep="X_pca", n_neighbors=25)
    sc.tl.leiden(adata_merged, key_added="leiden_rawPCA")  # save to .obs
    sc.tl.umap(adata_merged, random_state=0)  # min_dist=0.4, spread=1.5,
    adata_merged.obsm['X_umap_rawPCA'] = adata_merged.obsm['X_umap']

    #####recheck data structure
    print(adata_merged.layers['log_all'].max())  # same dims with primary adata_merge.X
    print(adata_merged.layers['scaled'].max())  # same dims with primary adata_merge.X
    print(adata_merged.layers['counts'].max())  # same dims with primary adata_merge.X
    print(adata_merged.X.max())  # same with scaled data, at the last step
    print(adata_merged.raw.X.max())  # same with log data
    print(adata_merged.obsm['X_pca'].shape)  # cell number * n_comps
    print(adata_merged.obsm['X_pca'].max())

    # ---------------- 1. Harmony ----------------
    print("Running Harmony...")

    torch.cuda.reset_peak_memory_stats()
    start_time = time.time()
    cpu_start = time.process_time()
    mem_before = get_memory_gb()

    harmony_out = harmonypy.run_harmony(
        adata_merged.obsm['X_pca'].T,
        adata_merged.obs,
        batch_key,
        verbose=True)
    if harmony_out.Z_corr.shape[0] == pcaNumber:  # PCA=50
        adata_merged.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
    else:
        adata_merged.obsm['X_pca_harmony'] = harmony_out.Z_corr

    runtime_sec = time.time() - start_time
    Runtime_Mins = round(runtime_sec / 60, 3)
    CPU_Time_Seconds = round(time.process_time() - cpu_start, 2)
    Memory_Used_GB = round(get_memory_gb() - mem_before, 3)
    GPU_Memory_Used_GB = round(get_gpu_memory_gb(), 3)

    resource_logs.append({"Step": "Harmony","Runtime_Mins": Runtime_Mins,"CPU_Time_Seconds": CPU_Time_Seconds,
        "Memory_Used_GB": Memory_Used_GB,"GPU_Memory_Used_GB": GPU_Memory_Used_GB})

    sc.pp.neighbors(adata_merged, use_rep="X_pca_harmony", n_neighbors=35)
    sc.tl.leiden(adata_merged, key_added="leiden_harmony")  # save to .obs
    sc.tl.umap(adata_merged, random_state=0)  # min_dist=0.4, spread=1.5,
    adata_merged.obsm['X_umap_harmony'] = adata_merged.obsm['X_umap']

    # ---------------- 2. BBKNN ----------------
    print("Running BBKNN...")
    torch.cuda.reset_peak_memory_stats()
    start_time = time.time()
    cpu_start = time.process_time()
    mem_before = get_memory_gb()

    sce.pp.bbknn(adata_merged,batch_key=batch_key)  # BBKNN reconstructs the neighbor graph, no need to run sc.pp.neighbors again

    runtime_sec = time.time() - start_time
    Runtime_Mins = round(runtime_sec / 60, 3)
    CPU_Time_Seconds = round(time.process_time() - cpu_start, 2)
    Memory_Used_GB = round(get_memory_gb() - mem_before, 3)
    GPU_Memory_Used_GB = round(get_gpu_memory_gb(), 3)

    resource_logs.append({"Step": "BBKNN","Runtime_Mins": Runtime_Mins,"CPU_Time_Seconds": CPU_Time_Seconds,
    "Memory_Used_GB": Memory_Used_GB,"GPU_Memory_Used_GB": GPU_Memory_Used_GB})

    sc.tl.leiden(adata_merged, key_added="leiden_bbknn")
    sc.tl.umap(adata_merged, random_state=0)
    adata_merged.obsm['X_umap_BBKNN'] = adata_merged.obsm['X_umap']

    # ---------------- 3. Scanorama ----------------
    print("Running Scanorama...")
    asca = adata_merged[:, adata_merged.var["highly_variable"]].copy()

    torch.cuda.reset_peak_memory_stats()
    start_time = time.time()
    cpu_start = time.process_time()
    mem_before = get_memory_gb()

    sce.pp.scanorama_integrate(asca, key=batch_key)

    runtime_sec = time.time() - start_time
    Runtime_Mins = round(runtime_sec / 60, 3)
    CPU_Time_Seconds = round(time.process_time() - cpu_start, 2)
    Memory_Used_GB = round(get_memory_gb() - mem_before, 3)
    GPU_Memory_Used_GB = round(get_gpu_memory_gb(), 3)

    resource_logs.append({"Step": "Scanorama","Runtime_Mins": Runtime_Mins,"CPU_Time_Seconds": CPU_Time_Seconds,
        "Memory_Used_GB": Memory_Used_GB,"GPU_Memory_Used_GB": GPU_Memory_Used_GB})

    adata_merged.obsm["X_scanorama"] = asca.obsm["X_scanorama"]
    del (asca)
    sc.pp.neighbors(adata_merged, use_rep="X_scanorama", n_neighbors=35)
    sc.tl.leiden(adata_merged, key_added="leiden_scanorama")
    sc.tl.umap(adata_merged, random_state=0)
    adata_merged.obsm["X_umap_scanorama"] = adata_merged.obsm["X_umap"]

    # ---------------- 4. scVI ----------------
    print("Running scVI...")
    ad_v = adata_merged[:, adata_merged.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(ad_v, layer="counts", batch_key="samples")
    vae = scvi.model.SCVI(ad_v, n_layers=2, n_latent=30, gene_likelihood="nb")

    torch.cuda.reset_peak_memory_stats()
    start_time = time.time()
    cpu_start = time.process_time()
    mem_before = get_memory_gb()

    vae.train(accelerator="gpu", devices=1, max_epochs=200, batch_size=20000,early_stopping=True)  # Using early_stopping parameter to prevent the overfitting of model

    runtime_sec = time.time() - start_time
    Runtime_Mins = round(runtime_sec / 60, 3)
    CPU_Time_Seconds = round(time.process_time() - cpu_start, 2)
    Memory_Used_GB = round(get_memory_gb() - mem_before, 3)
    GPU_Memory_Used_GB = round(get_gpu_memory_gb(), 3)

    resource_logs.append({"Step": "scVI","Runtime_Mins": Runtime_Mins,"CPU_Time_Seconds": CPU_Time_Seconds,
        "Memory_Used_GB": Memory_Used_GB,"GPU_Memory_Used_GB": GPU_Memory_Used_GB})


    adata_merged.obsm["X_scVI"] = vae.get_latent_representation()
    del (ad_v)
    sc.pp.neighbors(adata_merged, use_rep="X_scVI")
    sc.tl.leiden(adata_merged, key_added="leiden_scVI")
    sc.tl.umap(adata_merged, random_state=0)
    adata_merged.obsm["X_umap_scVI"] = adata_merged.obsm["X_umap"]

    # ---------------- 5. SCALEX ----------------
    print("Running SCALEX...")
    with open(pklfile, 'rb') as f:
        adata_list = pickle.load(f)

    torch.cuda.reset_peak_memory_stats()
    start_time = time.time()
    cpu_start = time.process_time()
    mem_before = get_memory_gb()

    adata_scalex = SCALEX(data_list=adata_list, batch_key=batch_key, n_top_features=2000,
        batch_size=1024, max_iteration=200, gpu=0,
        seed=124, verbose=True, ignore_umap=True,
        min_features=0, min_cells=0)

    runtime_sec = time.time() - start_time
    Runtime_Mins = round(runtime_sec / 60, 3)
    CPU_Time_Seconds = round(time.process_time() - cpu_start, 2)
    Memory_Used_GB = round(get_memory_gb() - mem_before, 3)
    GPU_Memory_Used_GB = round(get_gpu_memory_gb(), 3)

    resource_logs.append({"Step": "SCALEX","Runtime_Mins": Runtime_Mins,"CPU_Time_Seconds": CPU_Time_Seconds,
        "Memory_Used_GB": Memory_Used_GB,"GPU_Memory_Used_GB": GPU_Memory_Used_GB})


    adata_scalex_aligned = adata_scalex[adata_merged.obs_names].copy()
    print("cell order consistency:", (adata_scalex_aligned.obs_names == adata_merged.obs_names).all())
    adata_merged.obsm["X_scalex"] = adata_scalex_aligned.obsm["latent"]
    sc.pp.neighbors(adata_merged, use_rep="X_scalex")
    sc.tl.leiden(adata_merged, key_added="scalex_leiden")

    # Summary resource usage
    res_df = pd.DataFrame(resource_logs)
    res_df.to_excel(f"{pcaNumber}_ResourceUsage.xlsx", index=False)

    # Saving adata_merged after batch correction
    adata_merged.write_h5ad(f"{pcaNumber}_final_merged_object.h5ad")
    print(f"\n The final anndata after batch correction has been saved in:{pcaNumber}_final_merged_object.h5ad")
    export_all_for_R(adata_merged, f"{pcaNumber}_AllUMAPLeiden_annMerge.csv")


    pipeline_end_time = time.time()
    total_elapsed_sec = pipeline_end_time - pipeline_start_time
    total_elapsed_min = total_elapsed_sec / 60
    print("\n========================================")
    print(f"Total integration pipeline time: {total_elapsed_sec:.2f} s / {total_elapsed_min:.2f} min")
    print("========================================")


run_one_integration_pipeline(raw_adata_path="2.h5ad",pcaNumber=50,
    batch_key="samples",  n_top_genes=4000,
    hvg_strategy="individual",
    use_hvg=True,
    pklfile="adata_list_for_scalex.pkl")


