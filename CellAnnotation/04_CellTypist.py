import os
import time
import celltypist
import scanpy as sc
import pandas as pd

start_time = time.time()

os.chdir("/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50")
adata = sc.read_h5ad("716_pbmc1.h5ad")

model_high = "/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/Immune_All_High.pkl"
model_low = "/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/Immune_All_Low.pkl"


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

results = celltypist.annotate(adata, model=model_low, majority_voting=True)
adata.obs["celltypist_label"] = results.predicted_labels["majority_voting"]

meta_df = adata.obs
meta_df.to_csv("pbmc1_celltypist_metadata.csv")

end_time = time.time()
total_cost = round(end_time - start_time, 2)
print("metadata has been done! Total running time:", total_cost, "seconds")