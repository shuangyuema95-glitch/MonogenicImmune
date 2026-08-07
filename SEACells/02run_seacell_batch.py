import scanpy as sc
import os
import pandas as pd
from seacell_core import run_seacell_single


def main():
    """
    Run SEACells metacell computation for each sample separately from integrated h5ad object.

    Parameters
    ----------
    main_adata_path : str
        Path of merged preprocessed h5ad file, obs must contain 'samples' column.
    work_dir : str
        Working directory for SEACells output, script will cd into this path.
    cell_per_seacell : int
        Target number of cells per SEACells metacell.
    """
    main_adata_path = "/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/716_pbmc1.h5ad"
    work_dir = "/public/home/yuxiaomingroup/msy/work/AID/juzhen/clean/result/50/SEACELLS"
    cell_per_seacell = 50
    os.chdir(work_dir)

    adata_main = sc.read_h5ad(main_adata_path)
    sample_id_list = adata_main.obs["samples"].unique().tolist()

    for sid in sample_id_list:
        print(f"Process sample: {sid}")
        ad_sample = adata_main[adata_main.obs["samples"] == sid].copy()
        run_seacell_single(ad_sample, sc_per_cell=cell_per_seacell)
        del ad_sample

    print("All samples finished.")


if __name__ == "__main__":
    main()
