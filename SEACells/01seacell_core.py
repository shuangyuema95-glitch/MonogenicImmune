import pandas as pd
import scanpy as sc
import numpy as np
import math
import SEACells


def cumulative_explained_variance(adata, expVarThr=None, n_pc=None):
    """
    Calculate cumulative explained variance from PCA result.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing PCA result in adata.uns['pca'].
    expVarThr : float, optional
        Target cumulative explained variance threshold (0‑1).
    n_pc : int, optional
        Reserved parameter, not used in current implementation.

    Returns
    -------
    dict
        Dictionary with 'n_pcs' (number of PCs reaching threshold) and 'expVar' (actual cumulative variance).
    """
    cs = np.cumsum(adata.uns['pca']['variance'] / np.sum(adata.uns['pca']['variance']))
    res = {'n_pcs': None, 'expVar': None}
    if expVarThr is not None:
        nC = np.array(np.where(cs >= expVarThr))[0][0] + 1
        res['n_pcs'] = nC
        res['expVar'] = cs[nC - 1]
    return res


def get_soft_assignments_modified(model):
    """
    Modified function to extract top‑5 soft assignment labels and weights from SEACells model.

    Parameters
    ----------
    model : SEACells.core.SEACells
        Fitted SEACells model instance.

    Returns
    -------
    soft_labels : pandas.DataFrame
        Top5 SEACell assignment labels for each cell.
    weights : numpy.ndarray
        Corresponding assignment weights for top5 SEACells per cell.
    """
    import copy
    A = copy.deepcopy(model.A_.T)
    labels = []
    weights = []
    for _ in range(5):
        l = A.argmax(1)
        seacell_name = [f"SEACell-{i}" for i in l]
        labels.append(seacell_name)
        weights.append(A[np.arange(A.shape[0]), l])
        A[np.arange(A.shape[0]), l] = -1
    weights = np.vstack(weights).T
    labels = np.vstack(labels).T
    soft_labels = pd.DataFrame(labels)
    soft_labels.index = model.ad.obs_names
    return soft_labels, weights


def run_seacell_single(adata_sub, sc_per_cell=50, seed=42):
    np.random.seed(seed)
    sample_id_list = adata_sub.obs["samples"].unique()
    assert len(sample_id_list) == 1
    sample_id = sample_id_list[0]

    raw_adata = adata_sub.raw.to_adata()
    adata_sub.X = raw_adata.X.copy().astype(float)
    del raw_adata

    sc.pp.normalize_per_cell(adata_sub)
    sc.pp.log1p(adata_sub)

    HVG_top = 3000
    if len(adata_sub.obs["chemistry"].unique()) > 1:
        batch_key = "chemistry"
    else:
        batch_key = None

    sc.pp.highly_variable_genes(
        adata=adata_sub,
        batch_key=batch_key,
        flavor="seurat",
        n_top_genes=HVG_top
    )

    sc.tl.pca(adata_sub, svd_solver="arpack", n_comps=50, use_highly_variable=True)
    pc_result = cumulative_explained_variance(adata_sub, expVarThr=0.9)
    n_pcs_keep = pc_result["n_pcs"]
    adata_sub.obsm["X_pca_npcs"] = adata_sub.obsm["X_pca"][:, :n_pcs_keep]

    n_SEACells = math.ceil(adata_sub.n_obs / sc_per_cell)
    build_kernel_on = "X_pca_npcs"
    n_waypoint_eigs = 7

    print(f"Cell count: {adata_sub.n_obs}, target SEACells: {n_SEACells}")

    model = SEACells.core.SEACells(
        adata_sub,
        n_neighbors=10,
        build_kernel_on=build_kernel_on,
        n_SEACells=n_SEACells,
        n_waypoint_eigs=n_waypoint_eigs,
        convergence_epsilon=1e-5
    )
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=10000)

    soft_labels, weights = get_soft_assignments_modified(model)
    max_weight = np.max(weights, axis=1).tolist()
    non_trivial_count = (model.A_.T > 0.1).sum(axis=1)

    adata_sub.obs["SEACell"] = f"{sample_id}_" + adata_sub.obs["SEACell"]
    adata_sub.obs["SEACell_weight"] = max_weight
    adata_sub.obs["SEACell_non_trivial_assig"] = non_trivial_count

    seacell_df = adata_sub.obs[["SEACell", "SEACell_weight", "SEACell_non_trivial_assig"]].copy()
    seacell_df["cell_barcode"] = adata_sub.obs_names

    compactness_df = SEACells.evaluate.compactness(adata_sub, build_kernel_on)
    separation_df = SEACells.evaluate.separation(adata_sub, build_kernel_on, nth_nbr=1)
    seacell_df = seacell_df.merge(compactness_df, on="SEACell", how="left")
    seacell_df = seacell_df.merge(separation_df, on="SEACell", how="left")

    seacell_df = seacell_df.set_index("cell_barcode")

    out_pkl = f"{sample_id}.seacell.pkl"
    seacell_df.to_pickle(out_pkl)
    print(f"Sample {sample_id} finished, output: {out_pkl}")
    return seacell_df
