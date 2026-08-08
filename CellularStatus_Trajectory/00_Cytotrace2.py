import anndata
import pandas as pd
import os
from cytotrace2_py.cytotrace2_py import cytotrace2

h5ad_path = "716_pbmc1.h5ad"
adata = anndata.read_h5ad(h5ad_path)

raw_adata = adata.raw.to_adata()
raw_adata.var_names = adata.var_names.copy()
del adata

out_root = "cytotrace2"
os.makedirs(out_root, exist_ok=True)
target_celltype = "Monocyte"
cell_mask = raw_adata.obs["Level1"] == target_celltype
target_cells = raw_adata[cell_mask, :].copy()

raw_counts = target_cells.X.toarray()
lib_size = raw_counts.sum(axis=1, keepdims=True)
cpm = (raw_counts / lib_size) * 1e6

expr_df = pd.DataFrame(cpm.T, index=target_cells.var_names, columns=target_cells.obs_names)
annot_df = pd.DataFrame({"phenotype": target_cells.obs["Level1"]}, index=target_cells.obs_names)

expr_path = os.path.join(out_root, "expr_cpm.tsv")
annot_path = os.path.join(out_root, "annot.tsv")
expr_df.to_csv(expr_path, sep="\t", index=True, header=True)
annot_df.to_csv(annot_path, sep="\t", index=True)

res_df = cytotrace2(
    expr_path,
    annotation_path=annot_path,
    species="human",
    output_dir=out_root,
    disable_plotting=True
)

all_ct2 = res_df.copy()
all_ct2.to_csv(os.path.join(out_root, "cytoTRACE2_all_Monocyte.tsv"), sep="\t", index=False)

all_ct2 = all_ct2.set_index("cell_id")
raw_adata.obs.loc[all_ct2.index, "cytotrace2_score"] = all_ct2["cytotrace2_score"]
raw_adata.obs.loc[all_ct2.index, "cytotrace2_potency"] = all_ct2["cytotrace2_potency"]
print("Done")
