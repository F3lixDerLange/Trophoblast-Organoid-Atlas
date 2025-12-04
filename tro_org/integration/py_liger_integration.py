# https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html

import argparse
import numpy as np
import scanpy as sc
import tro_org.utils.plot_utils as pu
try:
    import pyliger as pl
except ImportError:
    print("install pyliger first on Linux")



def pyl_integration(adata, out_dir, batch_key, label_key):

    batch_cats = adata.obs[batch_key].astype(str).unique()
    print("Found batches:", batch_cats)

    pu.plot_umap_before_integration(adata, "pyliger", out_dir, batch_key, label_key)

    # tmp = adata.copy()

    # sc.pp.highly_variable_genes(tmp, batch_key=batch_key)# , n_top_genes=4000, flavor="seurat_v3")

    # Keep only HVGs
    # adata = adata[:, tmp.var["highly_variable"]].copy()

    adata_list = [adata[adata.obs[batch_key] == b].copy() for b in batch_cats]
    for i, ad in enumerate(adata_list):
        ad.uns["sample_name"] = batch_cats[i]
        # Hack to make sure each method uses the same genes
        ad.uns["var_gene_idx"] = np.arange(adata.n_vars)
        ad.obs.index.name = "barcodes"
        ad.var.index.name = "gene_names"

    for i, ad in enumerate(adata_list):
        print("===== AnnData", i, "=====")
        print("shape:", ad.shape)
        print("obs_names empty?", ad.obs_names is None or len(ad.obs_names) == 0)
        print("var_names empty?", ad.var_names is None or len(ad.var_names) == 0)
        print(ad.obs_names.names)
        print(ad.var_names.names)

    liger_data = pl.create_liger(adata_list, remove_missing=False, make_sparse=False)
    liger_data.var_genes = adata.var_names
    pl.normalize(liger_data)
    pl.select_genes(liger_data)
    pl.scale_not_center(liger_data)
    pl.optimize_ALS(liger_data, k=30)
    pl.quantile_norm(liger_data)

    adata.obsm["LIGER"] = np.zeros((adata.shape[0], liger_data.adata_list[0].obsm["H_norm"].shape[1]))
    for i, b in enumerate(batch_cats):
        adata.obsm["LIGER"][adata.obs[batch_key] == b] = liger_data.adata_list[i].obsm["H_norm"]

    pl.leiden_cluster(liger_data, resolution=0.25)
    pl.run_umap(liger_data, distance='cosine', n_neighbors=30, min_dist=0.3)
    all_plots = pl.plot_by_dataset_and_cluster(liger_data, axis_labels=['UMAP 1', 'UMAP 2'], return_plots=True)
    print(all_plots)
    print(f"liger_data: {liger_data}")
    print(f"adata: {adata}")

    pu.plot_umap_after_integration(adata, "LIGER", "liger", out_dir, batch_key, label_key)

    return adata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True, help="h5ad file")
    parser.add_argument("-output", required=True)
    parser.add_argument("-bk", "--batch_key", required=True, help="batch_key")
    parser.add_argument("-lk", "--label_key", required=True, help="label_key")
    args = parser.parse_args()
    input_file = args.input
    out_dir = args.output
    batch_key = args.batch_key
    label_key = args.label_key
    pyl_integration(sc.read_h5ad(input_file), out_dir, batch_key, label_key)

if __name__ == '__main__':
    main()