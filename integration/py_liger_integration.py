# https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html

import argparse

import numpy as np
import scanpy as sc
try:
    import pyliger as pl
except ImportError:
    print("install pyliger first on Linux")



def pyl_integration(h5ad_file):
    adata = sc.read_h5ad(h5ad_file)
    batch_cats = []
    # Pyliger normalizes by library size with a size factor of 1
    # So here we give it the count data
    # bdata.X = bdata.layers["counts"]
    # List of adata per batch
    adata_list = [adata[adata.obs["sample"] == b].copy() for b in batch_cats]
    for i, ad in enumerate(adata_list):
        ad.uns["sample_name"] = batch_cats[i]
        # Hack to make sure each method uses the same genes
        ad.uns["var_gene_idx"] = np.arange(adata.n_vars)

    print(adata_list)

    liger_data = pl.create_liger(adata_list, remove_missing=False, make_sparse=False)
    liger_data.var_genes = adata.var_names
    pl.select_genes(liger_data)
    pl.scale_not_center(liger_data)
    pl.optimize_ALS(liger_data, k=30)
    pl.quantile_norm(liger_data)

    adata.obsm["LIGER"] = np.zeros((adata.shape[0], liger_data.adata_list[0].obsm["H_norm"].shape[1]))
    for i, b in enumerate(batch_cats):
        adata.obsm["LIGER"][adata.obs.batch == b] = liger_data.adata_list[i].obsm["H_norm"]

    pl.leiden_cluster(liger_data, resolution=0.25)
    pl.run_umap(liger_data, distance='cosine', n_neighbors=30, min_dist=0.3)
    all_plots = pl.plot_by_dataset_and_cluster(liger_data, axis_labels=['UMAP 1', 'UMAP 2'], return_plots=True)
    print(all_plots)
    print(f"liger_data: {liger_data}")
    print(f"adata: {adata}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True)
    # parser.add_argument("-output", required=True)
    args = parser.parse_args()
    input_file = args.input
    # out_dir = args.output
    pyl_integration(input_file)

if __name__ == '__main__':
    main()