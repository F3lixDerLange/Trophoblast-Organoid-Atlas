#Tutorial: https://www.nature.com/articles/s41596-024-00991-3
import argparse

import numpy as np
import scanpy as sc
import scanorama as scn
import tro_org.utils.plot_utils as pu

def scn_integration(adata, output_dir, batch_key, label_key):

    pu.plot_umap_before_integration(adata, "scanorama", output_dir)

    sc.pp.highly_variable_genes(adata, batch_key=batch_key)
    # Select all genes that are variable in at least 2 batches
    var_select = adata.var.highly_variable_nbatches > 2
    var_genes = var_select.index[var_select]
    adatas = [adata[adata.obs[batch_key] == batch_value][:, var_genes].copy() for batch_value in adata.obs[batch_key].unique()]

    #adatas = [adata[adata.obs[batch_key] == batch_value].copy() for batch_value in adata.obs[batch_key].unique()]
    scn.integrate_scanpy(adatas)

    adata_sc = adata.copy()
    scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]
    all_s = np.concatenate(scanorama_int)
    adata_sc.obsm['Scanorama'] = all_s

    pu.plot_umap_after_integration(adata_sc,'Scanorama', "scanorama", output_dir, batch_key, label_key)

    return adata_sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True)
    parser.add_argument("-output", required=True)
    parser.add_argument("-bk", "--batch_key", required=True)
    parser.add_argument("-lk", "--label_key", required=True)
    args = parser.parse_args()
    input_file = args.input
    out_dir = args.output
    batch_key = args.batch_key
    label_key = args.label_key
    scn_integration(sc.read_h5ad(input_file), out_dir, batch_key, label_key)

if __name__ == '__main__':
    main()