#Tutorial: https://www.nature.com/articles/s41596-024-00991-3
import argparse

import numpy as np
import scanpy as sc
import scanorama as scn

def scn_integration(h5ad_file, output_dir):
    adata = sc.read_h5ad(h5ad_file)

    batch_key = "sample"
    label_key = "celltype"


    sc.pp.highly_variable_genes(adata, batch_key=batch_key)
    # Select all genes that are variable in at least 2 batches
    var_select = adata.var.highly_variable_nbatches > 2
    var_genes = var_select.index[var_select]
    adatas = [adata[adata.obs[batch_key] == batch_value][:, var_genes].copy() for batch_value in adata.obs[batch_key].unique()]
    print("save")

    #adatas = [adata[adata.obs[batch_key] == batch_value].copy() for batch_value in adata.obs[batch_key].unique()]
    scn.integrate_scanpy(adatas)

    adata_sc = adata.copy()
    scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]
    all_s = np.concatenate(scanorama_int)
    adata_sc.obsm['Scanorama'] = all_s

    save_file = f"{output_dir}/scanorama_integration.h5ad"
    adata.write_h5ad(save_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True)
    parser.add_argument("-output", required=True)
    args = parser.parse_args()
    input_file = args.input
    out_dir = args.output
    scn_integration(input_file, out_dir)

if __name__ == '__main__':
    main()