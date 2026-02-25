import argparse

import scanpy as sc
import anndata as ad
import numpy as np

def fix(input_file, output_file):
    adata = sc.read_h5ad(input_file)

    """print("Before cleanup:")
    print("  obsm keys:", list(adata.obsm.keys()))
    print("  varm keys:", list(adata.varm.keys()))

    #Drop all PCA / UMAP / TSNE etc. from obsm
    for key in list(adata.obsm.keys()):
        if key.startswith("X_pca") or key.startswith("X_umap") or key.startswith("X_tsne"):
            print("Dropping obsm[%r]" % key)
            del adata.obsm[key]

    #Drop PCA-related loadings from varm
    for key in list(adata.varm.keys()):
        if "PC" in key or "pca" in key.lower():
            print("Dropping varm[%r]" % key)
            del adata.varm[key]

    print("\nAfter cleanup:")
    print("  obsm keys:", list(adata.obsm.keys()))
    print("  varm keys:", list(adata.varm.keys()))"""

    if "log1p" in adata.uns:
        # simplest: drop it entirely (safe for scdownstream QC)
        del adata.uns["log1p"]

    adata.write_h5ad(output_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True)
    parser.add_argument("-output", required=True)
    args = parser.parse_args()
    input_file = args.input
    out_file = args.output
    fix(input_file, out_file)

if __name__ == '__main__':
    main()