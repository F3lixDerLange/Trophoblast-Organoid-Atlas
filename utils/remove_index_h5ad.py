import argparse
import os.path

import scanpy as sc


def remove_index_h5ad(input_file):
    adata = sc.read_h5ad(input_file)
    print("Has _index in obs?", "_index" in adata.obs.columns)
    print("Has _index in var?", "_index" in adata.var.columns)

    if adata.raw is not None:
        print("Has _index in raw.var?", "_index" in adata.raw.var.columns)
    else:
        print("adata.raw is None")

    if "_index" in adata.raw.var.columns:
        adata.raw.var.rename(columns={"_index": "raw_old_index"}, inplace=True)
    print("raw.var:", adata.raw.var.columns)

    adata.write_h5ad(f"{input_file}_fixed.h5ad")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True)
    args = parser.parse_args()
    input_file = args.input
    remove_index_h5ad(input_file)

if __name__ == '__main__':
    main()