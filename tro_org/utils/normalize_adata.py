import argparse
import os

import scanpy as sc

def normalize(in_path, out_path):
    print("--- Normalization ---")
    adata = sc.read_h5ad(in_path)
    adata.X = adata.layers["raw_counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    basename = os.path.basename(in_path)
    basename, ext = os.path.splitext(basename)
    adata.write_h5ad(os.path.join(out_path, f"{basename}_normalized.h5ad"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True)
    parser.add_argument("-output", required=True)
    args = parser.parse_args()
    input_file = args.input
    out_file = args.output
    normalize(input_file, out_file)

if __name__ == '__main__':
    main()