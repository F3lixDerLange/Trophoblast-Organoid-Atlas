#!/usr/bin/env python3

from pathlib import Path
import argparse
import scanpy as sc


def filter_genes(adata):
    print(f"Adata shape pre filter: {adata.shape}")
    allgenes = set(adata.var_names)

    genes_to_remove = ["5S-rRNA", "snoU13"]

    adata = adata[:, ~adata.var_names.isin(genes_to_remove)].copy()
    adata = adata[:, ~adata.var_names.str.startswith("hsa-mir")].copy()
    adata = adata[:, ~adata.var_names.str.startswith("7SK")].copy()

    print(f"Adata shape post filter: {adata.shape}")

    filtered_genes = allgenes.difference(adata.var_names)
    print(f"Filtered genes: {sorted(filtered_genes)}")

    return adata


def filter_invivo_cells(adata, obs_key="Model"):
    if obs_key not in adata.obs.columns:
        raise KeyError(f"Column '{obs_key}' not found in adata.obs")

    print(f"Adata shape before cell filtering: {adata.shape}")

    model_values = adata.obs[obs_key].astype(str).str.lower()
    keep_mask = model_values != "in vivo"
    adata = adata[keep_mask]

    print(f"Adata shape after cell filtering: {adata.shape}")
    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", help="Path to input .h5ad file")
    parser.add_argument("-output", help="Path to output .h5ad file")
    args = parser.parse_args()
    input_path = args.input
    output_path = args.output


    print(f"Reading: {input_path}")
    adata = sc.read_h5ad(input_path)

    adata = filter_genes(adata)
    adata = filter_invivo_cells(adata, obs_key="Model")

    print(f"Writing filtered file to: {output_path}")
    adata.write_h5ad(output_path)
    print("Done")


if __name__ == "__main__":
    main()