import pandas as pd
import scanpy as sc
from typing import Literal


def generate_hvg_per_batch(hvg_file, sc_flavor: Literal["seurat_v3_paper","seurat_v3", "cell_ranger"] = "cell_ranger"):
    adata = sc.read_h5ad(hvg_file)

    sc.pp.highly_variable_genes(
        adata,
        batch_key="batch",
        n_top_genes=None,  # keep all genes
        flavor=sc_flavor
    )
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    print(adata)
    print(adata.obs.head())
    print(adata.var.head())
    #print(adata.uns.values())

    identify_conserved_genes(adata, sc_flavor=sc_flavor)
    identify_batch_specific_genes(adata, sc_flavor=sc_flavor)

def identify_conserved_genes(processed_adata, sc_flavor):
    print("Conserved genes:")
    n_batches = processed_adata.obs["batch"].nunique()
    if sc_flavor == "cell_ranger":
        conserved_hvgs = processed_adata.var.index[processed_adata.var["highly_variable_intersection"]].tolist()
    else:
        conserved_hvgs = processed_adata.var.index[
            (processed_adata.var["highly_variable"]) &
            (processed_adata.var["highly_variable_nbatches"] == n_batches)
        ].tolist()
    print(conserved_hvgs)
    print(len(conserved_hvgs))

def identify_batch_specific_genes(processed_adata, sc_flavor):
    print("\nBatch Specific Genes (not in all batches)")
    n_batches = processed_adata.obs["batch"].nunique()

    batch_specific = processed_adata.var.index[
        (processed_adata.var["highly_variable"]) &
        (processed_adata.var["highly_variable_nbatches"] < n_batches)
        ].tolist()
    print(batch_specific)
    print(len(batch_specific))

    print("\nBatch Specific Genes (only in one batche)")
    strong_batch_specific = processed_adata.var.index[
        (processed_adata.var["highly_variable"]) & (processed_adata.var["highly_variable_nbatches"] == 1)
        ].tolist()

    print(strong_batch_specific)
    print(len(strong_batch_specific))

    print("\nBatch Specific Genes per Batch")
    batches = set(processed_adata.obs["batch"])
    batch_hvgs = {}
    for batch in batches:
        adata_b = processed_adata[processed_adata.obs["batch"] == batch].copy()
        sc.pp.highly_variable_genes(adata_b, flavor=sc_flavor)
        batch_hvgs[batch] = set(adata_b.var.index[adata_b.var["highly_variable"]])


    batch_to_genes = {
        batch: sorted(list(hvgs.intersection(strong_batch_specific)))
        for batch, hvgs in batch_hvgs.items()
    }

    for key, val in batch_to_genes.items():
        print(key, val)

def main():
    sc_flavor: Literal["seurat_v3_paper","seurat_v3", "cell_ranger"] = "seurat_v3"
    hvg_file = "processed_data/merged_hvg.h5ad"
    generate_hvg_per_batch(hvg_file, sc_flavor)

if __name__ == '__main__':
    main()