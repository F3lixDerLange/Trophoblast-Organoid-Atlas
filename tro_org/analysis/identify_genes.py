import numpy as np
import pandas as pd
import scanpy as sc
from typing import Literal
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors


import plot_analysis


def generate_hvg_per_batch(hvg_file, sc_flavor: Literal["seurat_v3_paper","seurat_v3", "cell_ranger"] = "cell_ranger"):
    adata = sc.read_h5ad(hvg_file)

    identify_batch_seqcific_genes(adata)

    sc.pp.highly_variable_genes(
        adata,
        batch_key="batch",
        n_top_genes=None,  # keep all genes
        flavor=sc_flavor,
        # min_mean=0.1
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
    np.random.seed(42)
    print("Conserved genes:")
    n_batches = len(set(processed_adata.obs["batch"]))
    if sc_flavor == "cell_ranger":
        conserved_hvgs = processed_adata.var.index[processed_adata.var["highly_variable_intersection"]].tolist()
    else:
        conserved_hvgs = processed_adata.var.index[
            (processed_adata.var["highly_variable"]) &
            (processed_adata.var["highly_variable_nbatches"] == n_batches)
        ].tolist()
    print(conserved_hvgs)
    print(len(conserved_hvgs))

    subsample_common_genes = conserved_hvgs if len(conserved_hvgs) <= 50 else np.random.choice(conserved_hvgs, size=50, replace=True)
    mean_expr_gene_batch_common_genes = processed_adata[: , subsample_common_genes]
    mean_expression_per_gene_per_batch(mean_expr_gene_batch_common_genes)

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


def mean_expression_per_gene_per_batch(processed_adata):
    batches = list(set(processed_adata.obs["batch"]))
    batch_means = {}

    for b in batches:
        subset_batch = processed_adata[processed_adata.obs["batch"] == b, :]
        X = subset_batch.X.A if hasattr(subset_batch.X, "A") else subset_batch.X
        batch_means[b] = np.asarray(X.mean(axis=0)).ravel()

    batch_mean_df = pd.DataFrame(batch_means, index=processed_adata.var_names).T

    plot_analysis.common_gene_heatmap(batch_mean_df, f"Common Genes ({batch_mean_df.shape[1]}) per dataset")

    print(batch_mean_df)

def identify_batch_seqcific_genes(adata):
    batches = list(set(adata.obs["batch"]))
    unique_batch_genes = {}

    for batch in batches:
        adata_b = adata[adata.obs["batch"] == batch]

        mean_expr_b = np.asarray(adata_b.X.mean(axis=0)).ravel()

        adata_other = adata[adata.obs["batch"] != batch]
        mean_expr_other = np.asarray(adata_other.X.mean(axis=0)).ravel()

        unique = adata.var_names[
            (mean_expr_b > 0.1) & (mean_expr_other < 0.01)
            ]

        unique_batch_genes[batch] = list(unique)

    print("true")
    for key, val in unique_batch_genes.items():
        print(key, val)

    all_genes = []
    for genes in unique_batch_genes.values():
        for gene in genes:
            all_genes.append(gene)

    expr_df = pd.DataFrame(index=unique_batch_genes.keys(), columns=all_genes) # df: rows = batches, cols = genes

    for batch in unique_batch_genes:
        adata_b = adata[adata.obs["batch"] == batch]

        X = adata_b[:, all_genes].X

        if not isinstance(X, np.ndarray):
            X = X.toarray()

        expr_df.loc[batch] = X.mean(axis=0)

    # plot heatmap
    plt.figure(figsize=(14, 6))
    sns.heatmap(expr_df.astype(float),
                cmap="viridis",
                #annot=True
                )
    plt.title("Mean expression of batch-specific genes")
    plt.xlabel("Genes")
    plt.ylabel("Batch")
    plt.tight_layout()
    plt.savefig("figures/batch_specific_genes.png")
    plt.show()

def main():
    sc_flavor: Literal["seurat_v3_paper","seurat_v3", "cell_ranger"] = "seurat_v3"
    # hvg_file = "processed_data/Shibata_Karvas_Shannon_Baltayeva_merged_hvg.h5ad"
    hvg_file = "processed_data/Shibata_Arutyunyan_merged_hvg.h5ad"
    generate_hvg_per_batch(hvg_file, sc_flavor)

if __name__ == '__main__':
    main()