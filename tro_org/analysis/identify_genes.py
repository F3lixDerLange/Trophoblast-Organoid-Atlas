import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

FS = 20

def identify_batch_specific_genes(adata, batch_key="batch", hi=0.1, lo=0.01):
    batches = list(pd.unique(adata.obs[batch_key]))
    genes = adata.var_names

    # mean expression matrix: rows=batches, cols=genes
    mean_mat = []
    for b in batches:
        Xb = adata[adata.obs[batch_key] == b].X
        mean_mat.append(np.asarray(Xb.mean(axis=0)).ravel())

    mean_df = pd.DataFrame(mean_mat, index=batches, columns=genes)

    unique_batch_genes = {}
    for b in batches:
        mean_b = mean_df.loc[b]
        max_other = mean_df.drop(index=b).max(axis=0)  # <-- key change!
        unique = genes[(mean_b > hi) & (max_other < lo)]
        unique_batch_genes[b] = list(unique)

    print("true")
    # for key, val in unique_batch_genes.items():
        # print(key, val)

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

    print("batch specific genes")
    print(expr_df.index.tolist())

    mat = expr_df.astype(float)
    batch_labels = [str(b).replace("Arutyunyan_", "") for b in mat.index]

    plt.figure(figsize=(12, 6))

    ax = sns.heatmap(
        mat,
        cmap="viridis",
        cbar=True
    )

    ax.set_yticklabels(batch_labels)
    ax.set_xticks(np.arange(mat.shape[1]) + 0.5)
    ax.set_xticklabels(mat.columns, rotation=45, ha="right")

    plt.ylabel("Batch" ,fontsize=FS)
    plt.xlabel("Genes", fontsize=FS)
    plt.title(f"Mean expression of batch-specific genes (n={mat.shape[1]})", fontsize=FS)
    plt.tight_layout()
    plt.savefig("figures/batch_specific_genes.png", dpi=300)
    plt.show()

def identify_shared_genes_all_batches(adata, batch_key="batch", min_mean=0.2, top_n=80):
    batches = list(pd.unique(adata.obs[batch_key]))
    var_names = adata.var_names

    means = []
    for b in batches:
        adata_b = adata[adata.obs[batch_key] == b]
        mean_b = np.asarray(adata_b.X.mean(axis=0)).ravel()
        means.append(mean_b)

    mean_mat = np.vstack(means)
    mean_df = pd.DataFrame(mean_mat, index=batches, columns=var_names)

    shared_mask = (mean_df >= min_mean).all(axis=0)
    shared_genes = mean_df.columns[shared_mask]

    if top_n is not None and len(shared_genes) > top_n:
        overall = mean_df[shared_genes].mean(axis=0).sort_values(ascending=False)
        shared_genes = overall.head(top_n).index

    expr_df = mean_df[shared_genes]

    print(f"Shared genes in all batches: {len(shared_genes)}")
    #print(expr_df.columns.tolist())
    #print(expr_df)

    mat = expr_df.T.astype(float)

    batch_labels = [str(b).replace("Arutyunyan_", "") for b in mat.columns]
    plt.figure(figsize=(8, 14))
    ax = sns.heatmap(
        mat,
        cmap="viridis",
    )
    plt.title("Mean expression of genes expressed in all batches", fontsize=FS, pad=15)
    plt.ylabel("Genes", fontsize=FS)
    plt.xlabel("Batch", fontsize=FS)
    plt.xticks(rotation=45, ha="right")
    ax.set_yticks(np.arange(len(expr_df.columns)) + 0.5)
    ax.set_yticklabels(expr_df.columns)
    ax.set_xticklabels(batch_labels)
    plt.tight_layout()
    plt.savefig("figures/genes_expressed_all_batches.png", dpi=300)
    plt.show()

def main():
    hvg_file = "/Users/felixlang/Documents/Uni/Master/master-thesis/tro_org/benchmark/benchmark_plots/final/merged_integration/merged_integration_V2_integrated.h5ad"
    adata = sc.read_h5ad(hvg_file)

    identify_batch_specific_genes(adata)
    identify_shared_genes_all_batches(adata)

if __name__ == '__main__':
    main()