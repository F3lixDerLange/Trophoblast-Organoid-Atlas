import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp

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

def identify_shared_genes_all_batches(adata, batch_key="batch", min_mean=0.2, min_pct=0.4,top_n=80):
    batches = list(pd.unique(adata.obs[batch_key]))
    var_names = adata.var_names

    means = []
    pcts = []

    for b in batches:
        adata_b = adata[adata.obs[batch_key] == b]
        mean_b = np.asarray(adata_b.X.mean(axis=0)).ravel()
        means.append(mean_b)
        pct_b = np.asarray((adata_b.X > 0).mean(axis=0)).ravel()
        pcts.append(pct_b)


    mean_df = pd.DataFrame(np.vstack(means), index=batches, columns=var_names)
    pct_df = pd.DataFrame(np.vstack(pcts), index=batches, columns=var_names)

    shared_mask = (mean_df >= min_mean).all(axis=0) & (pct_df >= min_pct).all(axis=0)
    shared_genes = mean_df.columns[shared_mask]

    print(len(shared_genes))

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

def batch_genes(adata):
    batch_means = []

    for batch in adata.obs["batch"].unique():
        sub = adata[adata.obs["batch"] == batch]

        if sp.issparse(sub.X):
            mean_expr = np.array(sub.X.mean(axis=0)).ravel()
        else:
            mean_expr = sub.X.mean(axis=0)

        batch_means.append(
            pd.DataFrame({
                "mean_expression": mean_expr,
                "batch": batch
            })
        )

    df = pd.concat(batch_means)

    plt.figure(figsize=(6, 10))
    sns.violinplot(data=df, x="batch", y="mean_expression", cut=0)
    plt.ylim(0, 0.2)
    plt.xticks(rotation=90, ha="right")
    plt.title("Mean gene expression per batch")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))

    ax = sns.histplot(
        data=df,
        x="mean_expression",
        hue="batch",
        bins=1000,
        element="step",  # clean overlay
        #stat="density",  # normalize per batch
        common_norm=False
    )

    plt.xlim(0, 0.2)  # zoom
    plt.xlabel("Mean gene expression")
    plt.ylabel("frequency")
    plt.title("Distribution of mean gene expression per batch (0–0.2)")
    sns.move_legend(ax, "center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()

def main():
    hvg_file = "/Users/felixlang/Downloads/merged_integration_final/merged_integration_final_integrated.h5ad"
    adata = sc.read_h5ad(hvg_file)

    #identify_batch_specific_genes(adata)
    identify_shared_genes_all_batches(adata)
    #batch_genes(adata)
    #test

if __name__ == '__main__':
    main()