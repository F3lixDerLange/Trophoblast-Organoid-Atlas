import os.path

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def plot_data(h5ad_file):
    adata = sc.read_h5ad(h5ad_file)

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    print(adata.obs.columns)

    print("x in obs:", "total_counts" in adata.obs)
    print("y in obs:", "n_genes_by_counts" in adata.obs)
    print("color in obs:", "percent_mito" in adata.obs)
    print("color in var:", "percent_mito" in adata.var)

    sc.pl.scatter(
        adata,
        x='total_counts',
        y="n_genes_by_counts",
        color='pct_counts_in_top_500_genes',
        title="QC: Total counts vs Number of genes"
    )

    sc.pl.violin(
        adata,
        keys=["total_counts"],
        groupby="sample",
        jitter=0.4,
        rotation=45
    )

    dir_name = os.path.dirname(h5ad_file)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    sns.histplot(adata.obs["total_counts"], ax=axs[0])
    sns.histplot(adata.obs["n_genes_by_counts"], ax=axs[1])
    fig.suptitle(dir_name)
    plt.tight_layout()
    plt.show()


def main():
    h5ad_file = "database/Shibata/GSE241052_ari_org_annotated_fixed_normalized.h5ad"
    # h5ad_file = "database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene_fixed.h5ad"
    plot_data(h5ad_file)


if __name__ == '__main__':
    main()