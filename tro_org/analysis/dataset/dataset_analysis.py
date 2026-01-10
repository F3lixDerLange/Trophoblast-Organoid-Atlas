import os.path
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def plot_data(h5ad_file, outdir):
    sc.set_figure_params(dpi_save=300, figsize=(10, 8), fontsize=14, vector_friendly=True)
    save_dir = Path(f"{outdir}/{h5ad_file[1]}")
    sc.settings.figdir = save_dir

    adata = sc.read_h5ad(h5ad_file[0])

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
        title="QC: Total counts vs Number of genes",
        save=f"_total_counts_vs_n_genes_by_counts.png"
    )

    sc.pl.violin(
        adata,
        keys=["total_counts"],
        groupby="sample",
        jitter=0.4,
        rotation=45,
        save=f"_total_counts_per_batch.png"
    )

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    sns.histplot(adata.obs["total_counts"], ax=axs[0])
    sns.histplot(adata.obs["n_genes_by_counts"], ax=axs[1])
    fig.suptitle(h5ad_file[1])
    plt.tight_layout()
    plt.savefig(f"{outdir}/{h5ad_file[1]}/total_counts_per_batch.png")
    plt.show()


def main():
    h5ad_files = [["database/Shibata/GSE241052_ari_org_annotated_fixed_normalized.h5ad", "Shibata"],
                  ["database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene.h5ad", "Arutyunyan_PTO"],
                  ["database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene.h5ad", "Arutyunyan_TSC"]]
    # h5ad_file = "database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene_fixed.h5ad"
    savedir = "figures"

    for dataset in h5ad_files:
        plot_data(dataset, savedir)


if __name__ == '__main__':
    main()