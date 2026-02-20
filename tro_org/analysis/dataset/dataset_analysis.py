from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def manage_data(h5ad_file, outdir):
    sc.set_figure_params(dpi_save=300, figsize=(10, 8), fontsize=14, vector_friendly=True)

    dataset_dict = merge_dataset_from_same_study(h5ad_file)
    for dataset_name in dataset_dict.keys():
        if len(dataset_dict[dataset_name]) == 1:
            save_dir = Path(f"{outdir}/{dataset_name}")
            sc.settings.figdir = save_dir
            adata = sc.read_h5ad(dataset_dict[dataset_name][0])
            #plot_data(adata, save_dir, dataset_name)

        elif len(dataset_dict[dataset_name]) >= 2:
            save_dir = Path(f"{outdir}/{dataset_name}")
            sc.settings.figdir = save_dir
            temp_adata = []
            temp_names = []
            for i, ds in enumerate(dataset_dict[dataset_name]):
                print(ds)
                temp_adata.append(sc.read_h5ad(ds))
                temp_names.append(f"{dataset_name}_{i}")
            adata_qc = sc.concat(
                temp_adata,
                axis=0,
                join="outer",  # recommended for QC
                label="dataset",
                keys=temp_names,  # automatically handles any number
                merge="same"
            )
            plot_data(adata_qc, save_dir, dataset_name)




def plot_data(adata, save_dir, dataset_name):

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    if "percent_mito" in adata.obs.columns:
        mt_key = "percent_mito"
    elif "percent_mt" in adata.obs.columns:
        mt_key = "percent_mt"
    else:
        mt_key = "percent.mt"

    mt = adata.obs[mt_key]
    if mt.max() <= 1:
        adata.obs[mt_key] = mt * 100
        print("MT percent fixed")

    fig, axs = plt.subplots(1, 4, figsize=(18, 4))

    # n_genes_by_counts → cell-level metric
    # total_counts → cell-level metric
    # percent_mt → cell-level metric
    # n_cells_by_counts → gene-level metric

    sc.pl.violin(
        adata,
        keys="n_genes_by_counts",
        ax=axs[0],
        #jitter=0.4,
        color="#383a6b",
        show=False
    )
    axs[0].set_title("Detected Genes per Cell")
    axs[0].set_ylabel("Number of genes")

    sc.pl.violin(
        adata,
        keys="total_counts",
        ax=axs[1],
        #jitter=0.4,
        color="#cb1f73",
        show=False
    )
    axs[1].set_title("Total UMI Counts per Cell")
    axs[1].set_ylabel("UMI counts")

    sc.pl.violin(
        adata,
        keys=mt_key,
        ax=axs[2],
        #jitter=0.4,
        color="#e03a3c",
        show=False
    )
    axs[2].set_title("Mitochondrial Content")
    axs[2].set_ylabel("Mitochondrial (%)")

    sns.violinplot(
        y=adata.var["n_cells_by_counts"],
        ax=axs[3],
        cut=0,
        color="#ea6d3d",
    )
    axs[3].set_title("Cells per Gene")
    axs[3].set_ylabel("Number of cells")
    axs[3].text(0.5, -0.03, "n_cells_by_counts", ha='center', va='top', transform=axs[3].transAxes)

    fig.suptitle(f"QC metrics {dataset_name}")
    plt.tight_layout()
    plt.savefig(f"{save_dir}/qcmetrics_{dataset_name}.png")
    plt.show()

    sc.pl.scatter(
        adata,
        x='total_counts',
        y="n_genes_by_counts",
        color='pct_counts_in_top_500_genes',
        title="QC: Total counts vs Number of genes",
        save=f"_total_counts_vs_n_genes_by_counts.png"
    )


def merge_dataset_from_same_study(datasets: list):
    dataset_dict = {}
    for dataset in datasets:
        dataset_name = dataset[1].split("_")[0]
        if dataset_name not in dataset_dict:
            dataset_dict[dataset_name] = [dataset[0]]
        else:
            dataset_dict[dataset_name].append(dataset[0])

    print(dataset_dict)
    return dataset_dict


def main():
    h5ad_files = [["database/Shibata/Shibata_fixed_normalized.h5ad", "Shibata"],
                  ["database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene.h5ad", "Arutyunyan_PTO"],
                  ["database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene.h5ad", "Arutyunyan_TSC"],
                  ["database/Shannon_McNeil/Seurat/shannon_trophoblast.h5ad", "Shannon"]]
    savedir = "figures"

    manage_data(h5ad_files, savedir)


if __name__ == '__main__':
    main()