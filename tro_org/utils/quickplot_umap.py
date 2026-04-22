import math
import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt
import  tro_org.utils.utils as tutils


def quickplot_umap(adata, output):
    sc.set_figure_params(dpi_save=300, fontsize=16, vector_friendly=True)

    print(adata)
    adata = tutils.filter_invivo_cells(adata, "batch")
    embeddings = [ "Unintegrated", "BBKNN", "Combat", "Harmony", "Scanorama", "scVI", "LIGER", "scGlue"]
    # embeddings = ["scVI"]

    """
    n_cells = 2000
    idx = np.random.choice(adata.n_obs, n_cells, replace=False)

    adata = adata[idx].copy()
    """

    for emb in embeddings:
        sc.pp.neighbors(adata, use_rep=emb, key_added=f"neighbors_{emb}")
        sc.tl.umap(adata, neighbors_key=f"neighbors_{emb}", min_dist=0.3, random_state=0)
        adata.obsm[f"X_umap_{emb}"] = adata.obsm["X_umap"].copy()

    umap_by_dataset(adata, embeddings, output)
    umap_by_cell_type(adata, embeddings, output)


def umap_by_dataset(adata, embeddings, output):
    fig, axes = plt.subplots(4, 2, figsize=(11, 14))
    axes = axes.flatten()

    for i, emb in enumerate(embeddings):
        ax = axes[i]
        sc.pl.embedding(
            adata,
            basis=f"X_umap_{emb}",
            color="sample",
            ax=ax,
            show=False,
            title=emb.replace("X_", ""),
        )

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        if i == 6:
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            ax.annotate(
                '', xy=(0.3, -0.05),
                xytext=(-0.06, -0.05),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2)
            )

            ax.annotate(
                '', xy=(-0.06, 0.3),
                xytext=(-0.06, -0.05),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2),
                annotation_clip=False
            )

    plt.tight_layout()
    plt.savefig(f"{output}/umap_all_tools_by_dataset.png")
    plt.show()


def umap_by_cell_type(adata, embeddings, output):
    n = len(embeddings)
    ncols = 3
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 14))
    axes = axes.flatten()

    for i, emb in enumerate(embeddings):
        ax = axes[i]
        legend_loc = "right margin" if i == 7 else None

        sc.pl.embedding(
            adata,
            basis=f"X_umap_{emb}",
            color="label",
            ax=ax,
            show=False,
            title=emb.replace("X_", ""),
            legend_loc=legend_loc,
            legend_fontsize = 11
        )

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        if i == 6:
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            ax.annotate(
                '', xy=(0.3, -0.05),
                xytext=(-0.06, -0.05),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2)
            )

            ax.annotate(
                '', xy=(-0.06, 0.3),
                xytext=(-0.06, -0.05),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2),
                annotation_clip = False
            )

    for j in range(n, len(axes)):
        fig.delaxes(axes[j])

    plt.savefig(f"{output}/umap_all_tools_by_cell_type.png", dpi=300, bbox_inches="tight")
    plt.show()

def quick_umap():
    adata = sc.read_h5ad("/Users/felixlang/Documents/Uni/Master/master-thesis/figures/cytetype/cytetype_annotated_adata.h5ad")

    fig, ax = plt.subplots(figsize=(10, 6))
    sc.pl.embedding(adata,
                    basis='umap',
                    color="label",
                    title=f'Cell Type annotation Final Atlas',
                    legend_loc='on data',
                    legend_fontsize=12,
                    show=False,
                    ax=ax,
                    size=3)

    legend = ax.get_legend()
    if legend is not None:
        handles = legend.legend_handles
        labels = [t.get_text() for t in legend.get_texts()]
        ax.legend(
            handles,
            labels,
            loc='center left',
            bbox_to_anchor=(1.0, 0.5),
            ncol=2,
            fontsize=12,
            frameon=False,
            markerscale=1.2
        )

    plt.tight_layout()
    plt.savefig(
        f'/Users/felixlang/Documents/Uni/Master/master-thesis/figures/cytetype/umap_celltype_final_atals.png',
        dpi=300, bbox_inches='tight')
    plt.show()

def main():
    output = "/Users/felixlang/Documents/Uni/Master/master-thesis/figures/integration_umaps"
    adata_file = "/Users/felixlang/Downloads/merged_integration_final/merged_integration_final_integrated.h5ad"
    adata = sc.read_h5ad(adata_file)
    quickplot_umap(adata, output)
    quick_umap()

if __name__ == '__main__':
    main()