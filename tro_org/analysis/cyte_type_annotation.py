import scanpy as sc
import cytetype
from dataclasses import dataclass
from tro_org.utils.utils import filter_invivo_cells
from matplotlib import pyplot as plt


@dataclass
class Dataset:
    path: str
    annotation_key: str
    name: str

def cytetype_annotation(datasets, contexts, final_atlas):
    sc.set_figure_params(dpi_save=300, vector_friendly=True)
    for ds in datasets:
        adata = sc.read_h5ad(ds.path)
        if "Type" in adata.obs.columns:
            adata = adata[adata.obs["Type"] != "in vivo"].copy()
        print(adata)
        adata.layers["norm"] = adata.X.copy()

        print("Number of clusters: ", adata.obs[ds.annotation_key].nunique())
        print(adata.obs[ds.annotation_key].value_counts())
        adata.var['gene_symbols'] = adata.var_names # shibata "features"

        if final_atlas:
            cytetype_annotation_final_atlas(adata, ds, contexts)
        else:
            cytetype_annotation_dataset(adata, ds, contexts)


def cytetype_annotation_dataset(datasets, contexts):
    sc.set_figure_params(dpi_save=300, vector_friendly=True)
    for ds in datasets:
        adata = sc.read_h5ad(ds.path)
        if "Type" in adata.obs.columns:
            adata = adata[adata.obs["Type"] != "in vivo"].copy()
        print(adata)
        adata.layers["norm"] = adata.X.copy()

        print("Number of clusters: ", adata.obs[ds.annotation_key].nunique())
        print(adata.obs[ds.annotation_key].value_counts())
        adata.var['gene_symbols'] = adata.var_names  # shibata "features"

        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        sc.tl.rank_genes_groups(adata, groupby=ds.annotation_key, method='t-test', layer="norm", use_raw=False,
                                key_added='rank_genes_' + ds.annotation_key)
        print(adata)
        print(adata.var['gene_symbols'])

        annotator = cytetype.CyteType(adata, group_key=ds.annotation_key, rank_key='rank_genes_' + ds.annotation_key)

        adata = annotator.run(study_context=contexts[ds.name])

        sc.pl.embedding(adata, basis='umap', color=f'cytetype_annotation_{ds.annotation_key}', )
        sc.pl.umap(adata,
                   color=ds.annotation_key,
                   title=f'Cell Type Annotation - {ds.name}',
                   legend_fontsize=12,
                   save=f"_{ds.name}_celltype_annotation.png")
        sc.pl.umap(adata,
                   color=f'cytetype_annotation_{ds.annotation_key}',
                   title=f'CyteType annotation {ds.name}',
                   legend_fontsize=12,
                   save=f"_{ds.name}_cytetype.png")

        cols = [f'{ds.annotation_key}',
                f'cytetype_annotation_{ds.annotation_key}',
                f'cytetype_cellOntologyTerm_{ds.annotation_key}',
                f'cytetype_cellOntologyTermID_{ds.annotation_key}',
                f'cytetype_cellState_{ds.annotation_key}'
                ]
        print(adata.obs[cols])



def cytetype_annotation_final_atlas(ds, context):
    sc.set_figure_params(dpi_save=300, vector_friendly=True)

    adata = sc.read_h5ad(ds.path)
    adata = filter_invivo_cells(adata, "batch")
    print(adata)
    adata.layers["norm"] = adata.X.copy()

    print("Number of clusters: ", adata.obs[ds.annotation_key].nunique())
    print(adata.obs[ds.annotation_key].value_counts())
    adata.var['gene_symbols'] = adata.var_names  # shibata "features"

    sc.pp.neighbors(adata, use_rep="X_scvi")
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scvi")
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden_scvi",
        method="t-test",
        layer="norm",
        use_raw=False,
        key_added="rank_genes_leiden_scvi"
    )
    print(adata)

    annotator = cytetype.CyteType(adata, group_key="leiden_scvi", rank_key='rank_genes_leiden_scvi')

    adata = annotator.run(study_context=context)

    sc.pl.embedding(adata, basis='umap', color=f'cytetype_annotation_leiden_scvi', )
    plot_final_atlas(adata,ds, "Cell Type", ds.annotation_key)
    plot_final_atlas(adata, ds, "CyteType", 'cytetype_annotation_leiden_scvi')

def plot_final_atlas(adata,ds, title, color):
    fig, ax = plt.subplots(figsize=(10, 6))
    sc.pl.embedding(adata,
                    basis='umap',
                    color=color,
                    title=f'{title} annotation {ds.name}',
                    legend_loc='right margin',
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
            ncol=1,
            fontsize=12,
            frameon=False,
            markerscale=1.2
        )

    plt.tight_layout()
    plt.savefig(
        f'/Users/felixlang/Documents/Uni/Master/master-thesis/figures/integration_umaps/umap_{ds.name}_{title.replace(" ", "")}.png',
        dpi=300, bbox_inches='tight')
    plt.show()

def main():
    datasets = [
        Dataset("database/final_data/Shibata_fixed_raw_filter_normalized.h5ad", "celltype", "shibata"),
        Dataset("database/final_data/Organoid_PTO_cellxgene_raw_filter_normalized.h5ad", "cell_annotation",
                "arutyunyan_PTO"),
        Dataset("database/final_data/Organoid_TSC_cellxgene_raw_filter_normalized.h5ad", "cell_annotation",
                "arutyunyan_TSC"),
        Dataset("database/final_data/shannon_trophoblast_raw_filter_normalized.h5ad", "celltype", "shannon"),
    ]

    contexts = {
        "shibata": """
                   Single-cell transcriptomics data, from trophoblast organoids from human samples.
                   The cells were from hormone-responsive endometrial organoids (EMO), termed apical-out (AO)–EMO.
                   The patients were either treated or untreated.
                   """,
        "arutyunyan_TSC": """
                   Single-cell (scRNA-seq) transcriptomics data from human trophoblast stem cell lines (TSC).
                   Cells were maintained in TSC medium (TSCM) or differentiated via a 2D protocol on collagen IV in EVT medium (EVTM) over 8 days.
                   The data was generated to analyze VCT-like states and differentiation efficiency compared to in vivo tissues.
                   """,
        "arutyunyan_PTO": """
                   Single-cell (scRNA-seq) transcriptomics data from self-renewing primary trophoblast organoids (PTOs).
                   The organoids were cultured in maintenance medium (TOM) or induced to differentiate into extravillous trophoblasts (EVT) using specific EVT medium (EVTM), with NRG1 withdrawn in later phases.
                   The data captures diverse trophoblast subtypes, including multinucleated syncytiotrophoblast (SCT) and early invasive EVT stages.
                   """,
        "shannon": """
                   Single-cell RNA sequencing dataset from healthy human trophoblast cells representing two distinct 3D placental in vitro models: 
                   primary trophoblast organoids (TBP-Orgs) and stem cell-derived organoids (TBS-Orgs).
                   Both organoid models were subjected to two specific treatment conditions: a regenerative maintenance medium to promote self-renewal and an extravillous trophoblast (EVT) differentiation medium supplemented with NRG1
                   """
    }

    final_atlas_context = """
                    Single-cell RNA sequencing data from human trophoblast organoids, integrated from four independent datasets.
                    The combined atlas comprises primary trophoblast organoids (PTOs), trophoblast stem cell-derived organoids (TSC organoids),
                    hormone-responsive endometrial organoids (EMOs), and both primary and stem cell-derived 3D placental organoid models.
                    Cells represent the major trophoblast lineages of the human placenta, including villous cytotrophoblasts (VCT),
                    extravillous trophoblasts (EVT) and their invasive subtypes, syncytiotrophoblast progenitors (SCTp) and
                    multinucleated syncytiotrophoblasts (SCT), cycling and column cytotrophoblasts (cCTB, CTB),
                    and trophoblast stem cells (TSC). Some datasets additionally contain non-trophoblast cell types from
                    endometrial organoids, including ciliated, stromal, glandular, luminal, and endothelial populations.
                    Organoids were cultured under maintenance conditions promoting self-renewal or under differentiation conditions
                    using EVT medium supplemented with NRG1 to induce extravillous trophoblast differentiation.
                    """


    final_atlas = Dataset("/Users/felixlang/Downloads/merged_integration_final/merged_integration_final_integrated.h5ad",
                          "label", "Final Atlas")

    cytetype_annotation_final_atlas(final_atlas, final_atlas_context)
    #cytetype_annotation_dataset(datasets, contexts)


if __name__ == '__main__':
    main()