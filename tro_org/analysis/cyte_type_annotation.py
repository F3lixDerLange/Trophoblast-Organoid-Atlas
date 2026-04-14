import scanpy as sc
import cytetype
from dataclasses import dataclass

@dataclass
class Dataset:
    path: str
    annotation_key: str
    name: str

def cytetype_annotation(datasets, contexts):
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
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        sc.tl.rank_genes_groups(adata, groupby=ds.annotation_key, method='t-test', layer="norm", use_raw=False, key_added='rank_genes_'+ds.annotation_key)
        print(adata)
        print(adata.var['gene_symbols'])

        annotator = cytetype.CyteType(adata, group_key=ds.annotation_key, rank_key='rank_genes_'+ds.annotation_key)

        adata = annotator.run(study_context=contexts[ds.name])

        sc.pl.embedding(adata, basis='umap', color=f'cytetype_annotation_{ds.annotation_key}',)
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



def main():
    datasets = [
        # Dataset("database/final_data/Shibata_fixed_raw_filter_normalized.h5ad", "celltype", "shibata"),
        Dataset("database/final_data/Organoid_PTO_cellxgene_raw_filter_normalized.h5ad", "cell_annotation",
                "arutyunyan_PTO"),
        Dataset("database/final_data/Organoid_TSC_cellxgene_raw_filter_normalized.h5ad", "cell_annotation",
                "arutyunyan_TSC"),
        # Dataset("database/final_data/shannon_trophoblast_raw_filter_normalized.h5ad", "celltype", "shannon"),
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

    cytetype_annotation(datasets, contexts)


if __name__ == '__main__':
    main()