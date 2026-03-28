import scanpy as sc
import cytetype
from matplotlib.pyplot import savefig


def cytetype_annotation(adata_file, clusters):
    adata = sc.read_h5ad(adata_file)
    adata = adata[adata.obs["Type"] != "in vivo"].copy()
    print(adata)
    adata.layers["norm"] = adata.X.copy()

    print("Number of clusters: ", adata.obs[clusters].nunique())
    print(adata.obs[clusters].value_counts())
    adata.var['gene_symbols'] = adata.var_names # shibata "features"
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    sc.tl.rank_genes_groups(adata, groupby=clusters, method='t-test', layer="norm", use_raw=False, key_added='rank_genes_'+clusters)
    print(adata)
    print(adata.var['gene_symbols'])

    annotator = cytetype.CyteType(adata, group_key=clusters, rank_key='rank_genes_'+clusters)

    context_shibata = """
    Single-cell transcriptomics data, from trophoblast organoids from human samples.
    The cells were from hormone-responsive endometrial organoids (EMO), termed apical-out (AO)–EMO.
    The patients were either treated or untreated.
    """

    context_TSC = """
    Single-cell (scRNA-seq) transcriptomics data from human trophoblast stem cell lines (TSC).
    Cells were maintained in TSC medium (TSCM) or differentiated via a 2D protocol on collagen IV in EVT medium (EVTM) over 8 days.
    The data was generated to analyze VCT-like states and differentiation efficiency compared to in vivo tissues.
    """

    context_PTO = """
    Single-cell (scRNA-seq) transcriptomics data from self-renewing primary trophoblast organoids (PTOs).
    The organoids were cultured in maintenance medium (TOM) or induced to differentiate into extravillous trophoblasts (EVT) using specific EVT medium (EVTM), with NRG1 withdrawn in later phases.
    The data captures diverse trophoblast subtypes, including multinucleated syncytiotrophoblast (SCT) and early invasive EVT stages.
    """

    context_shannon = """
    Single-cell RNA sequencing dataset from healthy human trophoblast cells representing two distinct 3D placental in vitro models: 
    primary trophoblast organoids (TBP-Orgs) and stem cell-derived organoids (TBS-Orgs).
    Both organoid models were subjected to two specific treatment conditions: a regenerative maintenance medium to promote self-renewal and an extravillous trophoblast (EVT) differentiation medium supplemented with NRG1
    """

    adata = annotator.run(study_context=context_PTO)

    sc.pl.embedding(adata, basis='umap', color=f'cytetype_annotation_{clusters}',)
    sc.pl.umap(adata, color=f'cytetype_annotation_{clusters}', save=f"_cytetype")

    cols = ['cell_annotation',
            f'cytetype_annotation_{clusters}',
            f'cytetype_cellOntologyTerm_{clusters}',
            f'cytetype_cellOntologyTermID_{clusters}',
            f'cytetype_cellState_{clusters}'
            ]
    print(adata.obs[cols])



def main():
    input_file_shibata = "database/final_data/Shibata_fixed_raw_filter_normalized.h5ad"
    input_file_PTO = "database/final_data/Organoid_PTO_cellxgene_raw_filter_normalized.h5ad"
    input_file_TSC = "database/final_data/Organoid_TSC_cellxgene_raw_filter_normalized.h5ad"
    input_file_shannon = "database/final_data/shannon_trophoblast_raw_filter_normalized.h5ad"
    clusters_shi = "celltype"
    clusters_aru = "cell_annotation"
    cytetype_annotation(input_file_shannon, clusters_shi)


if __name__ == '__main__':
    main()