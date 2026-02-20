import scanpy as sc
import cytetype

def cytetype_annotation(adata_file, clusters):
    adata = sc.read_h5ad(adata_file)
    print(adata)
    adata.layers["norm"] = adata.X.copy()

    print("Number of clusters: ", adata.obs[clusters].nunique())
    print(adata.obs[clusters].value_counts())
    adata.var['gene_symbols'] = adata.var_names # shibata "features"

    sc.tl.rank_genes_groups(adata, groupby=clusters, method='t-test', layer="norm", use_raw=False, key_added='rank_genes_'+clusters)
    print(adata)
    print(adata.var['gene_symbols'])

    annotator = cytetype.CyteType(adata, group_key=clusters, rank_key='rank_genes_'+clusters)

    context_shibata = """
    Single-cell transcriptomics data, from trophoblast organoids from human samples.
    The cells were from hormone-responsive endometrial organoids (EMO), termed apical-out (AO)–EMO.
    The patients were either treated or untreated.
    The data was generated using 10X genomics Chromium Controller.
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

    adata = annotator.run(study_context=context_PTO)

    sc.pl.embedding(adata, basis='umap', color=f'cytetype_annotation_{clusters}',)

    cols = ['cell_annotation',
            f'cytetype_annotation_{clusters}',
            f'cytetype_cellOntologyTerm_{clusters}',
            f'cytetype_cellOntologyTermID_{clusters}',
            f'cytetype_cellState_{clusters}'
            ]
    print(adata.obs[cols])



def main():
    input_file_shibata = "database/Shibata/Shibata_fixed_normalized.h5ad"
    input_file_PTO = "database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene.h5ad"
    input_file_TSC = "database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene.h5ad"
    clusters_shi = "celltype"
    clusters_aru = "cell_annotation"
    cytetype_annotation(input_file_PTO, clusters_aru)


if __name__ == '__main__':
    main()