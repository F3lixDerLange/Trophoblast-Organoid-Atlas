from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import scanpy as sc

def create_venn(datasets):
    venn_dict = {}
    after_pre_pros = False
    for dataset in datasets.keys():
        genes = set()
        for ds in datasets[dataset]:
            if "preprocessed" in ds:
                after_pre_pros = True
            adata = sc.read_h5ad(ds)
            genes.update(adata.var_names)
        venn_dict[dataset] = genes

    state = "after" if after_pre_pros == True else "before"

    v = venn3(tuple(venn_dict.values()), tuple(venn_dict.keys()), set_colors=("#383a6b", "#cb1f73", "#fcc72d"), alpha=0.8)
    plt.title(f"Venn diagram for genes in datasets {state} preprocessing")
    plt.savefig("/Users/felixlang/Documents/Uni/Master/master-thesis/tro_org/analysis/analysis_plots/venn_diagram_datasets_after.png")
    plt.show()


def main():
    datasets = {'Shibata': ['database/Shibata/Shibata_fixed_normalized.h5ad'],
                'Arutyunyan': ['database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene.h5ad', 'database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene.h5ad'],
                'Shannon': ['database/Shannon_McNeil/Seurat/shannon_trophoblast.h5ad']}


    dataset_prepros = {'Shibata': ['database/prepros_test/Shibata_fixed_normalized_preprocessed.h5ad'],
     'Arutyunyan': ['database/prepros_test/Organoid_PTO_cellxgene_preprocessed.h5ad', 'database/prepros_test/Organoid_TSC_cellxgene_preprocessed.h5ad'],
     'Shannon': ['database/prepros_test/shannon_trophoblast_preprocessed.h5ad']}

    create_venn(datasets)

if __name__ == '__main__':
    main()