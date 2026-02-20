from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import scanpy as sc

def create_venn(datasets):
    venn_dict = {}
    for dataset in datasets:
        genes = set()
        for ds in dataset[0]:
            adata = sc.read_h5ad(ds)
            genes.update(adata.var_names)
        venn_dict[dataset[1]] = genes

    v = venn3(tuple(venn_dict.values()), tuple(venn_dict.keys()), set_colors=("#383a6b", "#cb1f73", "#fcc72d"), alpha=0.8)
    plt.title("Venn diagram for genes in datasets after preprocessing")
    plt.savefig("/Users/felixlang/Documents/Uni/Master/master-thesis/tro_org/analysis/analysis_plots/venn_diagram_datasets_after.png")
    plt.show()


def main():
    datasets = [[["/Users/felixlang/Documents/Uni/Master/master-thesis/database/Shibata/Shibata_fixed_normalized.h5ad"], "Shibata"],
               [["/Users/felixlang/Documents/Uni/Master/master-thesis/database/Shannon_McNeil/Seurat/shannon_trophoblast.h5ad"], "Shannon"],
               [["/Users/felixlang/Documents/Uni/Master/master-thesis/database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene_fixed.h5ad",
                 "/Users/felixlang/Documents/Uni/Master/master-thesis/database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene_fixed.h5ad"], "Arutyunyan"]]
    datasets_processed = [
        [["/Users/felixlang/Documents/Uni/Master/master-thesis/database/prepros_test/GSE241052_ari_org_annotated_fixed_normalized_preprocessed.h5ad"],
         "Shibata"],
        [["/Users/felixlang/Documents/Uni/Master/master-thesis/database/prepros_test/shannon_trophoblast_preprocessed.h5ad"],
         "Shannon"],
        [["/Users/felixlang/Documents/Uni/Master/master-thesis/database/prepros_test/Organoid_PTO_cellxgene_fixed_preprocessed.h5ad",
          "/Users/felixlang/Documents/Uni/Master/master-thesis/database/prepros_test/Organoid_TSC_cellxgene_fixed_preprocessed.h5ad"],
         "Arutyunyan"]]
    create_venn(datasets_processed)

if __name__ == '__main__':
    main()