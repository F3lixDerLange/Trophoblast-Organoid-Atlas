import pandas as pd
import scanpy as sc
import scmappy as smp
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
import plot_analysis

from tro_org.analysis import utils


def scmap_comparison(datasets):
    print("Running scmap comparison")

    for i, ref_ds in enumerate(datasets):
        for j, target_ds in enumerate(datasets):
            if i <= j:
                continue  # skip self mapping

            ref_name = ref_ds["name"]
            target_name = target_ds["name"]
            ref_label_col = ref_ds["label_col"]
            target_label_col = target_ds["label_col"]
            ref_adata = ref_ds["adata"]
            target_adata = target_ds["adata"]
            tar_feature_key = next(k for k in target_adata.var.keys() if "feature" in k)

            print(f"\n--- Map {ref_name} -> {target_name} ---")

            ref_adata_common_genes, target_adata_common_genes = utils.common_genes(ref_ds, target_ds)
            sc.pp.highly_variable_genes(ref_adata_common_genes, flavor="cell_ranger")

            # print(f"ref = {ref_adata_common_genes}")
            # print(f"target = {target_adata_common_genes}")
            # print(target_label_col)



            smp.scmap_annotate(target_adata_common_genes,
                               ref_adata_common_genes,
                               ref_label_col,
                               algorithm_flavor="centroid",
                               gene_selection_flavor="HVGs",
                               #n_genes_selected=4000
                               )

            # print(f"ref = {ref_adata_common_genes}")
            # print(f"target = {target_adata_common_genes}")
            label_col = "cell_annotation" if "cell_annotation" in target_adata_common_genes.obs else "celltype"
            # print(target_adata_common_genes.obs[[f"{label_col}", 'scmap_annotation']])
            # print(target_adata_common_genes.obsm['X_umap'])
            # print("\n")

            map_celltype(target_adata_common_genes, label_col, ref_name, target_name)

            plt.rcParams['figure.figsize'] = (12, 6)
            sc.pl.umap(target_adata_common_genes,
                       color=label_col,
                       title=f"{target_name} Umap with {target_name} annotation",
                       legend_loc="right margin"
                       )
            sc.pl.umap(target_adata_common_genes,
                       color="scmap_annotation",
                       title=f"{target_name} with {ref_name} annotation",
                       legend_loc="right margin"
                       )


            smp.scmap_projection(target_adata_common_genes,
                                 ref_adata_common_genes,
                                 "X_umap",
                                 gene_selection_flavor="HVGs",
                                 inplace=True,
                                 key_added='scmap_annotation')

            print(target_adata_common_genes.obsm['X_umap'])


def map_celltype(target_adata_common_genes, label_col, ref_name, target_name):
    mapping, ct = optimal_celltype_mapping(target_adata_common_genes.obs[f"{label_col}"],
                                           target_adata_common_genes.obs["scmap_annotation"])
    print("cross-tabulation (contingency table) - frequency table")
    print(ct)
    print("\n--- Normalized Table (Row-wise: Percentage of Real Type) ---")
    normalized_table = ct.div(ct.sum(axis=1), axis=0).round(2)
    plot_analysis.similarity_df_heatmap(normalized_table, f"ref:{ref_name}_target:{target_name}")
    print(normalized_table)
    print("hungarian algorithm")
    print(mapping)



def optimal_celltype_mapping(real, transferred): # hungarian algorithm
    ct = celltype_contingency(real, transferred)

    # cost matrix = high cost for low counts
    cost = ct.max().max() - ct.values
    row_ind, col_ind = linear_sum_assignment(cost)

    mapping = pd.Series(
        index=ct.index[row_ind],      # real labels
        data=ct.columns[col_ind],     # transferred labels
        name="transferred_match"
    )
    return mapping, ct


def celltype_contingency(real, transferred):
    df = pd.DataFrame({"real": real, "trans": transferred})
    ct = pd.crosstab(df["real"], df["trans"])
    return ct


def main():
    print("install git+https://github.com/gatocor/scmappy.git@affa0c378ad274e27279f60c48a58644f01a8772")
    configs = utils.load_config("analysis/cosine_comp_config.yaml")
    datasets = utils.load_data(configs)
    scmap_comparison(datasets)

if __name__ == '__main__':
    main()