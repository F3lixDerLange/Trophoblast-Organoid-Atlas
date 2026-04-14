import pandas as pd
import tro_org.differential_expresison_analysis.plot_utils as plot_utils
import tro_org.differential_expresison_analysis.GO_analysis as go_analysis

LFC_THR = 2
PADJ_THR = 0.001
BASE_MEAN_THR = 0

def merge_datasets(df_dict, column):
    return pd.concat({k: v[column] for k, v in df_dict.items()},axis=1)

def define_ds_specific_genes(df_dict, lfc_df, padj_df):
    dataset_specific = {}

    for ds in df_dict.keys():
        others = [d for d in df_dict.keys() if d != ds]
        sig_in_X = (padj_df[ds] < PADJ_THR) & (lfc_df[ds].abs() > LFC_THR) & (df_dict[ds]["baseMean"] > BASE_MEAN_THR)
        not_sig_others = (padj_df[others] > PADJ_THR).all(axis=1)

        mask = sig_in_X & not_sig_others

        print(f"-- {ds} mask: {len(mask)}")
        print(mask.sum())

        top_genes = (
            lfc_df.loc[mask, ds]
            .sort_values(ascending=False)
            .head(400)
            .index
        )
        dataset_specific[ds] = top_genes

        print(f"** {ds} mask: {len(top_genes)}")

    print(dataset_specific["Shibata"])
    return dataset_specific


def dataset_specific_genes_analysis(df_dict, plot_dir, go_plot_dir):
    """
    Genes of Dataset are specific if:
    - significant in Dataset X
    - not significant in other datasets
    """

    lfc_df = merge_datasets(df_dict, "log2FoldChange")
    padj_df = merge_datasets(df_dict, "padj")
    ds_specific_genes = define_ds_specific_genes(df_dict, lfc_df, padj_df)
    print(ds_specific_genes)
    for key, val in ds_specific_genes.items():
        print(f"{key}: {len(val)}")

    plot_utils.plot_lfc_heatmap(lfc_df, ds_specific_genes, plot_dir)
    plot_utils.markergenes_venn(ds_specific_genes, plot_dir)
    plot_utils.upsetplot_dataset_specific_genes(ds_specific_genes, plot_dir)
    go_analysis.go_analysis(ds_specific_genes, go_plot_dir)
