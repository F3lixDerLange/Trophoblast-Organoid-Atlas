import pandas as pd
import scanpy as sc
import tro_org.differential_expresison_analysis.dea_deseq2 as dea_deseq2
import tro_org.differential_expresison_analysis.analysis_dataset_genes as analysis_dataset_genes

CANON = {
    "batch": ["batch", "Batch", "sample"],
    "label": ["label", "celltype", "cell_annotation"],
    "n_counts": ["n_counts", "nCount_RNA", "nUMI"],
    "n_genes": ["n_genes"],
    "n_genes_by_counts": ["n_genes_by_counts", "nFeature_RNA"],
    "pct_counts_mt": ["pct_counts_mt", "percent_mt", "percent.mt", "percent_mito"],
    "dataset": ["dataset"],
    "condition": ["condition"],
}

def map_and_keep(adata, canon_dict):
    old_obs = adata.obs.copy()
    new_obs = pd.DataFrame(index=old_obs.index)

    for canon_name, variants in canon_dict.items():
        existing = [v for v in variants if v in old_obs.columns]
        if not existing:
            continue

        combined = old_obs[existing].bfill(axis=1).iloc[:, 0]
        new_obs[canon_name] = combined

    adata.obs = new_obs

    print(adata.obs.nunique())
    print(adata.obs["batch"].unique().tolist())

    metadata = adata.obs[["condition", "dataset", "label"]].copy()
    for col in metadata.columns:
        metadata[col] = metadata[col].astype("category")

    return adata

def merge_adata(adata_files, filter_in_vivo=False):
    adatas = []

    for path, name in adata_files:
        ad = sc.read_h5ad(path)
        if filter_in_vivo and "Type" in ad.obs.columns:
            ad = ad[ad.obs["Type"] != "in vivo"].copy()
        """    if "Model" in ad.obs.columns:
                print(ad.obs["Model"].unique().tolist())
                print(ad.obs["Batch"].unique().tolist())
        else:
            continue"""

        ad.obs["dataset"] = name
        print(name, "layers:", ad.layers.keys())
        adatas.append(ad)
        print(len(ad.var_names))

    adata_raw = sc.concat(
        adatas,
        join="outer",  # keep all genes
        label="dataset",
        keys=[name for _, name in adata_files]
    )

    print(adata_raw)
    print(adata_raw.layers)

    adata_raw.obs["condition"] = (
        adata_raw.obs["Type"]
        .astype(str)
        .str.lower()
        .apply(lambda x: "in_vivo" if "in vivo" in x else "organoid")
    )

    adata_raw.obs["condition"] = adata_raw.obs["condition"].astype("category")

    return map_and_keep(adata_raw, CANON)

def prepare_dataframes(adata, out, case):
    print(adata.layers)
    bulk_adata = sc.get.aggregate(
        adata,
        by=["label", "dataset", "condition"],
        func="sum",
        layer="raw_counts"
    )

    print(bulk_adata)
    print(bulk_adata.X)
    print(bulk_adata.obs.head(n=100))

    counts_df = pd.DataFrame(
        bulk_adata.layers["sum"],
        index=bulk_adata.obs_names,
        columns=bulk_adata.var_names
    ).astype(int)

    metadata = bulk_adata.obs[["condition", "dataset", "label"]].copy()
    for col in metadata.columns:
        metadata[col] = metadata[col].astype("category")

    print(metadata.index == counts_df.index)
    print(counts_df)
    print(metadata)

    return counts_df, metadata

def handle_deseq_for_datasets(adata, out, plot_dir, go_plot_dir):
    deseq_result = {}

    for dataset in adata.obs["dataset"].unique():
        adata_mask, contrast_cond = mask_condition_for_datasets(adata, dataset)
        count_matrix_comp, metadata_comp = prepare_dataframes(adata_mask, out, "comp")
        result = dea_deseq2.run_dea_deseq2(count_matrix_comp, metadata_comp, plot_dir, "condition", contrast_cond)
        deseq_result[dataset] = result

        del count_matrix_comp, metadata_comp

    for key, value in deseq_result.items():
        print(f"{key}: {value}")

    analysis_dataset_genes.dataset_specific_genes_analysis(deseq_result, plot_dir, go_plot_dir)


def mask_condition_for_datasets(adata_raw, dataset):
    """
    Another idea is a One vs the rest approach:
    1. I divide my dataset into "A" (only dataset A) and "not A" (other datasets: B, C, D) saved in condition "is A"
    2. Perform DeSeq2 with condition "is A"
    3. compute stats: DeseqStats(dds, contrast=["is_A", "A", "notA"])
    4. In the result table I can see with the LFC and the adj Pvalue which genes are expressed more in A and which are more expressed in other
    5. do the same for all datasets
    """

    adata_raw.obs["condition"] = adata_raw.obs["dataset"].apply(lambda x: f"{dataset}" if x == dataset else f"not_{dataset}")

    print("--------- MASK CONDITON ----------")
    print(adata_raw.obs["condition"].value_counts())

    return adata_raw, [f"{dataset}", f"not_{dataset}"]


def differential_expression_analysis(h5ad_files, out, plot_dir, go_plot_dir):
    adata_comp = merge_adata(h5ad_files, filter_in_vivo=True)
    handle_deseq_for_datasets(adata_comp, out, plot_dir, go_plot_dir)


    """adata_cond = merge_adata(h5ad_files)
    count_matrix_cond, metadata_cond = prepare_dataframes(adata_cond, out, "condition")
    _ = dea_deseq2.run_dea_deseq2(count_matrix_cond,
                                  metadata_cond,
                                  plot_dir,
                                  "condition",
                                  ["organoid", "in_vivo"])"""


def main():
    h5ad_files = [
        ["database/Shibata/Shibata_fixed_normalized_raw_filter.h5ad", "Shibata"],
        ["database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene_raw_filter.h5ad", "Arutyunyan_PTO"],
        ["database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene_raw_filter.h5ad", "Arutyunyan_TSC"],
        ["database/Shannon_McNeil/Seurat/shannon_trophoblast_raw_filter.h5ad", "Shannon"]
    ]
    out = "tro_org/differential_expresison_analysis/matrix"
    plot_dir = "tro_org/differential_expresison_analysis/dea_plots"
    go_plot_dir = "tro_org/differential_expresison_analysis/GO"

    differential_expression_analysis(h5ad_files, out, plot_dir, go_plot_dir)


if __name__ == '__main__':
    main()