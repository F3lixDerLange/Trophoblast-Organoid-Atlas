from pathlib import Path

import scanpy as sc
import pandas as pd
import os
import scib
import numpy as np
from scipy.sparse import issparse

from matplotlib import pyplot as plt
from scib.metrics import isolated_labels, isolated_labels_f1, lisi_graph


def print_h5ad_info(h5ad_file):

    adata = sc.read_h5ad(h5ad_file)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)

    # --- Überprüfung der Hauptstruktur ---
    print("===================================")
    print("AnnData Objekt Struktur:")
    print(adata)

    # --- Überprüfung der Zell-Metadaten (.obs) ---
    print("\nZell-Metadaten (ersten 5 Zeilen):")
    print(adata.obs.head())
    print(adata.obs["pct_counts_mt"].min(), adata.obs["pct_counts_mt"].max(), adata.obs["pct_counts_mt"].mean())




    # --- Überprüfung der Gen-Metadaten (.var) ---
    print("\nGen-Metadaten (ersten 5 Zeilen):")
    print(adata.var.head())

    # print(adata.uns['merged_characteristic_genes'])

    # --- Überprüfung der Zählmatrix-Dimensionen ---
    print(f"\nDimensionen der Matrix: {adata.shape}")
    print(f"Anzahl Zellen (obs): {adata.n_obs}")
    print(f"Anzahl Gene (var): {adata.n_vars}")
    print("===================================")

    print(adata.obs['sample'].value_counts())
    print("Unique cell type annotation:")
    if "cell_annotation" in adata.obs:
        label_col = "cell_annotation"
    elif "gene_annotation" in adata.obs:
        label_col = "gene_annotation"
    elif "celltype" in adata.obs:
        label_col = "celltype"
    elif "label" in adata.obs:
        label_col = "label"

    print(adata.obs[label_col].value_counts())
    print("X_pca:", adata.obsm.get("X_pca", None).shape if "X_pca" in adata.obsm else None)
    print("PCs:", adata.varm.get("PCs", None).shape if "PCs" in adata.varm else None)

    # utils.normalize_data(adata)

    print(diagnose_adata_x(adata))

    print("check for raw counts and normalized counts:")
    raw_count_matrix_df = pd.DataFrame(
        adata.raw.X.toarray().T,
        index=adata.var_names,  # Set gene names as the rows names
        columns=adata.obs_names  # Set cell barcodes as the column names
    )

    print("Raw counts:")
    print(raw_count_matrix_df.iloc[:5, :5])
    print(f"\nDataFrame-Dimensionen: {raw_count_matrix_df.shape}")
    print("Global min:", raw_count_matrix_df.values.min())
    print("Global max:", raw_count_matrix_df.values.max())

    count_matrix_df = pd.DataFrame(
        adata.X.toarray().T,
        index=adata.var_names,  # Set gene names as the rows names
        columns=adata.obs_names  # Set cell barcodes as the column names
    )

    adata.layers["counts"] = adata.raw.X.copy()
    adata.copy().write(f"database/Shibata/GSE241052_ari_org_annotated_fixed_normalized_with_raw.h5ad")
    print("saved")

    print("Normalized counts:")
    print(count_matrix_df.iloc[:5, :5])
    print(f"\nDataFrame-Dimensionen: {count_matrix_df.shape}")
    print("Global min:", count_matrix_df.values.min())
    print("Global max:", count_matrix_df.values.max())


    if "cellxgene" in os.path.basename(h5ad_file):
        for i in range(14):
            cols = [f"gene_ids-{i}", f"feature_types-{i}", f"genome-{i}", f"n_cells-{i}"]
            print(f"\n--- Dataset {i} --- Shape: {adata.var[cols].shape}")
            print(adata.var[cols].head())

    plot(adata, "celltype", "sample", "X_umap")

def plot(adata, label_key, batch_key, osbmkey):
    save_dir = Path(f"ztest_folder")
    sc.settings.figdir = save_dir
    sc.set_figure_params(dpi_save=300, figsize=(10, 10), fontsize=14,  vector_friendly=True)
    sc.pp.neighbors(adata, use_rep=osbmkey)
    sc.pp.pca(adata)
    sc.tl.umap(adata)

    print("----------------------")
    print(adata)
    score = scib.metrics.ari(adata, label_key=label_key, cluster_key="clusters")
    print(score)
    test = scib.metrics.metrics(adata,
                                    adata_int=adata,
                                    batch_key=batch_key,
                                    label_key=label_key,
                                    embed="X_pca",
                                    cluster_nmi=True,
                                    ari_=True,
                                    nmi_=True,
                                    silhouette_=True,
                                    pcr_=True,
                                    hvg_score_=True,
                                    isolated_labels_=True,
                                    isolated_labels_f1_=True,
                                    isolated_labels_asw_=True,
                                    graph_conn_=True,
                                    trajectory_=True,
                                    kBET_=True,
                                    lisi_graph_=True,
                                    ilisi_=True,
                                    clisi_=True
                                    )
    print(test)
    sil = scib.metrics.silhouette_batch(adata, label_key=label_key, batch_key=batch_key, embed="X_pca")
    print(sil)

    for batch in adata.obs[batch_key].unique():
        ax = sc.pl.umap(adata, show=False)
        sc.pl.umap(adata[adata.obs[batch_key] == batch].copy(),
                   color=label_key,
                   ax=ax,
                   title=f"Umap {batch} batch highlighted",
                   save=f"_after_integration_{batch}_highlighted.png")


def diagnose_adata_x(adata):
    X = adata.X

    if issparse(X):
        data = X.data  # only non-zero entries
        has_negative = np.any(data < 0)
        min_val = data.min() if data.size > 0 else 0.0
        max_val = data.max() if data.size > 0 else 0.0
        cell_sums = np.array(X.sum(axis=1)).flatten()
        gene_means = np.array(X.mean(axis=0)).flatten()
    else:
        data = X
        has_negative = np.any(data < 0)
        min_val = data.min()
        max_val = data.max()
        cell_sums = data.sum(axis=1)
        gene_means = data.mean(axis=0)

    return {
        "min": float(min_val),
        "max": float(max_val),
        "has_negative": bool(has_negative),
        "cell_sum_range": (float(cell_sums.min()), float(cell_sums.max())),
        "gene_mean_range": (float(gene_means.min()), float(gene_means.max())),
        "is_sparse": issparse(X),
    }

def main():
    # h5ad_file = "processed_data/Shibata_Karvas_Shannon_Baltayeva_merged_hvg.h5ad"
    # h5ad_file = "database/Shibata/Shibata_EMO6_hESC/GSM7714458_EMO6_hor_merged.h5ad"
    # h5ad_file = "database/Shibata/GSE241052_ari_org.annotated.h5ad"
    h5ad_file = "database/Shibata/GSE241052_ari_org_annotated_fixed_normalized.h5ad"
    # h5ad_file = "/Users/felixlang/Documents/Uni/Master/master-thesis/tro_org/trajectory_analysis/figures/adata/subcluster_Stromal_integrated.h5ad"
    # h5ad_file = "database/Shannon_McNeil/Shannon_McNeil_TBp_EVT_D/GSM6664615_DPT_merged.h5ad"
    h5ad_file = "tro_org/GRN/data/adata/subcluster_Trophoblast_integrated.h5ad"
    # h5ad_file = "tro_org/GRN/data/adata/subcluster_Epithelial_integrated.h5ad"
    # h5ad_file = "tro_org/GRN/data/adata/subcluster_Syncytiotrophoblast_integrated.h5ad"
    #h5ad_file = "tro_org/GRN/data/adata/subcluster_Stromal_integrated.h5ad"
    #h5ad_file = "/Users/felixlang/Downloads/shannon_trophoblast.h5ad"
    #h5ad_file = "database/Shannon_McNeil/divided_files/shannon_in_vivo.h5ad"
    #h5ad_file = "database/Shannon_McNeil/divided_files/shannon_TBpOrg.h5ad"
    #h5ad_file = "database/Shannon_McNeil/divided_files/shannon_TBsOrg.h5ad"
    h5ad_file = "database/integrated_data/TBp_shannon_integrated.h5ad"
    h5ad_file = "processed_data/Shibata_Arutyunyan_merged_hvg.h5ad"
    # h5ad_file = "database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene.h5ad"
    # h5ad_file = "/Users/felixlang/Downloads/merged.h5ad"
    # h5ad_file = "/Users/felixlang/Downloads/pipline_fixed/annotated_samplesheet_scvi/finalized/merged.h5ad"
    # h5ad_file = "/Users/felixlang/Downloads/pipline_single/Shibata_single_samplesheet_scvi/finalized/merged.h5ad"
    print_h5ad_info(h5ad_file)

if __name__ == '__main__':
    main()
