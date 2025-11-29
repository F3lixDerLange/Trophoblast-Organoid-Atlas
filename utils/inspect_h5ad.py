import sys

import scanpy as sc
import pandas as pd
import os
import analysis.utils as utils

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
    label_col = "cell_annotation" if "cell_annotation" in adata.obs else "celltype"
    print(adata.obs[label_col].value_counts())
    print("X_pca:", adata.obsm.get("X_pca", None).shape if "X_pca" in adata.obsm else None)
    print("PCs:", adata.varm.get("PCs", None).shape if "PCs" in adata.varm else None)

    # utils.normalize_data(adata)

    count_matrix_df = pd.DataFrame(
        adata.X.toarray().T,
        index=adata.var_names,  # Set gene names as the rows names
        columns=adata.obs_names  # Set cell barcodes as the column names
    )

    print("--- DataFrame (Genes in Rows, Cells in Columns) ---")
    print(count_matrix_df.iloc[:40, :5])
    print(f"\nDataFrame-Dimensionen: {count_matrix_df.shape}")
    print("Global min:", count_matrix_df.values.min())
    print("Global max:", count_matrix_df.values.max())

    if "cellxgene" in os.path.basename(h5ad_file):
        for i in range(14):
            cols = [f"gene_ids-{i}", f"feature_types-{i}", f"genome-{i}", f"n_cells-{i}"]
            print(f"\n--- Dataset {i} --- Shape: {adata.var[cols].shape}")
            print(adata.var[cols].head())


def main():
    # h5ad_file = "processed_data/Shibata_Karvas_Shannon_Baltayeva_merged_hvg.h5ad"
    # h5ad_file = "database/Shibata/Shibata_EMO6_hESC/GSM7714458_EMO6_hor_merged.h5ad"
    # h5ad_file = "database/Shibata/GSE241052_ari_org.annotated.h5ad"
    h5ad_file = "database/Shibata/GSE241052_ari_org_annotated_fixed_normalized.h5ad"
    # h5ad_file = "database/Shannon_McNeil/Shannon_McNeil_TBp_EVT_D/GSM6664615_DPT_merged.h5ad"
    # h5ad_file = "database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene_fixed.h5ad"
    # h5ad_file = "/Users/felixlang/Downloads/merged.h5ad"
    print_h5ad_info(h5ad_file)

if __name__ == '__main__':
    main()
