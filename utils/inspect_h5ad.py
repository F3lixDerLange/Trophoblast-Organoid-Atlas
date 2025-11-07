import scanpy as sc
import pandas as pd

def print_h5ad_info(h5ad_file):

    # Laden der AnnData-Datei
    adata = sc.read_h5ad(h5ad_file)

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

    # --- Überprüfung der Zählmatrix-Dimensionen ---
    print(f"\nDimensionen der Matrix: {adata.shape}")
    print(f"Anzahl Zellen (obs): {adata.n_obs}")
    print(f"Anzahl Gene (var): {adata.n_vars}")
    print("===================================")

    count_matrix_df = pd.DataFrame(
        adata.X.toarray().T,
        index=adata.var_names,  # Set gene names as the rows names
        columns=adata.obs_names  # Set cell barcodes as the column names
    )

    print("--- DataFrame (Genes in Rows, Cells in Columns) ---")
    print(count_matrix_df.iloc[:40, :5])
    print(f"\nDataFrame-Dimensionen: {count_matrix_df.shape}")


def main():
    h5ad_file = "database/GSE272513_RAW/GSM840378/GSM8403781_merged.h5ad"
    print_h5ad_info(h5ad_file)

if __name__ == '__main__':
    main()
