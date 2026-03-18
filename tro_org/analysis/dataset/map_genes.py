import argparse
import os

import mygene
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import anndata as ad
from matplotlib import pyplot as plt



def handle_duplicate_gene_names(adata):
    if not adata.var_names.has_duplicates:
        print("No duplicate gene names")
        return adata

    func = "mean"
    adata.var["gene"] = adata.var_names

    agg_X = sc.get.aggregate(adata, by="gene", axis="var", func=func)
    X_new = agg_X.layers[func]
    agg_raw = sc.get.aggregate(adata, by="gene", axis="var", func=func, layer="raw_counts")
    raw_new = agg_raw.layers[func]

    adata_dedup = ad.AnnData(
        X=X_new,
        obs=adata.obs.copy(),
        var=agg_X.var.copy(),
        layers={"raw_counts": raw_new},
    )

    print("duplicates left:", adata_dedup.var_names.duplicated().sum())

    print(adata_dedup)
    if sp.issparse(adata_dedup.X):
        adata_dedup.X = adata_dedup.X.tocsr().astype("float32")
    else:
        adata_dedup.X = sp.csr_matrix(adata_dedup.X.astype("float32"))

    print("-------------")
    print(adata)
    print(adata_dedup.var_names)
    print(adata_dedup.var["gene"])
    print("--------------")
    for var_name, gene_col in list(zip(adata_dedup.var["gene"], adata_dedup.var_names))[:50]:
        print(var_name, " | ", gene_col)


    """gene = "RXFP4"

    idx = np.where(adata.var_names == gene)[0]

    print("Column indices:", idx)

    # Look at annotation differences
    print(adata.var.iloc[idx])

    X = adata.X[:, idx]

    if sp.issparse(X):
        X = X.toarray()

    print("First 10 cells:")
    print(X[:50, :])

    print("Are duplicates identical?",
          np.allclose(X[:, 0], X[:, 1]))"""

    return adata_dedup


def filter_genes(adata):
    print(f"Adata shape pre filter: {adata.shape}")
    allgenes = set(adata.var_names)
    genes_to_remove = ['5S-rRNA', "snoU13"] #, '126231', '55872'
    adata = adata[:, ~adata.var_names.isin(genes_to_remove)].copy()
    adata = adata[:, ~adata.var_names.str.startswith("hsa-mir")].copy()
    adata = adata[:, ~adata.var_names.str.startswith("7SK")].copy()
    print(f"Adata shape post filter: {adata.shape}")

    filtered_genes = allgenes.difference(adata.var_names)
    print(f"Filtered genes: {filtered_genes}")

    return adata


def raw_counts_2_layer(adata):
    adata.layers["raw_counts"] = adata.raw.X.copy()
    print("add raw_counts layer to adata")
    return adata


def generate_out(adata_path):
    adata_dir = os.path.dirname(adata_path)
    adata_base = os.path.splitext(os.path.basename(adata_path))[0]
    return os.path.join(adata_dir, f"{adata_base}_preprocessed.h5ad")

def check_expresison(adata, name):
    gene_sums = np.array(adata.X.sum(axis=0)).ravel()
    zero_genes = adata.var_names[gene_sums == 0]
    print("Number of genes with no expression:", len(zero_genes))
    print(zero_genes[:20])

    gene_cell_counts = np.array((adata.X > 0).sum(axis=0)).ravel()
    low_genes = adata.var_names[gene_cell_counts < 10]
    print("Number of genes expressed in <10 cells:", len(low_genes))
    print(low_genes[:20])

    gene_cell_counts = np.asarray((adata.X > 0).sum(axis=0)).ravel()

    plt.figure()
    plt.hist(gene_cell_counts, bins=50)
    plt.xlabel("Number of cells expressing a gene")
    plt.ylabel("Number of genes")
    plt.title(f"Distribution of Cells per Gene {name}")
    plt.show()

    plt.hist(gene_cell_counts, bins=100, range=(0, 200))
    plt.xlabel("Number of cells expressing a gene")
    plt.ylabel("Number of genes")
    plt.title(f"Distribution of Low Expressed Genes Across Cells {name}")
    plt.show()

    plt.figure()
    plt.violinplot(gene_cell_counts, showmedians=True)
    plt.ylabel("Number of cells expressing a gene")
    plt.title("Gene Detection Frequency Across Cells (Shannon)")
    plt.xticks([])
    plt.show()

    # Low-expression genes (<200 cells)
    low_counts = gene_cell_counts[gene_cell_counts < 200]

    plt.figure()
    plt.violinplot(low_counts, showmedians=True)
    plt.ylabel("Number of cells expressing a gene")
    plt.title("Low Gene Detection Frequency (<200 Cells) (Shannon)")
    plt.xticks([])
    plt.show()

    print(adata)
    print("Min cells per gene:", gene_cell_counts.min())
    print("Median cells per gene:", np.median(gene_cell_counts))
    print("Max cells per gene:", gene_cell_counts.max())

def genes_2_ens_id(adata, name):
    handle_duplicate_gene_names(adata)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", None)

    print(check_expresison(adata, name))

    print("mapping")

    mg = mygene.MyGeneInfo()

    genes = adata.var_names.tolist() #[:40]
    res = mg.querymany(genes,
                       scopes=["symbol", "alias", "name", "retired"],
                       fields=["symbol", "ensembl.gene", "entrezgene", "type_of_gene", "alias"],
                       species="human",
                       as_dataframe=True,
                       verbose=False)

    df_2 = res.reset_index()
    df_2 = df_2[(df_2["notfound"] == True)]
    print(df_2)
    print(f"Not mapped Genes {len(df_2)} out of {len(adata.var_names)} genes")


    df = res.reset_index()
    df = df[(df["symbol"] == df["query"]) & (df["ensembl.gene"].notna())]
    mapping = dict(zip(df["query"], df["ensembl.gene"]))
    print(f"Mapped genes {len(mapping.keys())}")

    new_names = [mapping.get(g, g) for g in genes]

    for old, new in list(zip(genes, new_names))[:20]:
        if old != new:
            print(old, "->", new)

    adata.var_names = new_names
    print(len(new_names))
    print(len(set(new_names)))
    new_adata = handle_duplicate_gene_names(adata)

    return new_adata