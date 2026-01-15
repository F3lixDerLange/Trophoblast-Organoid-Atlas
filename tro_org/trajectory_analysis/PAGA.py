import scanpy as sc

def sc_paga(adata, output_dir, label_key):
    embedding = "X_scvi" if "X_scvi" in adata.obsm else "X_pca"

    sc.pp.neighbors(adata, n_neighbors=30, use_rep=embedding)
    sc.tl.umap(adata)

    sc.tl.paga(adata, groups=label_key)
    sc.pl.paga(adata, threshold=0.03)
    sc.pl.umap(adata, color=[label_key], edges=True)
