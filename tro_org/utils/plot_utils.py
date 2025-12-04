from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import os

def plot_umap_before_integration(adata, method, outdir, batch_key, label_key):
    save_dir = Path(f"{outdir}/{method}")
    sc.settings.figdir = save_dir
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=label_key,save=f"umap_after_integration_{method}.png")
    plt.savefig()
    print("before")
    print(adata)

def plot_umap_after_integration(adata, key, method, outdir, batch_key, label_key):
    save_dir = Path(f"{outdir}/{method}")
    sc.settings.figdir = save_dir
    sc.pp.neighbors(adata, use_rep=key)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added=f"leiden_{method}", flavor="igraph", n_iterations=2)
    sc.pl.umap(adata, color=[label_key, batch_key, "leiden_scglue"], wspace=0.4, save=f"umap_after_integration_{method}.png")
    print("after")
    print(adata)