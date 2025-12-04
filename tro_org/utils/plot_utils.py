from pathlib import Path
import scanpy as sc

def plot_umap_before_integration(adata, method, outdir, batch_key, label_key):
    save_dir = Path(f"{outdir}/{method}")
    sc.settings.figdir = save_dir
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=label_key, save=f"_after_integration_{method}.png")
    print("before")
    print(adata)

def plot_umap_after_integration(adata, osbm_key, method, outdir, batch_key, label_key):
    #save_dir = Path(f"{outdir}/{method}")
    #sc.settings.figdir = save_dir
    sc.pp.neighbors(adata, use_rep=osbm_key)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added=f"leiden_{method}", flavor="igraph", n_iterations=2)
    sc.pl.umap(adata, color=[label_key, batch_key, f"leiden_{method}"], wspace=0.4, save=f"_after_integration_{method}.png")
    print("after")
    print(adata)