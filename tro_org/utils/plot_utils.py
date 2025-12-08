from pathlib import Path
import scanpy as sc
from matplotlib.pyplot import title


def plot_umap_before_integration(adata, method, outdir, batch_key, label_key):
    save_dir = Path(f"{outdir}/{method}")
    sc.settings.figdir = save_dir
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=[label_key, batch_key], wspace=0.4, save=f"_before_integration_{method}.png")
    print("before")
    print(adata)

def plot_umap_after_integration(adata, osbm_key, method, outdir, batch_key, label_key):
    save_dir = Path(f"{outdir}/{method}")
    sc.settings.figdir = save_dir
    sc.pp.neighbors(adata, use_rep=osbm_key)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added=f"leiden_{method}", flavor="igraph", n_iterations=2, resolution=0.4)
    sc.pl.umap(adata, color=[label_key, batch_key, f"leiden_{method}"], wspace=0.4, save=f"_after_integration_{method}.png")
    print("after")
    print(adata)

    print("highlight one batch")

    ax = sc.pl.umap(adata, show=False)
    for batch in adata.obs[batch_key].unique():
        sc.pl.umap(adata[adata.obs[batch_key] == batch],
                   color=label_key,
                   ax=ax,
                   title=f"{method} Umap {batch} batch highlighted",
                   save=f"_after_integration_{method}_{batch}_highlighted.png")