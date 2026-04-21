from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import tro_org.trajectory_analysis.plot_utils as pu
import tro_org.utils.utils as utils


def sc_paga(adata, emb_key, startcluster, cluster_key, output, n_neighbors=30):
    sc.set_figure_params(dpi_save=300, figsize=(10, 8), fontsize=14, vector_friendly=True)
    save_dir = Path(f"{output}/PAGA")
    sc.settings.figdir = save_dir
    utils.ensure_dir(f"{output}/PAGA")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=emb_key)
    sc.tl.umap(adata)

    sc.tl.diffmap(adata)

    sc.tl.paga(adata, groups=cluster_key)

    sc.set_figure_params(figsize=(10, 10))
    fig = sc.pl.paga(
        adata,
        color=cluster_key,
        #groups="label",
        layout="fr",
        threshold=0.05,  # tune this (0.03–0.1 typical)
        node_size_scale=3,
        node_size_power=0.7,
        edge_width_scale=2,
        fontsize=20,
        frameon=True,
        show=False
    )

    ax = plt.gca()
    for text in ax.texts:
        text.set_color("black")
    for line in ax.collections:
        if isinstance(line, LineCollection):
            line.set_color("gray")

    ax.set_title("PAGA connectivity of trophoblast organoid atlas", fontsize=18, pad=12)
    #ax.set_xlabel("PAGA dimension 1", fontsize=18)
    #ax.set_ylabel("PAGA dimension 2", fontsize=18)

    plt.tight_layout()
    plt.savefig(f"{output}/PAGA/paga_{emb_key}.png")

    sc.pl.umap(adata, color=[cluster_key], edges=True, save=f"_with_edges_{emb_key}.png")

    sc.tl.umap(adata, init_pos="paga") # https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-case_JUPYTER-trajectories/tutorial.html#re-draw-force-directed-graph
    sc.pl.umap(adata, color=[cluster_key, "batch"], wspace=0.4, save=f"_cluster_batch_{emb_key}.png")

    sc.pl.paga_compare(
        adata, threshold=0.03, title='Comparison', right_margin=0.2, size=10, edge_width_scale=0.5,
        legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=f"_comparison_{emb_key}.png")

    adata.uns['iroot'] = np.flatnonzero(adata.obs[cluster_key] == startcluster)[0]
    sc.tl.dpt(adata)
    sc.pl.umap(adata,
               color=['dpt_pseudotime', cluster_key],
               cmap="viridis",
               #legend_loc='on data',
               save=f"_dpt_pseudotime_{emb_key}.png"
               )

    #Debugging and decision making
    if "Cctb_2" in adata.obs[cluster_key].unique():
        ax = sc.pl.umap(adata, show=False)
        sc.pl.umap(adata[adata.obs[cluster_key] == "Cctb_2"].copy(),
                   color=cluster_key,
                   ax=ax,
                   )

    pu.plot_scatter_cluster_pseudotime(adata, 'dpt_pseudotime', output, "PAGA", emb_key, cluster_key)