import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection


def sc_paga(adata, output_dir, label_key):
    embedding = "X_scvi" if "X_scvi" in adata.obsm else "X_pca"

    sc.pp.neighbors(adata, n_neighbors=30, use_rep=embedding)
    sc.tl.umap(adata)

    sc.tl.diffmap(adata)

    sc.tl.paga(adata, groups=label_key)

    sc.set_figure_params(figsize=(15, 10))
    fig = sc.pl.paga(
        adata,
        color=label_key,
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
    plt.show()


    sc.pl.umap(adata, color=[label_key], edges=True)

    sc.tl.umap(adata, init_pos="paga") # https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-case_JUPYTER-trajectories/tutorial.html#re-draw-force-directed-graph
    sc.pl.umap(adata, color=[label_key, "batch"], wspace=0.4)

    sc.pl.paga_compare(
        adata, threshold=0.03, title='comparison', right_margin=0.2, size=10, edge_width_scale=0.5,
        legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)

    adata.uns['iroot'] = np.flatnonzero(adata.obs['label'] == 'Proliferative')[0]
    sc.tl.dpt(adata)
    sc.pl.umap(adata, color=['label', 'dpt_pseudotime'], legend_loc='on data')

    df = adata.obs[["label", "dpt_pseudotime"]].copy().dropna()

    order = df.groupby("label")["dpt_pseudotime"].median().sort_values().index.tolist()
    df["label"] = pd.Categorical(df["label"], categories=order, ordered=True)

    colors = ["#3fa7a3","#fcc72d","#ea6d3d","#e03a3c", "#cb1f73","#6a5fa8",
              "#383a6b","#f89c1c","#b33a2b","#7a1e3a","#1f2a44","#5c8d89"]

    color_map = dict(zip(order, colors))

    plt.figure(figsize=(10, 5))
    for lab, sub in df.groupby("label", observed=True):
        y = sub["label"].cat.codes.values
        y_jit = y + np.random.uniform(-0.15, 0.15, size=len(sub))
        plt.scatter(sub["dpt_pseudotime"], y_jit, s=6, color=color_map[lab], label=str(lab))

    plt.yticks(range(len(order)), order)
    plt.xlabel("DPT pseudotime")
    plt.ylabel("label")
    plt.title("Cells ordered by DPT pseudotime")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.show()