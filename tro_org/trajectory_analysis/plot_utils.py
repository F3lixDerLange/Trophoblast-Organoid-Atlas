import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def plot_scatter_cluster_pseudotime(adata, pseudotime_key, output_dir, method, emb_key):
    df = adata.obs[["label", pseudotime_key]].copy()
    df = df.dropna()

    order = df.groupby("label")[pseudotime_key].median().sort_values().index.tolist()
    df["label"] = pd.Categorical(df["label"], categories=order, ordered=True)

    colors = ["#3fa7a3", "#fcc72d", "#ea6d3d", "#e03a3c", "#cb1f73", "#6a5fa8",
              "#383a6b", "#f89c1c", "#b33a2b", "#7a1e3a", "#1f2a44", "#5c8d89",
              "#2b7f7a", "#ffe07a", "#f05a28", "#ff6f61", "#d64aa0", "#8a7bd1",
              "#2f2f5f", "#ffb347", "#8e2c2c", "#4b1630", "#0f1629"]

    color_map = dict(zip(order, colors))

    print(color_map)

    if len(color_map.keys()) <= 12:
        plt.figure(figsize=(10, 5))
    else:
        plt.figure(figsize=(10, 10))

    for lab, sub in df.groupby("label", observed=True):
        y = sub["label"].cat.codes.values
        y_jit = y + np.random.uniform(-0.15, 0.15, size=len(sub))
        plt.scatter(sub[pseudotime_key], y_jit, s=6, color=color_map[lab], label=str(lab))
    plt.yticks(range(len(order)), order)
    plt.xlabel(f"{emb_key} pseudotime")
    plt.ylabel("label")
    plt.title(f"Cells ordered by {emb_key} pseudotime")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{method}/scatter_cluster_pseudotime_{emb_key}.png")
    plt.show()