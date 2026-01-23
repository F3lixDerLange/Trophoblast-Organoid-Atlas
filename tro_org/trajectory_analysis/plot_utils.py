import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def plot_scatter_cluster_pseudotime(adata, pseudotime_key, output_dir, method, emb_key):
    df = adata.obs[["label", pseudotime_key]].copy()
    df = df.dropna()

    order = df.groupby("label")[pseudotime_key].median().sort_values().index.tolist()
    df["label"] = pd.Categorical(df["label"], categories=order, ordered=True)

    colors = ["#3fa7a3", "#fcc72d", "#ea6d3d", "#e03a3c", "#cb1f73", "#6a5fa8",
              "#383a6b", "#f89c1c", "#b33a2b", "#7a1e3a", "#1f2a44", "#5c8d89"]

    color_map = dict(zip(order, colors))

    plt.figure(figsize=(10, 5))
    for lab, sub in df.groupby("label", observed=True):
        y = sub["label"].cat.codes.values
        y_jit = y + np.random.uniform(-0.15, 0.15, size=len(sub))
        plt.scatter(sub[pseudotime_key], y_jit, s=6, color=color_map[lab], label=str(lab))
    plt.yticks(range(len(order)), order)
    plt.xlabel("Palantir pseudotime")
    plt.ylabel("label")
    plt.title("Cells ordered by Palantir pseudotime")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{method}/scatter_cluster_pseudotime_{emb_key}.png")
    plt.show()