import numpy as np
import palantir
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt


# based on https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb
# and https://github.com/quadbio/scRNAseq_analysis_python_vignette/blob/8a5d02e9d651fb130bad50e6c96575ec9f66086f/Tutorial.md#From-the-raw-counts-to-a-UMAP

def palantir_traj(adata, emb_key):
    diffusion_map(adata, emb_key)

def diffusion_map(adata, emb_key, n_neighbors=30):
    sc.pp.neighbors(adata, use_rep=emb_key, n_neighbors=n_neighbors)
    sc.tl.umap(adata)

    dm_res = palantir.utils.run_diffusion_maps(adata, pca_key=emb_key, n_components=30)

    print(adata)

    palantir.plot.plot_diffusion_components(adata, embedding_basis=emb_key)
    plt.show()

    ms_data = palantir.utils.determine_multiscale_space(adata, n_eigs=5)

    # if you have no biologically defined root
    """
    if "n_genes_by_counts" not in adata.obs:
        sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Quick heuristic: pick cell with highest n_genes_by_counts (CytoTRACE-like)
    start_cell = adata.obs["n_genes_by_counts"].idxmax()
    print("Start cell:", start_cell)
"""
    # If you have a biologically defined root
    start_cell = adata.obs_names[adata.obs["label"] == "Proliferative"][0]
    print("Start cell:", start_cell)

    terminal_states, excluded_boundaries = palantir.core.identify_terminal_states(ms_data=ms_data,
                                                                                  early_cell=start_cell,
                                                                                  n_jobs=1,
                                                                                  knn=20)

    print("Terminal states:", terminal_states)
    print("Excluded boundaries:", excluded_boundaries)

    palantir.plot.highlight_cells_on_umap(adata, terminal_states, embedding_basis="X_scvi-global_umap")
    plt.show()

    pr_res = palantir.core.run_palantir(adata,
                                        start_cell,
                                        num_waypoints=500,
                                        terminal_states=terminal_states,
                                        n_jobs=1,
                                        knn=20,
                                        use_early_cell_as_start=False,

    )

    palantir.plot.plot_palantir_results(adata, s=3)
    plt.show()

    masks = palantir.presults.select_branch_cells(adata, q=.02, eps=.02)
    palantir.plot.plot_branch_selection(adata, s=1)
    plt.show()

    print(adata)

    sc.pl.umap(
        adata,
        color="palantir_pseudotime",  # change key if needed
        cmap="viridis",
        size=10
    )

    plt.show()

    df = adata.obs[["label", "palantir_pseudotime"]].copy()
    df = df.dropna()

    # 3) define order (optional): keep your desired order or sort by median pseudotime
    order = df.groupby("label")["palantir_pseudotime"].median().sort_values().index.tolist()
    df["label"] = pd.Categorical(df["label"], categories=order, ordered=True)

    colors = ["#3fa7a3", "#fcc72d", "#ea6d3d", "#e03a3c", "#cb1f73", "#6a5fa8",
              "#383a6b", "#f89c1c", "#b33a2b", "#7a1e3a", "#1f2a44", "#5c8d89"]

    color_map = dict(zip(order, colors))

    plt.figure(figsize=(10, 5))
    for lab, sub in df.groupby("label", observed=True):
        y = sub["label"].cat.codes.values
        y_jit = y + np.random.uniform(-0.15, 0.15, size=len(sub))
        plt.scatter(sub["palantir_pseudotime"], y_jit, s=6, color=color_map[lab], label=str(lab))
    plt.yticks(range(len(order)), order)
    plt.xlabel("Palantir pseudotime")
    plt.ylabel("label")
    plt.title("Cells ordered by Palantir pseudotime")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.show()
