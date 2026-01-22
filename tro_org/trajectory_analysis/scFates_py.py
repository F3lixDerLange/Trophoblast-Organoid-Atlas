# based on https://scfates.readthedocs.io/en/latest/Tree_Analysis_Bone_marrow_fates.html
import numpy as np
import palantir
import pandas as pd
import scFates as scf
import scanpy as sc
from matplotlib import pyplot as plt


def scfates_traj(adata, osbm_key):
    sc.set_figure_params(figsize=(10, 10))

    sc.pp.neighbors(adata, use_rep=osbm_key)
    sc.tl.umap(adata, min_dist=0.3)
    adata.obsm["X_scvi_umap"] = adata.obsm["X_umap"].copy() # TODO check if necessary

    # Run Palantir to obtain multiscale diffusion space
    dm_res = palantir.utils.run_diffusion_maps(adata, pca_key="X_scvi", n_components=30)
    ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=10)

    # Generate embedding from the multiscale diffusion space
    adata.obsm["X_palantir"] = ms_data.values
    sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_palantir")

    # Tree learning with SimplePPT
    scf.tl.tree(adata, method="ppt", Nodes=150, use_rep="palantir",
                device="cpu", seed=1, ppt_lambda=100, ppt_sigma=0.025, ppt_nsteps=200)

    # projecting results onto ForceAtlas2 embedding
    scf.pl.graph(adata, basis="scvi_umap")
    sc.pl.umap(adata, color="label", size=20, )

    # Selecting a root and computing pseudotime
    root_cell = "Shibata_AAACGCTAGCATTGAA-1_1" # from palantir
    adata.obs["is_root_cell"] = adata.obs_names == root_cell

    scf.tl.cleanup(adata)

    scf.tl.root(adata,"is_root_cell")  # TODO get correct root

    scf.tl.pseudotime(adata, n_jobs=20, n_map=100, seed=42)
    scf.pl.trajectory(adata)

    #as a dendrogram representation
    scf.tl.dendrogram(adata)
    scf.pl.dendrogram(adata, color="seg")
    scf.pl.dendrogram(adata, color="t", show_info=False, cmap="viridis")
    scf.pl.dendrogram(adata, color="label", legend_loc="on data", color_milestones=True, legend_fontoutline=True)

    df = adata.obs[["label", "t"]].copy()
    df = df.dropna()
    order = df.groupby("label")["t"].median().sort_values().index.tolist()
    df["label"] = pd.Categorical(df["label"], categories=order, ordered=True)

    colors = ["#3fa7a3", "#fcc72d", "#ea6d3d", "#e03a3c", "#cb1f73", "#6a5fa8",
              "#383a6b", "#f89c1c", "#b33a2b", "#7a1e3a", "#1f2a44", "#5c8d89"]

    color_map = dict(zip(order, colors))

    plt.figure(figsize=(10, 5))
    for lab, sub in df.groupby("label", observed=True):
        y = sub["label"].cat.codes.values
        y_jit = y + np.random.uniform(-0.15, 0.15, size=len(sub))
        plt.scatter(sub["t"], y_jit, s=6, color=color_map[lab], label=str(lab))
    plt.yticks(range(len(order)), order)
    plt.xlabel("scFates pseudotime")
    plt.ylabel("label")
    plt.title("Cells ordered by scFates pseudotime")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()
    plt.show()

