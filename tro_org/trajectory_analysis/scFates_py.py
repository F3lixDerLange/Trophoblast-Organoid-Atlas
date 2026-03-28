# based on https://scfates.readthedocs.io/en/latest/Tree_Analysis_Bone_marrow_fates.html
from pathlib import Path
import palantir
import scFates as scf
import scanpy as sc
from matplotlib import pyplot as plt

import tro_org.trajectory_analysis.plot_utils as pu
import tro_org.utils.utils as utils


def scfates_traj(adata, emb_key, startcluster, cluster_key, output, n_neighbors=30):
    sc.set_figure_params(dpi_save=300, figsize=(10, 8), fontsize=14, vector_friendly=True)
    save_dir = Path(f"{output}/scFates")
    sc.settings.figdir = save_dir
    utils.ensure_dir(f"{output}/scFates")
    sc.set_figure_params(figsize=(10, 10))

    print(f"Running scFates trajectory analysis")
    sc.pp.neighbors(adata, use_rep=emb_key)
    sc.tl.umap(adata, min_dist=0.3)
    adata.obsm[f"X_{emb_key}_umap"] = adata.obsm["X_umap"].copy() # TODO check if necessary

    # Run Palantir to obtain multiscale diffusion space
    print("diffusion")
    dm_res = palantir.utils.run_diffusion_maps(adata, pca_key=emb_key, n_components=30)
    ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=10)
    start_cell = adata.obs_names[adata.obs[cluster_key] == startcluster][0]
    print("Start cell:", start_cell)

    # Generate embedding from the multiscale diffusion space
    adata.obsm["X_palantir"] = ms_data.values
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_palantir")

    # Tree learning with SimplePPT
    scf.tl.tree(adata, method="ppt", Nodes=150, use_rep="palantir",
                device="cpu", seed=42, ppt_lambda=100, ppt_sigma=0.025, ppt_nsteps=200)

    # projecting results onto ForceAtlas2 embedding
    scf.pl.graph(adata, save=f"_scf_tree_{emb_key}.png") #, basis=f"{emb_key}_umap", save=f"_scf_tree_{emb_key}.png")
    # Selecting a root and computing pseudotime
    adata.obs["is_root_cell"] = adata.obs_names == start_cell

    scf.tl.cleanup(adata)

    try:
        scf.tl.root(adata,"is_root_cell")  # TODO get correct root

        scf.tl.pseudotime(adata, n_jobs=20, n_map=100, seed=42)
        scf.pl.trajectory(adata, save=f"_pseudotime_tree_{emb_key}.png")

        #as a dendrogram representation
        scf.tl.dendrogram(adata)

        scf.pl.dendrogram(adata, color="seg", save=f"_dendrogram_seg_{emb_key}.png")
        scf.pl.dendrogram(adata, color="t", show_info=False, cmap="viridis", save=f"_dendrogram_pseudotime_{emb_key}.png")
        scf.pl.dendrogram(adata,
                          color=cluster_key,
                          legend_loc="on data",
                          color_milestones=True,
                          legend_fontoutline=True,
                          save=f"_dendrogram_label_{emb_key}.png"
                          )

        sc.pl.umap(
            adata,
            color=["t", cluster_key],
            cmap="viridis",
            #legend_loc='on data',
            save=f"_scFates_pseudotime_{emb_key}.png"
        )

        # Debugging and decision making
        if "Cctb_2" in adata.obs[cluster_key].unique():
            ax = sc.pl.umap(adata, show=False)
            sc.pl.umap(adata[adata.obs[cluster_key] == "Cctb_2"].copy(),
                       color=cluster_key,
                       ax=ax,
                       )

        pu.plot_scatter_cluster_pseudotime(adata, 't', output, "scFates", emb_key)

    except Exception as e:
        print(f"Something went wrong {e}")
        print(f"Error in  scFates {output}")
