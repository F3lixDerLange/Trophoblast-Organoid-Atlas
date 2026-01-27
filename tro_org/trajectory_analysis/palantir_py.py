from pathlib import Path

import numpy as np
import palantir
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
import tro_org.trajectory_analysis.plot_utils as pu
import tro_org.utils.utils as utils

# based on https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb
# and https://github.com/quadbio/scRNAseq_analysis_python_vignette/blob/8a5d02e9d651fb130bad50e6c96575ec9f66086f/Tutorial.md#From-the-raw-counts-to-a-UMAP

def palantir_traj(adata, emb_key, startcluster, cluster_key, output, n_neighbors=30):
    sc.set_figure_params(dpi_save=300, figsize=(10, 8), fontsize=14, vector_friendly=True)
    save_dir = Path(f"{output}/palantir")
    sc.settings.figdir = save_dir
    utils.ensure_dir(f"{output}/palantir")

    sc.pp.neighbors(adata, use_rep=emb_key, n_neighbors=n_neighbors)
    sc.tl.umap(adata, min_dist=0.3)
    adata.obsm[f"{emb_key}_umap"] = adata.obsm["X_umap"].copy()

    dm_res = palantir.utils.run_diffusion_maps(adata, pca_key=emb_key, n_components=30)

    palantir.plot.plot_diffusion_components(adata, embedding_basis=emb_key)
    plt.savefig(f"{output}/palantir/diffusion_components_{emb_key}.png")
    plt.show()

    ms_data = palantir.utils.determine_multiscale_space(adata, n_eigs=5)

    # if you have no biologically defined root TODO Check for reasonable approach with Luis
    """
    if "n_genes_by_counts" not in adata.obs:
        sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Quick heuristic: pick cell with highest n_genes_by_counts (CytoTRACE-like)
    start_cell = adata.obs["n_genes_by_counts"].idxmax()
    print("Start cell:", start_cell)
    """
    # If you have a biologically defined root
    start_cell = adata.obs_names[adata.obs[cluster_key] == startcluster][0]
    print("Start cell:", start_cell)

    terminal_states, excluded_boundaries = palantir.core.identify_terminal_states(ms_data=ms_data,
                                                                                  early_cell=start_cell,
                                                                                  n_jobs=1,
                                                                                  knn=20)

    print("Terminal states:", terminal_states)
    print("Excluded boundaries:", excluded_boundaries)

    try:
        palantir.plot.highlight_cells_on_umap(adata, terminal_states, embedding_basis=f"{emb_key}_umap")
        plt.savefig(f"{output}/palantir/terminal_states_highlight_cells_on_umap{emb_key}.png")
        plt.show()

        pr_res = palantir.core.run_palantir(adata,
                                            early_cell=start_cell,
                                            num_waypoints=500,
                                            terminal_states=terminal_states,
                                            n_jobs=1,
                                            knn=20,
                                            use_early_cell_as_start=False,

        )

        palantir.plot.plot_palantir_results(adata, s=3)
        plt.savefig(f"{output}/palantir/palantir_results_{emb_key}.png")
        plt.show()

        masks = palantir.presults.select_branch_cells(adata, q=.02, eps=.02)
        palantir.plot.plot_branch_selection(adata, s=1)
        plt.savefig(f"{output}/palantir/branch_selection_{emb_key}.png")
        plt.show()

        print(adata)

        sc.pl.umap(
            adata,
            color=["palantir_pseudotime", cluster_key],  # change key if needed
            cmap="viridis",
            legend_loc='on data',
            save=f"_palantir_pseudotime_{emb_key}.png"
        )

        pu.plot_scatter_cluster_pseudotime(adata, "palantir_pseudotime", output, "palantir", emb_key)

    except Exception as e:
        print(f"Something went wrong {e}")
        print(f"Error in  scFates {output}")