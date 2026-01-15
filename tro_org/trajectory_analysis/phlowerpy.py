import phlower
from matplotlib import pyplot as plt
from collections import Counter




# Followed https://phlower.readthedocs.io/en/latest/notebooks/fib2neuron.html and
# https://phlower.readthedocs.io/en/latest/notebooks/kidney.html

def phlower_traj(adata):
    print(adata.obs[["clusters", "label"]])
    cluster_label_map = (
        adata.obs
        .groupby("clusters")["label"]
        .agg(lambda x: x.value_counts().idxmax())
    )

    print(cluster_label_map)

    # load kidney anndata with MOJITOO reduction and clustering
    phlower.ext.ddhodge(adata, basis="X_scvi", roots=(adata.obs.label=="Lgr5"), k=7, npc=100, ndc=40, s=2,
                        lstsq_method='cholesky', verbose=True)

    figs = []
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    phlower.pl.nxdraw_group(adata, group_name="label",node_size=5, show_edges=False, label=False, ax=ax)
    plt.show()
    phlower.pl.nxdraw_score(adata, color='u', node_size=10)
    plt.show()

    # Delaunay triangulation to construct graph with holes
    phlower.tl.construct_delaunay(adata, cluster_name='label', node_attr='u', start_n=10, end_n=10, circle_quant=0.1,
                                  calc_layout=True)

    print(adata)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3), constrained_layout=True)
    phlower.pl.nxdraw_group(adata, graph_name="X_pca_ddhodge_g_triangulation_circle", node_size=5, show_edges=True,
                            show_legend=True, label=False, ax=ax[0])
    phlower.pl.nxdraw_score(adata, graph_name="X_pca_ddhodge_g_triangulation_circle", node_size=5, ax=ax[1],
                            colorbar=True)
    plt.show()

    fig, ax = plt.subplots(1, 2, figsize=(10, 3), constrained_layout=True)
    phlower.pl.nxdraw_group(adata, graph_name='X_pca_ddhodge_g_triangulation_circle',
                            layout_name='X_pca_ddhodge_g_triangulation_circle', node_size=5, show_edges=True,
                            labelstyle='text', labelsize=8, show_legend=True, ax=ax[0])
    phlower.pl.nxdraw_score(adata, graph_name='X_pca_ddhodge_g_triangulation_circle',
                            layout_name='X_pca_ddhodge_g_triangulation_circle', colorbar=True, node_size=5, label=False,
                            ax=ax[1])

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    phlower.pl.plot_triangle_density(adata, "X_pca_ddhodge_g_triangulation_circle",
                                     "X_pca_ddhodge_g_triangulation_circle", colorbar=True, edge_color='gray', ax=ax,
                                     node_size=5)

    # Graph holdge laplacian
    phlower.tl.L1Norm_decomp(adata)
    figs = []

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    phlower.pl.nxdraw_holes(adata, vector_dim=0, arrows=False, width=4, node_size=2, ax=ax, colorbar=True)

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    phlower.pl.nxdraw_holes(adata, vector_dim=1, arrows=False, width=4, node_size=2, ax=ax, colorbar=True)

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    phlower.pl.plot_eigen_line(adata, n_eig=8, linewidth='2', markersize=8, show_legend=False, ax=ax)

    phlower.tl.knee_eigen(adata)

    print("-----------------")
    print(adata)
    print("-----------------")

    # Preference random walk
    phlower.tl.random_climb_knn(adata, n=10000)
    fig_width = 9
    fig, axs = plt.subplots(1, 2, figsize=(8, 3))
    phlower.pl.plot_traj(adata, trajectory=adata.uns['knn_trajs'][0], colorid=0, node_size=1, ax=axs[0])
    phlower.pl.plot_traj(adata, layout_name="X_pca_ddhodge_g_triangulation_circle",
                         trajectory=adata.uns['knn_trajs'][0], colorid=0, node_size=1, ax=axs[1])


    # Project trajectories to harmonic space
    phlower.tl.trajs_matrix(adata)

    # dbscan clustering on the culmulative trajectory space
    phlower.tl.trajs_clustering(adata, eps=0.5)
    Counter(adata.uns['trajs_clusters'])

    anno_dic = {
        0: "Nh",
        1: "Lgr5",
        2: "Surface",
        3: "Stromal",
        4: "Glandular",
        5: "Pv_like",
        6: "Luminal",
        7: "Proliferative",
        8: "Surface",
        9: "Ciliated",
        10: "10",
        11: "Endo",
        12: "Nh",
        13: "13",
        14: "Ciliated"
    }


    adata.uns['annotation'] = [anno_dic.get(int(i), i) for i in adata.uns['trajs_clusters']]
    phlower.tl.harmonic_stream_tree(adata,
                                trajs_clusters='annotation',
                                retain_clusters=list(set(anno_dic.values())),
                                min_bin_number=20,
                                cut_threshold=1.5,
                                verbose=True)
    adata.obs['group_str'] = [str(i) for i in adata.obs['label']]
    phlower.ext.plot_stream_sc(adata, fig_size=(8, 5), color=['group_str'], show_legend=False, dist_scale=1, s=10)
    plt.show()

