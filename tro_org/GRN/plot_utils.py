import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import seaborn as sns
import loompy as lp
import scanpy as sc
import anndata as ad
import tro_org.utils.utils as utils
from tro_org.analysis.plot_analysis import fontsize

COLORS = ["#fcc72d", "#383a6b", "#cb1f73", "#e03a3c", "#ea6d3d", "#6a5fa8",
            "#f89c1c", "#b33a2b", "#7a1e3a", "#1f2a44", "#5c8d89", "#8a7bd1",
              "#2b7f7a", "#ffe07a", "#3fa7a3","#f05a28", "#ff6f61", "#d64aa0",
              "#2f2f5f", "#ffb347", "#8e2c2c", "#4b1630", "#0f1629"]

def quantile_histplot(adata, dataset, outdir):
    utils.ensure_dir(f"{outdir}/{dataset}")

    n_genes_detected_per_cell = np.sum(adata.X > 0, axis=1)
    percentiles = pd.Series(n_genes_detected_per_cell.flatten().A.flatten()).quantile(
        [0.01, 0.05, 0.10, 0.50, 1]
    )
    print(percentiles)

    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
    sns.histplot(n_genes_detected_per_cell, palette=["#383a6b"], bins='fd', ax=ax, legend=False, alpha=1.0)
    for i, x in enumerate(percentiles):
        fig.gca().axvline(x=x, ymin=0, ymax=1, color='red')
        ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i] * 100}%)', color='red', rotation=30,
                size='x-small', rotation_mode='anchor')
    ax.set_xlabel('# of genes')
    ax.set_ylabel('# of cells')
    fig.suptitle(f"Distribution of expressed genes per cell in {dataset} (adata.X > 0)", y=0.95)
    fig.tight_layout()
    plt.savefig(f"{outdir}/{dataset}/{dataset}_quantile_histplot.png", dpi=150)
    plt.show()

def pyscenic_heatmaps(adata, data_dir, dataset, outdir, label_key):
    sc.set_figure_params(dpi_save=300, figsize=(10, 8), fontsize=14, vector_friendly=True)
    utils.ensure_dir(f"{outdir}/{dataset}")
    sc.settings.figdir = f"{outdir}/{dataset}"


    lf = lp.connect(f"{data_dir}/{dataset}_pyscenic_output.loom", mode="r+", validate=False)
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    lf.close()

    ad_auc_mtx = ad.AnnData(auc_mtx)
    sc.pp.neighbors(ad_auc_mtx, n_neighbors=10, metric="correlation")
    sc.tl.umap(ad_auc_mtx)

    adata.obsm["X_umap_aucell"] = ad_auc_mtx.obsm["X_umap"]

    sc.pl.embedding(adata, basis="X_umap_aucell", color=label_key)

    auc_mtx[label_key] = adata.obs[label_key]
    mean_auc_by_cell_type = auc_mtx.groupby(label_key).mean()

    top_n = 50
    top_tfs = mean_auc_by_cell_type.max(axis=0).sort_values(ascending=False).head(top_n)
    mean_auc_by_cell_type_top_n = mean_auc_by_cell_type[
        [c for c in mean_auc_by_cell_type.columns if c in top_tfs]
    ]

    #spezi_colors = ["#fcc72d", "#ea6d3d", "#e03a3c", "#cb1f73", "#383a6b"]
    #spezi_cmap = LinearSegmentedColormap.from_list("spezi", spezi_colors)

    fig = sns.clustermap(
        mean_auc_by_cell_type_top_n,
        figsize=(15, 6.5),
        cmap=LinearSegmentedColormap.from_list("single_color",["#ffffff", "#cb1f73"]),
        xticklabels=True,
        yticklabels=True,
    )
    fig.figure.suptitle(f"Top {top_n} tf-regulons per cell-type", y=1.0)
    plt.savefig(f"{outdir}/{dataset}/{dataset}_tf_regulon_celltpye_heatmap.png", dpi=150)
    plt.show()

    tf_names = top_tfs.index.str.replace(r"\(\+\)", "", regex=True)
    adata_batch_top_tfs = adata[:, adata.var_names.isin(tf_names)]

    sc.pl.matrixplot(
        adata_batch_top_tfs,
        tf_names,
        groupby=label_key,
        cmap=LinearSegmentedColormap.from_list("single_color",["#ffffff", "#383a6b"]),
        dendrogram=False,
        figsize=(15, 5.5),
        standard_scale="group",
        title=f"Top {top_n} TFs associated to cell types",
        save=f"{dataset}_tf_regulon_celltpye.png",
    )
    plt.show()

    auc_threshold = 0.1
    auc_mtx_top = auc_mtx[mean_auc_by_cell_type_top_n.columns]
    auc_mtx_top = auc_mtx_top.loc[adata.obs_names]
    auc_mtx_binary = (auc_mtx_top > auc_threshold).astype(np.float32)
    adata_auc_cells = ad.AnnData(
        X=auc_mtx_binary.values,
        obs=adata.obs[[label_key]].copy(),
        var=pd.DataFrame(index=auc_mtx_top.columns)
    )
    adata_auc_cells.layers["aucell"] = auc_mtx_top.values.astype(np.float32)

    sc.pl.dotplot(
        adata_auc_cells,
        var_names=mean_auc_by_cell_type_top_n.columns.tolist(),
        groupby=label_key,
        layer="aucell",
        cmap=LinearSegmentedColormap.from_list("single_color", ["#ffffff", "#cb1f73"]),
        standard_scale="var",  # normalise per TF/regulon, not per group
        dendrogram=False,
        figsize=(15, 5.5),
        title=f"Top {top_n} TF regulons per cell type (AUCell activity)",
        colorbar_title="Mean AUCell\nactivity (scaled)",
        save=f"{dataset}_tf_regulon_celltype_dotplot.png",
    )

    safe_dataset = dataset.replace(".", "_")
    save_path = f"{outdir}/{safe_dataset}"
    utils.ensure_dir(save_path)
    top_tf_per_celltype_scatter(adata_batch_top_tfs, safe_dataset, save_path, label_key)
    top_tf_per_celltype_scatter_allinone(adata_batch_top_tfs, safe_dataset, save_path, label_key)


def top_tf_per_celltype_scatter(adata_batch_top_tfs, dataset, outdir, label_key):
    celltypes = adata_batch_top_tfs.obs[label_key].astype(str).unique()

    for ct in celltypes:
        expr_df = adata_batch_top_tfs.to_df()
        mean_expr = expr_df.mean(axis=0).sort_values(ascending=False)
        tf_order = mean_expr.sort_values(ascending=True)
        tf_order_len = np.arange(len(tf_order))

        fig, ax = plt.subplots(figsize=(5, 10), dpi=150)
        ax.scatter(tf_order.values, tf_order_len, color="#383a6b")
        ax.set_yticks(tf_order_len)
        ax.set_yticklabels(tf_order.index)
        ax.set_xlabel("Mean TF expression")
        ax.set_title(f"Top 50 mean expressed TFs in \n{ct} - {dataset}")
        ax.grid(True, axis="x", alpha=0.3)
        fig.tight_layout()
        plt.savefig(f"{outdir}/{dataset}_tf_mean_expression_in_{ct}", dpi=150)
        plt.show()

def top_tf_per_celltype_scatter_allinone(adata_batch_top_tfs, dataset, outdir, label_key):
    expr_df = adata_batch_top_tfs.to_df()
    expr_df[label_key] = adata_batch_top_tfs.obs[label_key].values
    print(expr_df.head())

    mean_expr = expr_df.groupby(label_key).mean()
    tf_order = mean_expr.mean(axis=0).sort_values(ascending=False).index
    mean_expr = mean_expr[tf_order]

    plot_df = (mean_expr.reset_index().melt(id_vars=label_key, var_name="TF", value_name="mean_expression"))

    plt.figure(figsize=(3,12), dpi=150)
    sns.scatterplot(data=plot_df, x="mean_expression", y="TF", hue=label_key, s=90, palette=COLORS)
    plt.xlabel("Mean TF expression", fontsize=18)
    plt.ylabel("TF", fontsize=18)
    plt.title(f"TF expression across cell types \n{dataset}", fontsize=20)
    plt.legend(title="label", loc='upper left', bbox_to_anchor=(1, 1), fontsize=14)
    plt.grid(True, axis="x", alpha=0.3)
    #plt.tight_layout()
    #save_path = f"{outdir}/{dataset}/{dataset}_tf_mean_expression_celltype_allinone".replace(".", "_")
    plt.savefig(f"{outdir}/{dataset}_tf_mean_expression_celltype_allinone.png", dpi=150, bbox_inches="tight")
    plt.show()


def tf_target_importance(adjacencies, dataset, outdir):
    utils.ensure_dir(f"{outdir}/{dataset}")
    plt.figure(figsize=(7, 5), dpi=150)
    plt.hist(np.log10(adjacencies["importance"]), bins=100)
    plt.xlim([-10, 10])
    plt.xlabel("log10(importance)")
    plt.ylabel("frequency")
    plt.title(f"Distribution TF-target importance")
    plt.savefig(f"{outdir}/{dataset}/{dataset}_tf_target_importance.png", dpi=150)
    plt.show()
