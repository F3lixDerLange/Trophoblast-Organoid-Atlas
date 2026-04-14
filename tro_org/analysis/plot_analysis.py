import matplotlib.pyplot as plt
import seaborn as sns

figsize = (12,11)
fontsize = 20

def similarity_df_heatmap(result_df, title, method, savedir, ident):
    plt.figure(figsize=figsize)


    ax = sns.heatmap(
        result_df,
        annot=True,  # change to True if you want numbers inside
        cmap="magma",
        annot_kws={"size": fontsize-4},
        cbar_kws={"label": f"{method}"},
    )

    if title:
        ax.set_title(title, fontsize=fontsize+5)

    ax.set_xlabel("target cell types", fontsize=fontsize+2)
    ax.set_ylabel("reference cell types", fontsize=fontsize+2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=45, ha="right", fontsize=fontsize)
    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel(f"{method}", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize-2)

    plt.tight_layout()
    plt.savefig(f"{savedir}/{ident}_{method}_hvg.png")
    plt.show()


def common_gene_heatmap(df, title):
    plt.figure(figsize=(int(1.25*df.shape[1]),10))
    ax = sns.heatmap(
        df,
        annot=True,  # change to True if you want numbers inside
        cmap="viridis",
        cbar_kws={"label": "expression mean"},
    )

    if title:
        ax.set_title(title, fontsize=14)

    ax.set_xlabel("common genes")
    ax.set_ylabel("dataset")

    plt.tight_layout()
    plt.savefig(f"plots/{title}.png")
    plt.show()