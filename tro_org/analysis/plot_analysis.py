import matplotlib.pyplot as plt
import seaborn as sns

figsize = (10,8)

def similarity_df_heatmap(result_df, title, method, savedir, ident):
    plt.figure(figsize=figsize)

    ax = sns.heatmap(
        result_df,
        annot=True,  # change to True if you want numbers inside
        cmap="viridis",
        cbar_kws={"label": f"{method}"},
    )

    if title:
        ax.set_title(title, fontsize=14)

    ax.set_xlabel("target cell types")
    ax.set_ylabel("reference cell types")

    plt.tight_layout()
    plt.savefig(f"{savedir}/{ident}_{method}.png")
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