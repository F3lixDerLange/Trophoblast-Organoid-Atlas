import matplotlib.pyplot as plt
import seaborn as sns

figsize = (10,8)

def similarity_df_heatmap(result_df, title):
    plt.figure(figsize=figsize)

    ax = sns.heatmap(
        result_df,
        annot=True,  # change to True if you want numbers inside
        cmap="viridis",
        cbar_kws={"label": "Cosine similarity"},
    )

    if title:
        ax.set_title(title, fontsize=14)

    ax.set_xlabel("target cell types")
    ax.set_ylabel("reference cell types")

    plt.tight_layout()
    plt.show()
