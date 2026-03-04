import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text



def adjpvalue_hist(ds, plot_dir):
    plt.figure(figsize=(7, 4))
    plt.hist(ds["padj"].clip(0, 1), bins=50, color="#383a6b")
    plt.xlabel("Adjusted p-value (padj)")
    plt.ylabel("Number of genes")
    plt.title("Adjusted p-value distribution")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/adjpvalue_hist_organoid_vs_inVivo.png", dpi=300)
    plt.show()


def vulcano_plot(ds, lfc_thr, adjp_thr, plot_dir):
    df = ds[np.isfinite(ds["log2FoldChange"]) & np.isfinite(ds["padj"])]

    padj_safe = np.clip(df["padj"].to_numpy(), np.nextafter(0, 1), 1.0)
    neglog10 = -np.log10(padj_safe)

    df["_neglog10padj"] = neglog10

    filtered_ds = (df["padj"] < adjp_thr) & (np.abs(df["log2FoldChange"]) >= lfc_thr)

    plt.figure(figsize=(8, 6))
    plt.scatter(df.loc[~filtered_ds, "log2FoldChange"],
                df.loc[~filtered_ds, "_neglog10padj"],
                s=8, alpha=0.4, color="#383a6b")

    plt.scatter(df.loc[filtered_ds, "log2FoldChange"],
                df.loc[filtered_ds, "_neglog10padj"],
                s=8, alpha=0.8, color="#cb1f73")

    # Threshold lines
    plt.axhline(-np.log10(adjp_thr), linestyle="--", color="#fcc72d")
    plt.axvline(-lfc_thr, linestyle="--", color="#fcc72d")
    plt.axvline(+lfc_thr, linestyle="--", color="#fcc72d")

    n_each_side = 8
    n_top_pval = 20
    label_only_sig = True

    label_pool = df.loc[filtered_ds].copy() if label_only_sig else df.copy()
    if label_pool.shape[0] == 0:
        label_pool = df.copy()

    top_up = label_pool.sort_values("log2FoldChange", ascending=False).head(n_each_side)
    top_down = label_pool.sort_values("log2FoldChange", ascending=True).head(n_each_side)
    top_pval = label_pool.sort_values("_neglog10padj", ascending=False).head(n_top_pval)

    to_label = set(top_up.index) | set(top_down.index) | set(top_pval.index)
    texts = []

    for gene in to_label:
        x = float(df.loc[gene, "log2FoldChange"])
        y = float(df.loc[gene, "_neglog10padj"])

        texts.append(
            plt.text(x, y, str(gene), fontsize=7)
        )

    adjust_text(
        texts,
        ax=plt.gca(),
        expand_text=(1.2, 1.2),
        expand_points=(1.2, 1.2),
        force_text=(0.5, 0.5),
        force_points=(0.2, 0.2),
        arrowprops=dict(arrowstyle='-',
                        color='gray',
                        lw=0.5,
                        alpha=0.7,
                        shrinkA=10,
                        shrinkB=5)
    )

    plt.xlabel("log2 fold change (organoid vs in Vivo)")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano plot organoid vs in Vivo")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/volcano_plot_organoid_vs_inVivo.png", dpi=300)
    plt.show()