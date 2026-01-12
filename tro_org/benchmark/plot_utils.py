import matplotlib.pyplot as plt
import os

def plot_dotplot_benchmark(df, output_dir):
    methods = df.columns.drop("Metric Type", errors="ignore")
    batch_correction = df.loc["Batch correction", methods].astype(float)
    bio_conservation = df.loc["Bio conservation", methods].astype(float)

    plt.figure(figsize=(8, 6))
    plt.scatter(batch_correction, bio_conservation, c=["#cb1f73"])

    for method in methods:
        plt.text(
            batch_correction[method],
            bio_conservation[method],
            method,
            fontsize=10,
            ha="left",
            va="bottom"
        )

    basename = os.path.basename(output_dir)

    plt.xlabel("Batch correction (aggregate score)")
    plt.ylabel("Biological conservation (aggregate score)")
    plt.title(f"{basename} Benchmark: Batch correction vs Bio conservation")

    if "Unintegrated" in methods:
        plt.axhline(bio_conservation["Unintegrated"], linestyle="--", linewidth=1, color="#383a6b")
        plt.axvline(batch_correction["Unintegrated"], linestyle="--", linewidth=1, color="#383a6b")

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{basename}_dotplot_benchmark.png")
    plt.show()

