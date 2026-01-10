import matplotlib.pyplot as plt

def plot_dotplot_benchmark(benchmark, output_dir):
    df = benchmark.get_results(min_max_scale=False)
    df.transpose()

    methods = df.columns.drop("Metric Type", errors="ignore")

    batch_correction = df.loc["Batch correction", methods].astype(float)
    bio_conservation = df.loc["Bio conservation", methods].astype(float)

    plt.figure(figsize=(7, 6))
    plt.scatter(batch_correction, bio_conservation)

    for method in methods:
        plt.text(
            batch_correction[method],
            bio_conservation[method],
            method,
            fontsize=10,
            ha="left",
            va="bottom"
        )

    plt.xlabel("Batch correction (aggregate score)")
    plt.ylabel("Biological conservation (aggregate score)")
    plt.title("Integration methods: Batch correction vs Bio conservation")

    if "Unintegrated" in methods:
        plt.axhline(bio_conservation["Unintegrated"], linestyle="--", linewidth=1)
        plt.axvline(batch_correction["Unintegrated"], linestyle="--", linewidth=1)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/dotplot_benchmark.png")
    plt.show()

