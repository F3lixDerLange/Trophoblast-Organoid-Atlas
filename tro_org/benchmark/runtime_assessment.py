import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from matplotlib.pyplot import title

import utils

fontsize = 16

def parse_data(log_path):
    usage_data = {}
    log_files = utils.find_file_w_pattern(log_path, "*_usage.log")
    print("log_files:", log_files)

    for log_file in log_files:
        with open(log_file, "r") as file:
            next(file)
            tool = log_file.split("/")[-1].split("_")[0]
            print(tool)
            usage_data[tool] = []
            for line in file:
                usage_data[tool].append(line.strip().split(" "))
                # print(usage_data[tool])

    return usage_data

def calculate_metrics(usage_data, adata, plot_dir):
    FMT = "%Y-%m-%d %H:%M:%S"

    # bar plot
        # time
        # cpu max
        # RAM max
        # VRAM max

    # line plot
        # x = time in sec
        # y = RAM in MB

    max_usage = {}
    ram_usage = {}

    for tool in usage_data:
        time_tmp = [x[0] for x in usage_data[tool]]
        tdelta = datetime.strptime(utils.unit_conversion_time(max(time_tmp)), FMT) - datetime.strptime(utils.unit_conversion_time(min(time_tmp)), FMT)
        cpu_tmp = [float(x[2]) for x in usage_data[tool]]
        cpu_max = max(cpu_tmp)
        ram_tmp = [float(x[4]) / (1024 ** 2) for x in usage_data[tool]]
        ram_max = max(ram_tmp) # / 1e6  # MB = 1e3 but also change unit in plot
        vram_tmp = [float(x[5]) for x in usage_data[tool]]
        vram_max = max(vram_tmp) / (1024 ** 2)  # MB = 1e3 but also change unit in plot

        print(tool , tdelta)
        max_usage[tool] = {"time_min": tdelta.total_seconds() / 60, "cpu": cpu_max, "rss": ram_max, "vsz": vram_max}
        ram_usage[tool] = ram_tmp

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    resources_normalized_plot(max_usage, [n_cells, n_genes], plot_dir)
    resources_plot(max_usage, [n_cells, n_genes], plot_dir)
    mem_line_plot(ram_usage, [n_cells, n_genes], plot_dir)

def resources_normalized_plot(data_dir, adata_shape, plot_dir):
    methods = list(data_dir.keys())
    metrics = ["time_min", "cpu", "rss", "vsz"]
    metric_labels = ["Time", "max CPU", "max RSS", "max VSZ"]

    data = np.array([[data_dir[m][metric] for m in methods] for metric in metrics])
    data_norm = data / data.max(axis=1, keepdims=True)
    x = np.arange(len(metrics))
    width = 0.18

    plt.figure(figsize=(12, 6))

    for i, method in enumerate(methods):
        plt.bar(x + (i - 1.5) * width, data_norm[:, i], width, label=method)

    title = "Computational Requirement Metrics Normalized"

    if len(adata_shape) > 1:
        title = title + f"  for {adata_shape[0]} cells & {adata_shape[1]} genes"

    plt.xticks(x, metric_labels)
    plt.ylabel("Normalized value (0–1)")
    plt.title(title, fontsize=fontsize)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/computational_requirement_metrics_norm_comparison.png", dpi=300)
    plt.show()

def resources_plot(data_dir, adata_shape, plot_dir):
    dataset = {}
    for tool in data_dir:
        dataset[tool] = [x for x in data_dir[tool].values()]

    methods = list(dataset.keys())
    metrics = ["Time (min)", "CPU (%)", "Peak RSS (GB)", "Peak VMem (GB)"]

    # Convert dict → matrix for convenience
    data = np.array(list(dataset.values())).T

    fig, axes = plt.subplots(1, 4, figsize=(14, 5), sharey=False)

    for i, ax in enumerate(axes):
        ax.bar(methods, data[i], color="steelblue")
        ax.set_title(metrics[i], fontsize=12)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', linestyle='--', alpha=0.4)

    title = "Computational Requirement Metrics"

    if len(adata_shape) > 1:
        title = title + f" for {adata_shape[0]} cells & {adata_shape[1]} genes"

    fig.suptitle(title, fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/computational_requirement_metrics_comparison.png", dpi=300)
    plt.show()

def mem_line_plot(mem_data, adata_shape, plot_dir):
    plt.figure(figsize=(10, 6))

    colors = {
        "Harmony": "red",
        "BBKNN": "blue",
        "scVI": "green",
        "Combat": "purple",
    }

    for tool, ram_values in mem_data.items():
        ram_values = np.array(ram_values)
        x = np.arange(len(ram_values))
        color = colors.get(tool, None)
        plt.plot(x, ram_values, label=tool, color=color, linewidth=2)

    title = "RAM usage over time per integration tool"

    if len(adata_shape) > 1:
        title = title + f" for {adata_shape[0]} cells & {adata_shape[1]} genes"


    plt.xlabel("Time (s)")
    plt.ylabel("RAM usage (GB)")
    plt.title(title, fontsize=fontsize)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/computational_requirement_metrics_ram_usage.png", dpi=300)
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-h5ad", required=False, help="h5ad file")
    parser.add_argument("-log", required=True, help="resource_usage dir scdownstream")
    parser.add_argument("-out", required=True, help="plot dir")

    args = parser.parse_args()
    h5ad = args.h5ad
    resource_usage_dir = args.log
    output_dir = args.out

    data = parse_data(resource_usage_dir)
    calculate_metrics(data, sc.read_h5ad(h5ad), output_dir)

if __name__ == '__main__':
    main()