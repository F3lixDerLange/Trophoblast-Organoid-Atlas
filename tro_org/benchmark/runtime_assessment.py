import argparse
import os
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from tro_org.utils.utils import filter_invivo_cells
from matplotlib.pyplot import title

import tro_org.benchmark.utils as utils

fontsize = 16

def parse_data(log_path, tsv):
    usage_data = {}
    pattern = "*_usage.tsv" if tsv else "*_usage.log"
    log_files = utils.find_file_w_pattern(log_path, pattern)
    print("log_files:", log_files)

    for log_file in log_files:
        with open(log_file, "r") as file:
            next(file)
            tool = log_file.split("/")[-1].split("_")[0]
            # print(tool)
            usage_data[tool] = []
            for line in file:
                separator = "\t" if tsv else " "
                usage_data[tool].append(line.strip().split(separator))
                # print(usage_data[tool])

    return usage_data

def calculate_metrics(usage_data, tsv):
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

    if tsv:
        for tool in usage_data:
            time_tmp = [float(x[2]) for x in usage_data[tool]]
            tdelta = int(time_tmp[-1] + 1)
            ram_tmp = [float(x[3]) / 1024 for x in usage_data[tool]]
            ram_max = max(ram_tmp)
            vram_tmp = [float(x[4]) / 1024 for x in usage_data[tool]]
            vram_max = max(vram_tmp)
            cpu_tmp = [float(x[5]) for x in usage_data[tool]]
            cpu_max = max(cpu_tmp)

            max_usage[tool] = {"time_min": tdelta / 60, "cpu": cpu_max, "rss": ram_max, "vsz": vram_max}
            ram_usage[tool] = ram_tmp

            print(tool, tdelta)

    else:
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

    return max_usage, ram_usage


def resources_normalized_plot(data_dir, adata_shape, plot_dir, dataset_ident):
    methods = list(data_dir.keys())
    metrics = ["time_min", "cpu", "rss", "vsz"]
    metric_labels = ["Time", "max CPU", "max RSS", "max VSZ"]

    data = np.array([[data_dir[m][metric] for m in methods] for metric in metrics])
    data_norm = data / data.max(axis=1, keepdims=True)
    x = np.arange(len(metrics))
    n_methods = len(methods)
    group_width = 0.8
    width = group_width / n_methods

    plt.figure(figsize=(12, 6))

    for i, method in enumerate(methods):
        plt.bar(x + (i - (n_methods - 1) / 2) * width, data_norm[:, i], width, label=method, color="#383a6b")

    title = "Computational Resource Usage Metrics Normalized"

    if len(adata_shape) > 1:
        title = title + f"  for {adata_shape[0]} cells & {adata_shape[1]} genes"

    plt.xticks(x, metric_labels)
    plt.ylabel("Normalized value (0–1)")
    plt.title(title, fontsize=fontsize)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{dataset_ident}_computational_requirement_metrics_norm_comparison.png", dpi=300)
    plt.show()

def resources_plot(data_dir, adata_shape, plot_dir, dataset_ident):
    dataset = {}
    for tool in data_dir:
        dataset[tool] = [x for x in data_dir[tool].values()]

    methods = list(dataset.keys())
    metrics = ["Time (min)", "CPU (%)", "Peak RSS (GB)", "Peak VMem (GB)"]

    # Convert dict → matrix for convenience
    data = np.array(list(dataset.values())).T

    fig, axes = plt.subplots(1, 4, figsize=(14, 5), sharey=False)

    for i, ax in enumerate(axes):
        ax.bar(methods, data[i], color="#383a6b")
        ax.set_title(metrics[i], fontsize=12)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', linestyle='--', alpha=0.4)

    title = "Computational Resource Usage Metrics"

    if len(adata_shape) > 1:
        title = title + f" for {adata_shape[0]} cells & {adata_shape[1]} genes"

    fig.suptitle(title, fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{dataset_ident}_computational_requirement_metrics_comparison.png", dpi=300)
    plt.show()

def mem_line_plot(mem_data, adata_shape, plot_dir, dataset_ident):
    plt.figure(figsize=(10, 6))

    colors = {
        "harmony": "#3fa7a3",
        "bbknn": "#fcc72d",
        "scvi": "#ea6d3d",
        "combat": "#e03a3c",
        "liger": "#cb1f73",
        "scglue": "#6a5fa8",
        "scanorama": "#383a6b",
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
    plt.xscale("log")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{dataset_ident}_computational_requirement_metrics_ram_usage.png", dpi=300)
    plt.show()

def plot_metrics(max_usage, ram_usage, adata, plot_dir, dataset_ident):
    n_cells = adata.n_obs
    n_genes = adata.n_vars

    resources_normalized_plot(max_usage, [n_cells, n_genes], plot_dir, dataset_ident)
    resources_plot(max_usage, [n_cells, n_genes], plot_dir, dataset_ident)
    mem_line_plot(ram_usage, [n_cells, n_genes], plot_dir, dataset_ident)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-h5ad", required=False, help="h5ad file")
    parser.add_argument("-log", required=False, help="resource_usage dir scdownstream")
    parser.add_argument("-out", required=False, help="plot dir")
    parser.add_argument("-tsv", required=False, help="tsv dir")

    args = parser.parse_args()
    h5ad = args.h5ad
    resource_usage_dir = args.log
    output_dir = args.out
    tsv_path = args.tsv


    tsv_data = parse_data(tsv_path, True)
    log_data = parse_data(resource_usage_dir, False)
    log_usage, log_ram = calculate_metrics(log_data, False)
    tsv_usage, tsv_ram = calculate_metrics(tsv_data, True)


    if resource_usage_dir is None and tsv_path is None:
        print("You must specify a resource usage directory or a tsv file")
    elif resource_usage_dir is not None and tsv_path is None:
        max_usage, ram_usage = log_usage, log_ram
    elif resource_usage_dir is None and tsv_path is not None:
        max_usage, ram_usage = tsv_usage, tsv_ram
    elif resource_usage_dir is not None and tsv_path is not None:
        max_usage = tsv_usage | log_usage
        ram_usage = tsv_ram | log_ram

    dataset_identifier = os.path.basename(tsv_path)
    plot_metrics(max_usage, ram_usage, filter_invivo_cells(sc.read_h5ad(h5ad),"batch"), output_dir, dataset_identifier)

    """
    -h5ad
    /Users/felixlang/Downloads/pipline_fixed/annotated_samplesheet_scvi/finalized/merged.h5ad
    -log
    /Users/felixlang/Downloads/ann_integration_in_order
    -out
    tro_org/benchmark/benchmark_plots
    -tsv
    /Users/felixlang/Downloads/ann_merged_inOrder/resource_usage
    """

    """
    python3 -m tro_org.benchmark.runtime_assessment -h5ad ~/nf-core_out/merged_integration/merged_samplesheet_scvi/finalized/merged.h5ad -out ~/torg/runtime/final -log ~/nf-core_out/merged_integration/ -tsv ~/torg/benchmark/merged_integration_final/resource_usage    
    """

if __name__ == '__main__':
    main()