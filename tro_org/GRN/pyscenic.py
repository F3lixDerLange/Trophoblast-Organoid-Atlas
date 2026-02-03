# based on https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess

# all TFs from https://resources.aertslab.org/cistarget/tf_lists/

def adata2loom(adata, loom_path):
    row_attrs = {
        "Gene": np.array(adata.var_names),
    }
    col_attrs = {
        "CellID": np.array(adata.obs_names),
        "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
        "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
    }
    lp.create(loom_path, adata.X.transpose(), row_attrs, col_attrs)


def create_grn(adata, loom_path_scenic, data_dir, dataset, image=None):
    adata2loom(adata, loom_path_scenic)

    # STEP 1: Gene regulatory network inference, and generation of co-expression modules
    run_dir = Path(data_dir).resolve()
    scenic_grn = [
        "pyscenic", "grn",
        f"/data/{dataset}_filtered_scenic.loom",
        "/data/allTFs_hg38.txt",
        "--method", "grnboost2",
        "--num_workers", "20",
        # "--sparse",
        "-o", f"/data/{dataset}_adj.tsv",
    ]

    if image == "docker":
        img = "aertslab/pyscenic:0.12.1"
        cmd = [
                         "docker", "run", "--rm",
                         "--platform", "linux/amd64",
                         "-v", f"{run_dir}:/data",
                         img
                     ] + scenic_grn

    elif image == "singularity":
        cmd = [   "singularity", "run",
                              "-B", f"{run_dir}:/data",
                              f"{run_dir}/aertslab-pyscenic-0.12.1.sif",
                          ] + scenic_grn
    else:
        raise SystemExit(f"Unknown image {image}")

    print("Running pyscenic command")
    subprocess.run(cmd, check=True)
    print("Done")

    adjacencies = pd.read_csv("tro_org/GRN/data/adj.tsv", index_col=False, sep='\t')
    adjacencies.head()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--image", required=True, help="docker or sigularity image")
    parser.add_argument("-a", "--adata", required=True, help="h5ad file")
    parser.add_argument("-d", "--data_dir", required=True, help="data dir for docker")
    #parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    image = args.image
    data_path = args.data_dir
    adata_path = args.adata
    #out_dir = args.output

    #adata_path = "tro_org/trajectory_analysis/figures/adata/subcluster_Epithelial_integrated.h5ad"
    #data_path = "/Users/felixlang/Documents/Uni/Master/master-thesis/tro_org/GRN/data"
    # all_tfs = f"{adata_path}/allTFs_hg38.txt"
    dataset = "_".join(os.path.basename(adata_path).split("_")[:2])
    loom_path_scenic = f"{data_path}/{dataset}_filtered_scenic.loom"


    adata = sc.read_h5ad(adata_path)
    create_grn(adata, loom_path_scenic, data_path, dataset, image)



if __name__ == '__main__':
    main()