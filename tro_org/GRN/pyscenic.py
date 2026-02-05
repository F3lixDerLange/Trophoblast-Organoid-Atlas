# based on https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb
# https://www.sc-best-practices.org/mechanisms/gene_regulatory_networks.html#preparation-of-scenic
import argparse
import os
from pathlib import Path
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import subprocess
import time
import tro_org.GRN.plot_utils as pu

# all files downloaded from https://resources.aertslab.org/cistarget/
# there is also an explanation which files to use

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
    print("Loom created")


def filter_adata(adata):
    # filtering steps from https://www.nature.com/articles/s41596-020-0336-2
    if "percent_mito" in adata.obs.columns:
        percent_mito = "percent_mito"
    elif "pct_counts_mt" in adata.obs.columns:
        percent_mito = "pct_counts_mt"
    else:
        raise LookupError

    print(f"Adata shape pre filtering {adata.shape}")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs['n_genes'] < 4000, :]
    adata = adata[adata.obs[percent_mito] < 0.15, :]
    print(f"Adata shape post filtering {adata.shape}")

    return adata


def create_grn(adata, data_dir, dataset, image=None, num_workers=8, adata_filter=False):
    scenic_dir = os.path.split(data_dir)[0]
    save_dir = Path(f"{scenic_dir}/figure")

    if adata_filter:
        adata = filter_adata(adata)
        dataset = f"{dataset}_filtered"

    loom_path_scenic = f"{data_dir}/{dataset}_scenic.loom"

    if not os.path.exists(os.path.join(data_dir, f"{dataset}_scenic.loom")):
        adata2loom(adata, loom_path_scenic)
    else:
        print(f"Skip loom generation ---- {dataset}_scenic.loom already exists")


    run_dir = Path(data_dir).resolve()
    if image == "docker":
        cmd = ["docker","run","--rm","--platform","linux/amd64","-v",f"{run_dir}:/data","aertslab/pyscenic:0.12.1"]

    elif image == "singularity":
        cmd = ["singularity", "run", "--cleanenv","-B", f"{run_dir}:/data", f"{run_dir}/aertslab-pyscenic-0.12.1.sif"]

    else:
        raise SystemExit(f"Unknown image {image}")


    # STEP 1: Gene regulatory network inference, and generation of co-expression modules
    if not os.path.exists(os.path.join(data_dir, f"{dataset}_adj.tsv")):
        scenic_grn = [
            "pyscenic", "grn",
            f"/data/{dataset}_scenic.loom",
            "/data/allTFs_hg38.txt",
            "--method", "grnboost2",
            "--num_workers", f"{num_workers}",
            "--sparse",
            "-o", f"/data/{dataset}_adj.tsv",
        ]

        print("Running pyscenic GRN")
        subprocess.run(cmd + scenic_grn, check=True)
        print("Done")
    else:
        print(f"Skip adj generation ---- {dataset}_adj.tsv already exists")

    adjacencies = pd.read_csv(f"{data_dir}/{dataset}_adj.tsv", index_col=False, sep='\t')
    print(adjacencies.tail())
    print(f"shape: {adjacencies.shape}")

    pu.tf_target_importance(adjacencies, dataset, save_dir)


    # STEP 2-3: Regulon prediction aka cisTarget from CLI
    f_db_glob = f"{data_dir}/*.feather"
    f_db_names = ' '.join([os.path.join("/data", os.path.basename(i)) for i in glob.glob(f_db_glob)])
    f_motif_path = "/data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

    if not os.path.exists(os.path.join(data_dir, f"{dataset}_reg.csv")):
        scenic_ctx = [
            "pyscenic", "ctx",
            f"/data/{dataset}_adj.tsv",
            f"{f_db_names}",
            f"--annotations_fname", f"{f_motif_path}",
            "--expression_mtx_fname", f"/data/{dataset}_scenic.loom",
            "--output", f"/data/{dataset}_reg.csv",
            "--mask_dropouts",
            "--num_workers", f"{num_workers}",
        ]

        print("Running pyscenic CTX")
        subprocess.run(cmd + scenic_ctx, check=True)
        print("Done")
    else:
        print(f"Skip ctx generation ---- {dataset}_reg.csv already exists")


    #STEP 4: Cellular enrichment (aka AUCell) from CLI
    pu.quantile_histplot(adata, dataset, save_dir)

    if not os.path.exists(os.path.join(data_dir, f"{dataset}_pyscenic_output.loom")):
        scenic_aucell = [
            "pyscenic", "aucell",
            f"/data/{dataset}_scenic.loom",
            f"/data/{dataset}_reg.csv",
            "--output", f"/data/{dataset}_pyscenic_output.loom",
            "--num_workers", f"{num_workers}",
        ]

        print("Running pyscenic AUCELL")
        subprocess.run(cmd + scenic_aucell, check=True)
        print("Done")
    else:
        print(f"Skip ctx generation ---- {dataset}_pyscenic_output.loom already exists")

    pu.pyscenic_heatmaps(adata, data_dir, dataset, save_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--image", required=True, help="docker or sigularity image")
    parser.add_argument("-a", "--adata", required=True, help="h5ad file")
    parser.add_argument("-d", "--data_dir", required=True, help="data dir for docker")
    parser.add_argument("-f", "--adata_filter", required=False, action="store_true", help="filter adata")
    #parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    image = args.image
    data_path = args.data_dir
    adata_path = args.adata
    adata_f = args.adata_filter
    #out_dir = args.output

    dataset = "_".join(os.path.basename(adata_path).split("_")[:2])

    start = time.time()
    adata = sc.read_h5ad(adata_path)
    create_grn(adata, data_path, dataset, image, adata_filter=adata_f)
    end = time.time()
    print(f"{dataset} pyscenic GRN took {end-start} seconds -- {(end-start)//60} minutes")

    """
    -i docker
    -a tro_org/trajectory_analysis/figures/adata/subcluster_Stromal_integrated.h5ad
    -d /Users/felixlang/Documents/Uni/Master/master-thesis/tro_org/GRN/data
    """


if __name__ == '__main__':
    main()