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
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import anndata as ad

# all TFs from https://resources.aertslab.org/cistarget/tf_lists/
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


def create_grn(adata, loom_path_scenic, data_dir, dataset, image=None, num_workers=20):

    if not os.path.exists(os.path.join(data_dir, f"{dataset}_filtered_scenic.loom")):
        adata2loom(adata, loom_path_scenic)
    else:
        print(f"Skip loom generation ---- {dataset}_filtered_scenic.loom already exists")

    run_dir = Path(data_dir).resolve()
    if image == "docker":
        cmd = ["docker","run","--rm","--platform","linux/amd64","-v",f"{run_dir}:/data","aertslab/pyscenic:0.12.1"]

    elif image == "singularity":
        cmd = ["singularity", "run", "-B", f"{run_dir}:/data", f"{run_dir}/aertslab-pyscenic-0.12.1.sif"]

    else:
        raise SystemExit(f"Unknown image {image}")


    # STEP 1: Gene regulatory network inference, and generation of co-expression modules
    if not os.path.exists(os.path.join(data_dir, f"{dataset}_adj.tsv")):
        scenic_grn = [
            "pyscenic", "grn",
            f"/data/{dataset}_filtered_scenic.loom",
            "/data/allTFs_hg38.txt",
            "--method", "grnboost2",
            "--num_workers", f"{num_workers}",
            # "--sparse",
            "-o", f"/data/{dataset}_adj.tsv",
        ]

        print("Running pyscenic GRN")
        subprocess.run(cmd + scenic_grn, check=True)
        print("Done")
    else:
        print(f"Skip adj generation ---- {dataset}_adj.tsv already exists")

    adjacencies = pd.read_csv(f"{data_dir}/{dataset}_adj.tsv", index_col=False, sep='\t')
    print(adjacencies.tail())
    print(f"Number of associations: {adjacencies.shape[0]}")

    #plt.hist(np.log10(adjacencies["importance"]), bins=100)
    #plt.xlim([-10, 10])
    #plt.show()


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
            "--expression_mtx_fname", f"/data/{dataset}_filtered_scenic.loom",
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
    n_genes_detected_per_cell = np.sum(adata.X > 0, axis=1)
    percentiles = pd.Series(n_genes_detected_per_cell.flatten().A.flatten()).quantile(
        [0.01, 0.05, 0.10, 0.50, 1]
    )
    print(percentiles)

    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
    sns.histplot(n_genes_detected_per_cell, bins='fd', ax=ax, legend=False)
    for i, x in enumerate(percentiles):
        fig.gca().axvline(x=x, ymin=0, ymax=1, color='red')
        ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i] * 100}%)', color='red', rotation=30,
                size='x-small', rotation_mode='anchor')
    ax.set_xlabel('# of genes')
    ax.set_ylabel('# of cells')
    fig.suptitle(f"Distribution of genes per cell in {dataset}", y=0.95)
    fig.tight_layout()
    plt.show()

    if not os.path.exists(os.path.join(data_dir, f"{dataset}_pyscenic_output.loom")):
        scenic_aucell = [
            "pyscenic", "aucell",
            f"/data/{dataset}_filtered_scenic.loom",
            f"/data/{dataset}_reg.csv",
            "--output", f"/data/{dataset}_pyscenic_output.loom",
            "--num_workers", f"{num_workers}",
        ]

        print("Running pyscenic AUCELL")
        subprocess.run(cmd + scenic_aucell, check=True)
        print("Done")
    else:
        print(f"Skip ctx generation ---- {dataset}_pyscenic_output.loom already exists")

    lf = lp.connect(f"{data_dir}/{dataset}_pyscenic_output.loom", mode="r+", validate=False)
    auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    lf.close()

    ad_auc_mtx = ad.AnnData(auc_mtx)
    sc.pp.neighbors(ad_auc_mtx, n_neighbors=10, metric="correlation")
    sc.tl.umap(ad_auc_mtx)

    adata.obsm["X_umap_aucell"] = ad_auc_mtx.obsm["X_umap"]

    sc.pl.embedding(adata, basis="X_umap_aucell", color="label")

    auc_mtx["label"] = adata.obs["label"]
    mean_auc_by_cell_type = auc_mtx.groupby("label").mean()

    top_n = 50
    top_tfs = mean_auc_by_cell_type.max(axis=0).sort_values(ascending=False).head(top_n)
    mean_auc_by_cell_type_top_n = mean_auc_by_cell_type[
        [c for c in mean_auc_by_cell_type.columns if c in top_tfs]
    ]

    sns.clustermap(
        mean_auc_by_cell_type_top_n,
        figsize=[15, 6.5],
        cmap="Blues",
        xticklabels=True,
        yticklabels=True,
    )
    plt.show()

    tf_names = top_tfs.index.str.replace(r"\(\+\)", "", regex=True)
    print(tf_names)
    adata_batch_top_tfs = adata[:, adata.var_names.isin(tf_names)]

    sc.pl.matrixplot(
        adata_batch_top_tfs,
        tf_names,
        groupby="label",
        cmap="Reds",
        dendrogram=False,
        figsize=[15, 5.5],
        standard_scale="group",
    )

    plt.show()


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