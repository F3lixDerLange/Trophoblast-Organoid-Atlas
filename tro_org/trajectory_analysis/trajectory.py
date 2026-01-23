import argparse
import scanpy as sc
from tro_org.trajectory_analysis import PAGA
from tro_org.trajectory_analysis import phlowerpy
from tro_org.trajectory_analysis import palantir_py
from tro_org.trajectory_analysis import scFates_py

def analysis(adata_file, output_dir, label_key, emb_key, startcluster):
    adata = sc.read_h5ad(adata_file)
    print(adata)
    PAGA.sc_paga(adata, emb_key, startcluster, label_key, output_dir)
    phlowerpy.phlower_traj(adata, emb_key, startcluster, label_key, output_dir)
    palantir_py.palantir_traj(adata, emb_key, startcluster, label_key,output_dir)
    scFates_py.scfates_traj(adata, emb_key, startcluster, label_key, output_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-adata", required=True, help="h5ad file")
    parser.add_argument("-output", required=False)
    parser.add_argument("-lk", "--label_key", required=True, help="label_key")
    args = parser.parse_args()
    adata_file = args.adata
    out_dir = args.output
    label_key = args.label_key

    analysis(adata_file, out_dir, label_key, "X_scvi", "Proliferative")

    """
    -adata
    /Users/felixlang/Downloads/Shibata_samplesheet_scvi/finalized/merged.h5ad
    -output
    ztest_folder
    -lk 
    celltype
    """

if __name__ == '__main__':
    main()