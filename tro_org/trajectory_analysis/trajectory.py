import argparse
import scanpy as sc
import anndata as ad
from tro_org.trajectory_analysis import PAGA
from tro_org.trajectory_analysis import phlowerpy
from tro_org.trajectory_analysis import palantir_py
from tro_org.trajectory_analysis import scFates_py
from tro_org.analysis.utils import load_config

def setup_traj_analysis(adata_file, output_dir, label_key, dataset=None, startcluster=None):

    if adata_file.endswith(".yaml") or adata_file.endswith(".yml"):
        datasets_config = load_config(adata_file)
        for dc in datasets_config:
            name = dc["name"]
            path = dc["path"]
            root = dc["root"]

            print(f"Loading {name} from {path}")
            adata = sc.read_h5ad(path)
            traj_analysis(adata, output_dir, label_key, root, name)

    elif isinstance(adata_file, ad.AnnData):
        print(adata_file)
        traj_analysis(adata_file, output_dir, label_key, startcluster, dataset)


def traj_analysis(adata, output_dir, label_key, startcluster, dataset):
    methods = ['BBKNN', 'Combat', 'Harmony', 'LIGER', 'Scanorama', 'Unintegrated', 'scGlue', 'scVI']

    for emb_key in methods:
        if emb_key in adata.obsm.keys():
            output_dir_tmp = f"{output_dir}/{dataset}/{emb_key}"
            print(output_dir_tmp)
            PAGA.sc_paga(adata, emb_key, startcluster, label_key, output_dir_tmp)
            phlowerpy.phlower_traj(adata, emb_key, startcluster, label_key, output_dir_tmp)
            palantir_py.palantir_traj(adata, emb_key, startcluster, label_key,output_dir_tmp)
            scFates_py.scfates_traj(adata, emb_key, startcluster, label_key, output_dir_tmp)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-adata", required=True, help="h5ad file")
    parser.add_argument("-output", required=False)
    parser.add_argument("-lk", "--label_key", required=True, help="label_key")
    # parser.add_argument("-d", "--dataset", required=True, help="name of dataset")
    args = parser.parse_args()
    adata_file = args.adata
    out_dir = args.output
    label_key = args.label_key
    # dataset = args.dataset

    # setup_traj_analysis(adata_file, out_dir, label_key, dataset, "Vct_p")
    setup_traj_analysis(adata_file, out_dir, label_key)

    """
    -adata tro_org/trajectory_analysis/trajectory_datasets.yaml 
    -output tro_org/trajectory_analysis/figures 
    -lk label 
    # -d shibata
    """

if __name__ == '__main__':
    main()