import argparse
import scanpy as sc
try:
    import pyliger as pl
except ImportError:
    print("install pyliger first on Linux")



def pyl_integration(h5ad_file):
    adata = sc.read_h5ad(h5ad_file)
    ifnb_liger = pl.create_liger(adata)
    pl.select_genes(ifnb_liger)
    pl.scale_not_center(ifnb_liger)
    pl.optimize_ALS(ifnb_liger, k=25)
    pl.quantile_norm(ifnb_liger)
    pl.leiden_cluster(ifnb_liger, resolution=0.25)
    pl.run_umap(ifnb_liger, distance='cosine', n_neighbors=30, min_dist=0.3)
    all_plots = pl.plot_by_dataset_and_cluster(ifnb_liger, axis_labels=['UMAP 1', 'UMAP 2'], return_plots=True)
    print(all_plots)
    print(f"ifnb: {ifnb_liger}")
    print(f"adata: {adata}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True)
    # parser.add_argument("-output", required=True)
    args = parser.parse_args()
    input_file = args.input
    # out_dir = args.output
    pyl_integration(input_file)

if __name__ == '__main__':
    main()