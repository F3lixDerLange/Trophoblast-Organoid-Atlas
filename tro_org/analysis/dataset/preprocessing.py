import argparse
import os
import scanpy as sc
import tro_org.analysis.dataset.map_genes as map_genes

def map_gene_symbols(adata_path):
    adata = sc.read_h5ad(adata_path)

    raw_count_adata = map_genes.raw_counts_2_layer(adata)
    mapped_adata = map_genes.genes_2_ens_id(raw_count_adata, os.path.splitext(os.path.basename(adata_path))[0])
    processed_adata = map_genes.filter_genes(mapped_adata)

    #processed_adata.write(
    #    f"/Users/felixlang/Documents/Uni/Master/master-thesis/database/prepros_test/{os.path.splitext(os.path.basename(adata_path))[0]}_preprocessed.h5ad")
    #processed_adata.write_h5ad(map_genes.generate_out(adata_path))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-adata", required=True, help="h5ad file")
    args = parser.parse_args()
    adata_file = args.adata

    map_gene_symbols(adata_file)

if __name__ == '__main__':
    main()




