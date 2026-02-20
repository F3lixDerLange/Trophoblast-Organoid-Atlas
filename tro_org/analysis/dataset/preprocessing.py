import argparse
import os
import scanpy as sc
import tro_org.analysis.dataset.map_genes as map_genes
import tro_org.analysis.dataset.dataset_analysis as dataset_analysis
import tro_org.analysis.dataset.gene_venn as gene_venn

def map_gene_symbols(adata_path):
    adata = sc.read_h5ad(adata_path)

    raw_count_adata = map_genes.raw_counts_2_layer(adata)
    mapped_adata = map_genes.genes_2_ens_id(raw_count_adata, os.path.splitext(os.path.basename(adata_path))[0])
    processed_adata = map_genes.filter_genes(mapped_adata)


    # outpath = f"/Users/felixlang/Documents/Uni/Master/master-thesis/database/prepros_test/{os.path.splitext(os.path.basename(adata_path))[0]}_preprocessed.h5ad"
    outpath = map_genes.generate_out(adata_path)
    #processed_adata.write_h5ad(outpath)

    return outpath

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-adata", required=False, help="h5ad file")
    args = parser.parse_args()
    adata_file = args.adata

    if adata_file is None:
        h5ad_files = [["database/Shibata/Shibata_fixed_normalized.h5ad", "Shibata"],
                      ["database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene.h5ad", "Arutyunyan_PTO"],
                      ["database/Arutyunyan/Arutyunyan_TSC/Organoid_TSC_cellxgene.h5ad", "Arutyunyan_TSC"],
                      ["database/Shannon_McNeil/Seurat/shannon_trophoblast.h5ad", "Shannon"]]
        processed_adata_dict = {}
        for dataset in h5ad_files:
            print(dataset)
            dataset_name = dataset[1].split("_")[0]
            if dataset_name not in processed_adata_dict :
                processed_adata_dict [dataset_name] = [map_gene_symbols(dataset[0])]
            else:
                processed_adata_dict [dataset_name].append(map_gene_symbols(dataset[0]))

        print(processed_adata_dict)

        savedir = "figures"
        dataset_analysis.manage_data(h5ad_files, savedir)
        gene_venn.create_venn(dataset_analysis.merge_dataset_from_same_study(h5ad_files)) # turn list with multiple files form one study in dict #before
        gene_venn.create_venn(processed_adata_dict) # after
    else:
        map_gene_symbols(adata_file)


if __name__ == '__main__':
    main()




