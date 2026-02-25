import scanpy as sc

def divide_dataset(adata):
    adata_in_vivo = adata[adata.obs["Model"] == "in vivo"].copy()
    adata_TBsOrg = adata[adata.obs["Model"] == "TBsOrg"].copy()
    adata_TBpOrg = adata[adata.obs["Model"] == "TBpOrg"].copy()

    adata_in_vivo.write_h5ad("database/Shannon_McNeil/divided_files/shannon_in_vivo_raw_filter.h5ad")
    adata_TBsOrg.write_h5ad("database/Shannon_McNeil/divided_files/shannon_TBsOrg_raw_filter.h5ad")
    adata_TBpOrg.write_h5ad("database/Shannon_McNeil/divided_files/shannon_TBpOrg_raw_filter.h5ad")


def main():
    adata_path = "database/Shannon_McNeil/Seurat/shannon_trophoblast_raw_filter.h5ad"
    adata = sc.read(adata_path)
    divide_dataset(adata)

if __name__ == '__main__':
    main()