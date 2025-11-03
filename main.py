import scanpy as sc

adata = sc.read_10x_mtx(
    "database/GSE216244_RAW/UPT",         # folder containing the three files
    var_names="gene_symbols",  # or 'gene_ids' depending on the file
    make_unique=True
)
print(adata)
