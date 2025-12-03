import argparse
import os

import networkx as nx
import scanpy as sc
import scglue as scg

def scg_integration(adata, out_dir, batch_key, label_key):

    sc.pp.highly_variable_genes(adata, batch_key=batch_key) # , n_top_genes=4000, flavor="seurat_v3")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)

    print(adata)
    print(batch_key)

    scg.models.configure_dataset(
        adata,
        prob_model="Normal",
        use_highly_variable=True,
        use_rep="X_pca",
        use_batch=str(batch_key),
        #use_cell_type=label_key
    )

    genes = adata.var_names

    G = nx.Graph()
    for gene in genes:
        G.add_node(gene)
        G.add_edge(gene, gene, weight=1, sign=1)

    print(G)

    if os.path.exists(f"glue_batch/pretrain/pretrain.dill"):
        print("loading pretrain model")
        glue = scg.models.load_model("glue_batch/pretrain/pretrain.dill")
    else:
        print("pretrain model not found")
        glue =scg.models.fit_SCGLUE(
            {"rna": adata},
            G,
            skip_balance=True,
            fit_kws={
                "directory": "glue_batch",
                "max_epochs": 150,
            }
        )
        glue.save("glue_batch/pretrain/pretrain.dill")

    adata.obsm["X_scglue"] = glue.encode_data("rna", adata)

    sc.pp.neighbors(adata, use_rep="X_scglue")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="leiden_scglue", flavor="igraph", n_iterations=2)

    sc.pl.umap(adata, color=[label_key,batch_key ,"leiden_scglue"], wspace=0.4)

    print(adata)
    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True, help="h5ad file")
    parser.add_argument("-output", required=True)
    parser.add_argument("-bk", "--batch_key", required=True, help="batch_key")
    parser.add_argument("-lk", "--label_key", required=True, help="label_key")
    args = parser.parse_args()
    input_file = args.input
    out_dir = args.output
    batch_key = args.batch_key
    label_key = args.label_key
    scg_integration(sc.read_h5ad(input_file), out_dir, batch_key, label_key)

if __name__ == '__main__':
    main()