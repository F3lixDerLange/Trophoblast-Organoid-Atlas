import argparse
import os

import networkx as nx
import scanpy as sc
import scglue as scg
import tro_org.utils.plot_utils as pu

def scg_integration(adata, out_dir, batch_key, label_key, modeldir, gtf):

    # sc.pp.highly_variable_genes(adata, batch_key=batch_key) # , n_top_genes=4000, flavor="seurat_v3")
    # adata = adata[:, adata.var["highly_variable"]].copy()
    # if adata.raw.X is not None:
    #     adata.layers["counts"] = adata.raw.X.copy()

    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=100, svd_solver='auto')

    pu.plot_umap_before_integration(adata, "scglue", out_dir, batch_key, label_key)

    scg.data.get_gene_annotation(
        adata, gtf=gtf,
        gtf_by="gene_name"
    )


    coord_cols = ["chrom", "chromStart", "chromEnd"]
    print("NaNs per coord column:")
    print(adata.var[coord_cols].isna().sum())
    mask = adata.var[coord_cols].notna().all(axis=1)
    print(f"Keeping {mask.sum()} / {adata.n_vars} genes with valid coordinates")
    adata = adata[:, mask].copy()

    guidance = scg.genomics.rna_anchored_guidance_graph(adata, adata, propagate_highly_variable=False)
    scg.graph.check_graph(guidance, adata)

    scg.models.configure_dataset(
        adata,
        prob_model="Normal",
        use_highly_variable=False,
        use_rep="X_pca",
        use_batch=str(batch_key),
        #use_cell_type=label_key
    )

    """
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
                "max_epochs": 200,
            }
        )
        glue.save("glue_batch/pretrain/pretrain.dill")
    """

    glue = scg.models.fit_SCGLUE(
        {"rna": adata},
        guidance,   #G
        skip_balance=True,
        fit_kws={
            "directory": f"glue_batch/{modeldir}",
            "max_epochs": 200,
        }
    )

    adata.obsm["X_scglue"] = glue.encode_data("rna", adata)
    pu.plot_umap_after_integration(adata,"X_scglue", "scglue", out_dir, batch_key, label_key)

    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True, help="h5ad file")
    parser.add_argument("-output", required=True)
    parser.add_argument("-bk", "--batch_key", required=True, help="batch_key")
    parser.add_argument("-lk", "--label_key", required=True, help="label_key")
    parser.add_argument("-gtf", required=True, help="path to gtf annotation file")
    args = parser.parse_args()
    input_file = args.input
    out_dir = args.output
    batch_key = args.batch_key
    label_key = args.label_key
    gtf_path = args.gtf
    base_name = os.path.basename(input_file)
    filename = os.path.splitext(base_name)[0]
    scg_integration(sc.read_h5ad(input_file), out_dir, batch_key, label_key, filename, gtf_path)

if __name__ == '__main__':
    main()