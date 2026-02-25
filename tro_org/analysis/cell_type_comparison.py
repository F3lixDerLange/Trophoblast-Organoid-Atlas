from typing import Literal

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from tro_org.analysis import utils
import plot_analysis
import scanpy as sc

def compute_cosign_similarity(datasets, savedir):
    results = {}

    print("Running cosine comparison")

    for i, ref_ds in enumerate(datasets):
        for j, target_ds in enumerate(datasets):
            if i <= j:
                continue  # avoid duplicates

            ref_ds_name = ref_ds["name"]
            target_ds_name = target_ds["name"]
            ref_label_col = ref_ds["label_col"]
            target_label_col = target_ds["label_col"]

            print(f"\n--- Map {ref_ds_name} -> {target_ds_name} ---")

            ref_common_genes, target_common_genes = utils.common_genes(ref_ds, target_ds)

            ref_cell_profile = compute_celltype_profiles(ref_common_genes, ref_label_col)
            target_cell_profile = compute_celltype_profiles(target_common_genes, target_label_col)

            cosine_similarity = pairwise_cosine_similarity(ref_cell_profile, target_cell_profile)
            results[(ref_ds_name, target_ds_name)] = cosine_similarity

            print(cosine_similarity.shape)
            plot_analysis.similarity_df_heatmap(cosine_similarity, f"ref:{ref_ds_name}_target:{target_ds_name}", "Cosine Similarity", savedir, f"{ref_ds_name}_{target_ds_name}")


def pairwise_cosine_similarity(profiles1, profiles2):
    sim = cosine_similarity(profiles1.values, profiles2.values)
    return pd.DataFrame(sim, index=profiles1.index, columns=profiles2.index)

def compute_celltype_profiles(adata, label_col):   # mean expression vectors of each gene per celltype
    labels = adata.obs[label_col].unique()
    X = adata.X.A if hasattr(adata.X, "A") else adata.X  # ensure dense or csr

    profiles = []
    celltypes = []

    for ct in labels:
        idx = np.where(adata.obs[label_col].values == ct)[0]
        mean_profile = X[idx].mean(axis=0)
        mean_profile = np.asarray(mean_profile).ravel()
        profiles.append(mean_profile)
        celltypes.append(ct)

    return pd.DataFrame(profiles, index=celltypes)

def jaccard_similarity(datasets, method, top_n, savedir):
    for i, ref_ds in enumerate(datasets):
        for j, target_ds in enumerate(datasets):
            if i <= j:
                continue  # avoid duplicates

            ref_ds_name = ref_ds["name"]
            target_ds_name = target_ds["name"]
            ref_label_col = ref_ds["label_col"]
            target_label_col = target_ds["label_col"]
            ref_adata = ref_ds["adata"]
            target_adata = target_ds["adata"]

            # Jaccard with marker genes: Do the same biological cell identities show the same differential signatures
            # Jaccard with hvg: Do the datasets have similar variance structure / technical + biological heterogeneity
            markerA = generate_marker_genes(ref_adata, ref_label_col, method, top_n)
            markerB = generate_marker_genes(target_adata, target_label_col, method, top_n)

            print(f"Computing Jaccard similarity -> {ref_ds_name} -> {target_ds_name}")
            jaccard_scors = calculate_jaccard_score(markerA, markerB)

            plot_analysis.similarity_df_heatmap(jaccard_scors, f"ref:{ref_ds_name}_target:{target_ds_name}, gene select: {method}", "Jaccard_Similarity", savedir, f"{ref_ds_name}_{target_ds_name}_methods:{method}")




def calculate_jaccard_score(markerA, markerB):
    cell_types_A = list(markerA.keys())
    cell_types_B = list(markerB.keys())

    jaccard_matrix = pd.DataFrame(
        np.zeros((len(cell_types_A), len(cell_types_B))),
        index=cell_types_A,
        columns=cell_types_B
    )

    for typ_a in cell_types_A:
        for typ_b in cell_types_B:
            genesA = markerA[typ_a]
            genesB = markerB[typ_b]
            union = genesA | genesB
            intersection = genesA & genesB
            jaccard_index = len(intersection) / len(union) if len(union) > 0 else 0
            jaccard_matrix.loc[typ_a, typ_b] = jaccard_index

    return jaccard_matrix


def generate_marker_genes(adata, label_col, methods: Literal["wilcoxon","t-test", None] = None, top_n=200, pval_cutoff=0.05):
    markers = {}

    if methods is None:
        for ct in adata.obs[label_col].unique():
            subset = adata[adata.obs[label_col] == ct]
            mean_expr = np.asarray(subset.X.mean(axis=0)).ravel()
            print(mean_expr)
            gene_names = np.array(adata.var_names)
            top_genes = gene_names[np.argsort(mean_expr)[-top_n:]]
            markers[ct] = set(top_genes)

    elif methods == "wilcoxon" or methods == "t-test":
        print("warning: normalized data must be in adata.X")
        adata.layers["norm"] = adata.X.copy()
        sc.tl.rank_genes_groups(adata,
                                method=methods,
                                groupby=label_col,
                                layer="norm",
                                use_raw=False,
                                key_added=str(methods)
        )
        df = sc.get.rank_genes_groups_df(adata, group=None, key=str(methods))

        print(df)

        for ct, subset in df.groupby("group"):
            sub = subset[
                (subset["logfoldchanges"] > 1.0) &
                (subset["pvals_adj"] < pval_cutoff)
                ]
            top = sub.sort_values("scores", ascending=False)["names"].head(top_n)
            print(ct, top)
            markers[ct] = set(top)

    return markers


def main():
    save_dir = "tro_org/analysis/analysis_plots"
    config = utils.load_config("tro_org/analysis/cosine_comp_config.yaml")
    print("Loaded datasets:", [d["name"] for d in config])
    hvg = True
    datasets = utils.load_data(config, hvg=hvg)
    if hvg:
        save_dir = f"{save_dir}_hvg"
    print(datasets)
    compute_cosign_similarity(datasets, save_dir)
    method: Literal["wilcoxon", "t-test", None] = "t-test"
    top_n = 400
    jaccard_similarity(datasets, method, top_n, save_dir)


if __name__ == '__main__':
    main()

