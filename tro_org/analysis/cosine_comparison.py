import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from tro_org.analysis import utils
import plot_analysis

def compute_cosign_similarity(datasets):
    results = {}

    print("Running scmap comparison")

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
            plot_analysis.similarity_df_heatmap(cosine_similarity, f"ref:{ref_ds_name}_target:{target_ds_name}")


def pairwise_cosine_similarity(profiles1, profiles2):
    sim = cosine_similarity(profiles1.values, profiles2.values)
    return pd.DataFrame(sim, index=profiles1.index, columns=profiles2.index)

def compute_celltype_profiles(adata, label_col):   # mean expression vectors of each gene per cell
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


def main():
    config = utils.load_config("analysis/cosine_comp_config.yaml")
    print("Loaded datasets:", [d["name"] for d in config])
    datasets = utils.load_data(config)
    compute_cosign_similarity(datasets)


if __name__ == '__main__':
    main()

