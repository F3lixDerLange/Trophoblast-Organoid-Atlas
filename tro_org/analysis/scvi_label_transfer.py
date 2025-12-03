import os.path
import scanpy as sc
import scvi
import pandas as pd
from matplotlib import pyplot as plt


def label_transfer(ref_file, target_file):
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    ref = load_h5ad(ref_file)
    target = load_h5ad(target_file)
    # print(ref.obs)
    target.obs_names_make_unique()
    # print(target.obs_names)
    target = target.concatenate(ref)
    print(target.obs)

    target.layers["counts"] = target.X.copy()
    target.raw = target
    sc.pp.highly_variable_genes(target, flavor='seurat_v3', n_top_genes= 2000, layer= "counts", batch_key= "batch", subset= True)

    if os.path.exists("scvi_model/scvi_model"):
        print("loading scvi model")
        model = scvi.model.SCVI.load("scvi_model/scvi_model", adata=target)
    else:
        scvi.model.SCVI.setup_anndata(target, layer='counts', batch_key='batch')
        model = scvi.model.SCVI(target)
        model.train(accelerator="gpu")
        model.save("scvi_model/scvi_model", overwrite=True)

    target.obs["celltype"] = target.obs["celltype"].cat.add_categories("unknown")
    target.obs = target.obs.fillna(value={"celltype" : "unknown"})


    if os.path.exists("scvi_model/scvi_label_model"):
        print("loading scvi label_model")
        label_model = scvi.model.SCANVI.load("scvi_model/scvi_label_model", adata=target)
    else:
        label_model = scvi.model.SCANVI.from_scvi_model(model, adata=target, unlabeled_category="unknown", labels_key="celltype")
        label_model.train(accelerator="gpu", max_epochs=20, n_samples_per_label=100)
        label_model.save("scvi_model/scvi_label_model", overwrite=True)

    target.obs["predicted"] = label_model.predict(target)
    print(target)

    # target.obs["bc2"] = target.obs.index.map(lambda x: x[:-2])
    # cell_mapper = dict(zip(target.obs.bc2, target.obs.predicted))

    sc.tl.pca(target, svd_solver="arpack")
    sc.pp.neighbors(target)
    sc.tl.umap(target)
    print(target)

    plt.rcParams['figure.figsize'] = (12, 6)
    for ann in ["cell_annotation", "predicted"]:
        adata_batch = target[target.obs["batch"] == "0"]

        sc.pl.umap(
            adata_batch,
            color=ann,
            title=f"UMAP – Batch {ann}",
            frameon=False,
            legend_loc="right margin"
        )

def load_h5ad(h5ad_file):
    h5ad = sc.read_h5ad(h5ad_file)
    return h5ad


def main():
    ref_file = "database/Shibata/GSE241052_ari_org_annotated_fixed_normalized.h5ad"
    target_file = "database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene_fixed.h5ad"
    label_transfer(ref_file, target_file)

if __name__ == '__main__':
    main()