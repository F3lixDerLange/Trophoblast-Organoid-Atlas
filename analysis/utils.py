import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics.pairwise import cosine_similarity
import os
import yaml


def load_config(config_file):
    if config_file.endswith(".yaml") or config_file.endswith(".yml"):
        if yaml is None:
            raise ImportError("Config file must be YAML file. Install with: pip install pyyaml")
        with open(config_file) as f:
            cfg = yaml.safe_load(f)
    return cfg["datasets"]

def load_data(config):
    datasets_loaded = []

    for dataset in config:
        name = dataset["name"]
        path = dataset["path"]
        label_col = dataset["label_col"]

        print(f"Loading {name} from {path}")
        adata = sc.read_h5ad(path)

        datasets_loaded.append({
            "name": name,
            "label_col": label_col,
            "adata": adata,
        })

    return datasets_loaded

def common_genes(dataset_ref, dataset_query):
    ref_name = dataset_ref["name"]
    query_name =dataset_query["name"]
    ref_adata = dataset_ref["adata"]
    query_adata = dataset_query["adata"]

    common_genes_list = ref_adata.var_names.intersection(query_adata.var_names)
    if len(common_genes_list) == 0:
        print(f"WARNING: no common genes between {dataset_ref["name"]} and {dataset_query["name"]}")

    print(f"Number of common genes: {len(common_genes_list)}")

    ref_adata = ref_adata[:, common_genes_list].copy()
    query_adata = query_adata[:, common_genes_list].copy()
    return ref_adata, query_adata