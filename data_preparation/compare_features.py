import argparse
from xml.sax.handler import all_features

import scanpy as sc
import pandas as pd
import os
import sys


def read_docs(db_dir):
    all_features_per_study =[]
    for study in os.listdir(db_dir):
        study_folder = os.path.join(db_dir, study)
        if os.path.isfile(study_folder):
            continue
        for dataset in os.listdir(study_folder):
            dataset_folder = os.path.join(study_folder, dataset)
            if os.path.isfile(dataset_folder):
                continue
            for filename in os.listdir(dataset_folder):
                if filename.endswith(".h5ad") and "Patient" not in dataset_folder:
                    h5_data = sc.read_h5ad(os.path.join(dataset_folder, filename))
                    print(filename)
                    print(h5_data.shape)
                    gene_ids = h5_data.var['gene_ids-0'] if "cellxgene" in filename else h5_data.var['gene_ids']
                    gene_ids_list = gene_ids.tolist()
                    all_features_per_study.append(gene_ids_list)

    all_features_sets = [set(lst) for lst in all_features_per_study]

    core_genes = set.intersection(*all_features_sets)

    print(f"Core Genes in all sets: {len(core_genes)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study_dir", required=True)
    args = parser.parse_args()
    study_dir = args.study_dir
    read_docs(study_dir)


if __name__ == '__main__':
    main()