import argparse
import os
from pathlib import Path

import scanpy as sc
from scib_metrics.benchmark import *
import tro_org.integration.scGLUE_integration
import tro_org.integration.py_liger_integration
import tro_org.integration.scanorama_integration
import tro_org.benchmark.utils as utils
import tro_org.utils.plot_utils as pu
from tro_org.benchmark.usage_profiler import profile_resources
from tro_org.benchmark.plot_utils import plot_dotplot_benchmark

def run_integrations(base_adata, output_dir, batchkey, labelkey, modeldir, usage, gtf):
    osbm_keys = ["Unintegrated"]
    print("Running integrations")

    if usage:
        print("scGlue integration")
        run_scglue = (profile_resources("scGlue", f"{output_dir}/resource_usage/scglue_usage.tsv")
                      (tro_org.integration.scGLUE_integration.scg_integration))
        scGlue_adata = run_scglue(base_adata.copy(), output_dir, batchkey, labelkey, modeldir, gtf)
        osbm_keys.append("scGlue")

        print("PyLiger integration")
        run_liger = (profile_resources("LIGER", f"{output_dir}/resource_usage/liger_usage.tsv")
                     (tro_org.integration.py_liger_integration.pyl_integration))
        liger_adata = run_liger(base_adata.copy(), output_dir, batchkey, labelkey)
        osbm_keys.append("LIGER")

        print("scanorama_integration")
        run_scanorama = (profile_resources("Scanorama",
                                          f"{output_dir}/resource_usage/scanorama_usage.tsv")
                         (tro_org.integration.scanorama_integration.scn_integration))
        scanorama_adata = run_scanorama(base_adata.copy(), output_dir, batchkey, labelkey)
        osbm_keys.append("Scanorama")

    else:
        print("scGlue integration")
        scGlue_adata = tro_org.integration.scGLUE_integration.scg_integration(base_adata.copy(), output_dir, batchkey, labelkey, modeldir)
        osbm_keys.append("scGlue")

        print("PyLiger integration")
        liger_adata = tro_org.integration.py_liger_integration.pyl_integration(base_adata.copy(), output_dir, batchkey, labelkey)
        osbm_keys.append("LIGER")

        print("scanorama_integration")
        scanorama_adata = tro_org.integration.scanorama_integration.scn_integration(base_adata.copy(), output_dir, batchkey, labelkey)
        osbm_keys.append("Scanorama")

    integrated_adata = base_adata.copy()

    integrated_adata.obsm["scGlue"] = scGlue_adata.obsm["X_scglue"]
    integrated_adata.obsm["Scanorama"] = scanorama_adata.obsm["Scanorama"]
    integrated_adata.obsm["LIGER"] = liger_adata.obsm["LIGER"]

    if "X_pca" not in integrated_adata.obsm:
        print("Perform PCA")
        sc.pp.pca(integrated_adata)

    integrated_adata.obsm["Unintegrated"] = integrated_adata.obsm["X_pca"]

    return osbm_keys, integrated_adata


def get_data_from_scdownstream_merged(merged_dir, integrated_adata, obsm_keys, output_dir, batchkey, labelkey):
    merged_files = utils.find_file(merged_dir, "merged.h5ad" )
    print(merged_files)
    for tool, file in merged_files:
        merged_scdownstream = sc.read_h5ad(file)

        if tool == "scvi":
            integrated_adata.obsm["scVI"] = merged_scdownstream.obsm["X_scvi"]
            obsm_keys.append("scVI")
            pu.plot_umap_before_integration(merged_scdownstream, "scVI", output_dir, batchkey, labelkey)
            pu.plot_umap_after_integration(integrated_adata, "scVI", "scVI", output_dir, batchkey, labelkey)
        elif tool == "harmony":
            integrated_adata.obsm["Harmony"] = merged_scdownstream.obsm["X_harmony"]
            obsm_keys.append("Harmony")
            pu.plot_umap_before_integration(merged_scdownstream, "Harmony", output_dir, batchkey, labelkey)
            pu.plot_umap_after_integration(integrated_adata, "Harmony", "Harmony", output_dir, batchkey, labelkey)
        elif tool == "combat":
            integrated_adata.obsm['Combat'] = merged_scdownstream.obsm['X_combat']
            obsm_keys.append("Combat")
            pu.plot_umap_before_integration(merged_scdownstream, "Combat", output_dir, batchkey, labelkey)
            pu.plot_umap_after_integration(integrated_adata, "Combat", "Combat", output_dir, batchkey, labelkey)
        elif tool == "bbknn":
            integrated_adata.obsm["BBKNN"] = merged_scdownstream.obsm["X_bbknn-global_umap"]
            obsm_keys.append("BBKNN")
            pu.plot_umap_before_integration(merged_scdownstream, "BBKNN", output_dir, batchkey, labelkey)
            pu.plot_umap_after_integration(integrated_adata, "BBKNN", "BBKNN", output_dir, batchkey, labelkey)

    return obsm_keys, integrated_adata


def benchmark(h5ad_file, output_dir, batch_key, label_key, modeldir, merged_adata, usage, gtf):

    osbm_keys, integrated_adata = run_integrations(h5ad_file, output_dir, batch_key, label_key, modeldir, usage, gtf)
    if merged_adata is not None:
        osbm_keys, integrated_adata = get_data_from_scdownstream_merged(merged_adata, integrated_adata, osbm_keys, output_dir, batch_key, label_key)

    integrated_adata.obsm["Unintegrated"] = integrated_adata.obsm["X_pca"]

    print("benchmark")
    bm = Benchmarker(
        integrated_adata,
        batch_key=batch_key,
        label_key=label_key,
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        embedding_obsm_keys=osbm_keys,
        n_jobs=6
    )
    bm.benchmark()
    print("plotting benchmark")
    bm.plot_results_table(save_dir=output_dir)
    bm.plot_results_table(min_max_scale=False, save_dir=output_dir)
    df = bm.get_results(min_max_scale=False)
    plot_dotplot_benchmark(df.transpose(), output_dir)

    path = Path(output_dir)
    dataset = path.name
    integrated_adata.copy().write(f"{output_dir}/{dataset}_integrated.h5ad")

    print("Success")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", required=True, help="h5ad file")
    parser.add_argument("-output", required=True)
    parser.add_argument("-bk", "--batch_key", required=True, help="batch_key")
    parser.add_argument("-lk", "--label_key", required=True, help="label_key")
    parser.add_argument("-m", "--mergeddir", required=False, help="dir with merged scdownstream h5ad file")
    parser.add_argument("-u", "--usage", required=False, action="store_true" ,help="log resource usage")
    parser.add_argument("-gtf", required=True, help="path to gtf annotation file")
    args = parser.parse_args()
    input_file = args.input
    out_dir = args.output
    batch_key = args.batch_key
    label_key = args.label_key
    merged_dir = args.mergeddir
    usage = args.usage
    gtf_path = args.gtf
    base_name = os.path.basename(input_file)
    filename = os.path.splitext(base_name)[0]

    benchmark(sc.read_h5ad(input_file), out_dir, batch_key, label_key, filename, merged_dir, usage, gtf_path)

    """
    -input
    database/Shibata/GSE241052_ari_org_annotated_fixed_normalized.h5ad
    -output
    ztest_folder
    -bk
    sample
    -lk
    celltype
    -m
    /Users/felixlang/Downloads/pipline_single/
    -u
    -gtf
    path to gtf file
    """

if __name__ == '__main__':
    main()