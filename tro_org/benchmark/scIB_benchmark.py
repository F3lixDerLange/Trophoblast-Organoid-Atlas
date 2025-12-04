import argparse
from scib_metrics.benchmark import *
import tro_org.integration.scGLUE_integration
import tro_org.integration.py_liger_integration
import tro_org.integration.scanorama_integration
import scanpy as sc

def run_integrations(base_adata, output_dir, batchkey, labelkey):
    osbm_keys = ["Unintegrated"]
    print("Running integrations")

    print("PyLiger integration")
    liger_adata = tro_org.integration.py_liger_integration.pyl_integration(base_adata.copy(), output_dir, batchkey, labelkey)
    osbm_keys.append("LIGER")

    print("scGlue integration")
    scGlue_adata = tro_org.integration.scGLUE_integration.scg_integration(base_adata.copy(), output_dir, batchkey, labelkey)
    osbm_keys.append("X_scglue")

    print("scanorama_integration")
    scanorama_adata = tro_org.integration.scanorama_integration.scn_integration(base_adata.copy(), output_dir, batchkey)
    osbm_keys.append("Scanorama")

    integrated_adata = base_adata.copy()

    integrated_adata.obsm["X_scglue"] = scGlue_adata.obsm["X_scglue"]
    integrated_adata.obsm["Scanorama"] = scanorama_adata.obsm["Scanorama"]
    integrated_adata.obsm["LIGER"] = liger_adata.obsm["LIGER"]



    return osbm_keys, integrated_adata


def benchmark(h5ad_file, output_dir, batch_key, label_key):

    osbm_keys, integrated_adata = run_integrations(h5ad_file, output_dir, batch_key, label_key)
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
    benchmark(sc.read_h5ad(input_file), out_dir, batch_key, label_key)

if __name__ == '__main__':
    main()