import argparse
import h5py
import scanpy as sc
import os

def convert_mtx2h5(mtx_folder):
    for filename in os.listdir(mtx_folder):
        if filename.endswith(".mtx") or filename.endswith(".mtx.gz"):
            basename = os.path.basename(filename)
            prefix_mtx = basename.removesuffix("matrix.mtx.gz")
            print(prefix_mtx, filename)
    try:
        adata = sc.read_10x_mtx(mtx_folder,
                                var_names='gene_symbols',
                                cache=True,
                                prefix=prefix_mtx
        )
    except Exception as e:
        print(f"Error using sc.read_10x_mtx: {e}")

    adata.write(os.path.join(mtx_folder, f"{prefix_mtx}merged.h5ad"))
    print(f"Converted .h5")
    test = sc.read_h5ad(os.path.join(mtx_folder, f"{prefix_mtx}merged.h5ad"))
    print(f"Test shape: {test.shape}")


def main():
    parser = argparse.ArgumentParser(description='Convert Mtx to h5')
    parser.add_argument('-mtx', type=str, help='Mtx file')
    #parser.add_argument('h5', type=str, help='h5 file')
    args = parser.parse_args()

    mtx = args.mtx
    #h5 = args.h5

    convert_mtx2h5(mtx)

if __name__ == '__main__':
    main()