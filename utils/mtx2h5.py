import argparse
import scanpy as sc
import os
import inspect_h5ad

def convert_mtx2h5(mtx_folder):
    result_file_path_list = []
    for sub_entry in os.listdir(mtx_folder):
        if sub_entry.startswith('.'):
            continue
        print(f"Processing subfolder: {sub_entry}")
        sub_entry_path = os.path.join(mtx_folder, sub_entry)

        for filename in os.listdir(sub_entry_path):
            if filename.endswith(".mtx") or filename.endswith(".mtx.gz"):
                basename = os.path.basename(filename)
                if filename.endswith(".mtx.gz"):
                    prefix_mtx = basename.removesuffix("matrix.mtx.gz")
                elif filename.endswith(".mtx"):
                    prefix_mtx = basename.removesuffix("matrix.mtx")
                print(prefix_mtx)

                try:
                    adata = sc.read_10x_mtx(sub_entry_path,
                                            var_names='gene_symbols',
                                            cache=True,
                                            prefix=prefix_mtx,
                    )
                except Exception as e:
                    print(f"Error using sc.read_10x_mtx: {e}")

                h5ad_path = os.path.join(sub_entry_path, f"{prefix_mtx}merged.h5ad")
                adata.write(h5ad_path)
                print(f"Converted .h5")
                test = sc.read_h5ad(h5ad_path)
                print(f"Test shape: {test.shape}")
                result_file_path_list.append(h5ad_path)

    return result_file_path_list


def print_info(file2print):
    inspect_h5ad.print_h5ad_info(file2print)


def main():
    parser = argparse.ArgumentParser(description='Convert Mtx to h5')
    parser.add_argument('-mtx', type=str, help='Mtx file')
    parser.add_argument('-pinfo', action='store_true', help='print info about h5ad file')
    args = parser.parse_args()
    mtx = args.mtx
    pinfo = args.pinfo

    result_file_path = convert_mtx2h5(mtx)

    if pinfo:
        for file_path in result_file_path:
            print_info(file_path)

if __name__ == '__main__':
    main()
