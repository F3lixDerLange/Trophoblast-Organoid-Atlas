import os

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def filter_invivo_cells(adata, obs_key="Model"):
    if obs_key not in adata.obs.columns:
        raise KeyError(f"Column '{obs_key}' not found in adata.obs")

    print(f"Adata shape before cell filtering: {adata.shape}")

    model_values = adata.obs[obs_key].astype(str)
    keep_mask = ~model_values.isin(['1_Shannon', '2_Shannon', '3_Shannon'])
    adata = adata[keep_mask]

    print(f"Adata shape after cell filtering: {adata.shape}")
    return adata