from pathlib import Path

from distributed.utils import palette
from matplotlib.colors import LinearSegmentedColormap
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import scanpy as sc
import tro_org.differential_expresison_analysis.plot_utils as plot_utils
import tro_org.utils.utils as utils

# followed the tutorial: https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_minimal_pydeseq2_pipeline.html#

def run_dea_deseq2(count_df, metadata, plot_dir, condition, contrast_cond):
    sc.set_figure_params(dpi_save=300, figsize=(8, 6), fontsize=14, vector_friendly=True)
    save_dir = Path(f"{plot_dir}")
    sc.settings.figdir = save_dir
    utils.ensure_dir(f"{plot_dir}")

    count_df = count_df[count_df.sum(axis=1) > 0]
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=count_df,
        metadata=metadata,
        design=f"~{condition}",
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    ds = DeseqStats(dds, contrast=["condition", contrast_cond[0], contrast_cond[1]], inference=inference)
    ds.summary()
    print(ds.LFC.columns)
    ds.lfc_shrink(coeff=f"{ds.LFC.columns[1]}") # VERY IMPORTANT takes care of exploding lfc values

    stat_result = ds.results_df
    # maybe filter low expressed genes:
    print(stat_result.shape)
    stat_result = stat_result[stat_result.baseMean >= 10]
    print(stat_result.shape)

    plot_cond = f"{contrast_cond[0]}_vs_{contrast_cond[1]}"

    DeseqDataSet.plot_dispersions(dds, save_path=f"{plot_dir}/deseq2_dispersions_{plot_cond}.png")
    DeseqStats.plot_MA(ds, save_path=f"{plot_dir}/deseq2_ma_plot_{plot_cond}.png")
    plot_utils.adjpvalue_hist(stat_result, plot_dir, plot_cond)
    plot_utils.vulcano_plot(stat_result, 2.0, 0.001, plot_dir, plot_cond)
    sc.tl.pca(dds)
    spezi_color = ["#ea6d3d", "#383a6b"]
    sc.pl.pca(dds,
              color="condition",
              size=200,
              palette=spezi_color,
              title= f"PCA {plot_cond.replace('_', ' ')}",
              save=f"_deseq2_pca_{plot_cond}.png")

    return ds.results_df