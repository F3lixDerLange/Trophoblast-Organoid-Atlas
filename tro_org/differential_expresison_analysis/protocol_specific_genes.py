import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

def run_dea_sample_specific(count_df, metadata):

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=count_df,
        metadata=metadata,
        design="~dataset",
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()


    ds = DeseqStats(dds, contrast=["dataset", "Arutyunyan_TSC", "Arutyunyan_PTO", "Shannon", "Shibata"], inference=inference)
    ds.summary()

    print(dds)

