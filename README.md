# Trophoblast Organoid Atlas

## QC and Preprocessing
```
preprocessing.py 
```
 + Controls the following QC and preprocessing functions
 + Main function to call 
 + Dataset input: currently with list Future: use nfcore/scdownstream samplesheet
 + Just for mapping CLI: -adata: single adata file path
```
dataset_analysis.py 
```
  * generates  plots for Gene-Level and Cell-Level Quality control metrics

```
map_genes.py 
```
+ Map used gene symbols to Ensembl IDs 
+ return preprocessed dataset

```
gene_venn.py 
```
+ plot Venn diagram before and after mapping 


## Data Analysis
+ tro_org/analysis/dataset/dataset_analysis.py -> plot dataset related (QC) metrics #
+ tro-org/analysis/cell_type_comparison.py
    * The datasets to be compared are stored in the config: "cosine_comp_config.yaml" with name, path and label column
    * Calculate cosine similarity
    * Calculate Jaccard similarity with marker genes (sc.get.ranks_gene_groups) and with HVG (TODO)
    * Plot Results in Heatmaps


+ tro-org/analysis/identify_genes.py
  * Identify batch specific genes and return heatmap


+ tro-org/analysis/cyte_type.py
  * Compare celltypes with Cytetype



## Workflow
[scDownstream](https://nf-co.re/scdownstream/dev/) performs quality control, integration, dimensionality reduction and clustering <br>
Integration methods: 
* scVI
* Harmony
* Combat 
* BBKNN

### Run scDownstream on DaisyBio:

| Flag                                     |               input                |
|------------------------------------------|:----------------------------------:|
| nextflow run nf-core/scdownstream -r dev |          run scDownstream          |
| -input                                   |      input file: sample sheet      |
| -outdir                                  |           output folder            |
| -profile                                 |        daisybio, apptainer         |
| --integration_methods                    |     choose integration method      |
| --integration_hvgs                       |            0: all genes            |
| -c                                       | custom config (mem_profile.config) |
| -with-report                             |   html Nextflow workflow report    |
| -with-trace                              |        trace in txt formal         |
| -with-timeline                           | html Processes execution timeline  |

Example:
```
nextflow run nf-core/scdownstream -r dev --input test --outdir out/test_combat --integration_methods combat --integration_hvgs 0 --prep_cellxgene -profile daisybio,apptainer -c /nfs/home/students/f.lang/Trophoblast-Organoid-Atlas/nf-core/custom_config/mem_profile.config -with-report out/report/combat_report.html -with-trace out/report/combat_trace.txt -with-timeline out/report/combat_timeline.html
```
### Run all integration methods at once:
```
pipeline/scdownstream_pipeline.sh
```
| Flag                     |                     input                      |
|--------------------------|:----------------------------------------------:|
| scdownstream_pipeline.sh |           run scdownstream_pipeline            |
| $1                       |            input file: sample sheet            |
| $2                       |                 output folder                  |
| $3 --monitor             |            bool flag for monitoring            |
| $4        -bg            | bool flag run integration in order pr parallel |


### Run self implemented integration tools on DaisyBio server <br>
Integration methods: 
* pyLiger
* Scanorama
* scGlue 

Are automatically done in the benchmarking step <br>

### Benchmark on DaisyBio server
Run scdownstream before <br>
Example:
```
sbatch slurm/benchmark-cpu.slurm
```

| Flag in slurm script             |                        input                        |
|----------------------------------|:---------------------------------------------------:|
| tro_org.benchmark.scIB_benchmark |                    run benchmark                    |
| -input                           |             single or merged .h5ad file             |
| -output                          |                    output folder                    |
| -bk / -batch_key                 |                  key batch column                   |
| -lk / --label_key                |                  key label column                   |
| -m                               | result / output dir scdownstream (in pipeline_info) |
| -u                               |               bool log resource usage               |
| -gtf                             |           gtf annotation file for scGlue            |


### Evaluate resource usage
```
benchmark/runtime_assessment.py
```
| Flags                           |                                input                                |
|---------------------------------|:-------------------------------------------------------------------:|
| benchmark/runtime_assessment.py |                           run assessment                            |
| -h5ad                           |                      h5ad file integrated file                      |
| -out                            |                            output folder                            |
| -log                            |                resource_usage dir from scdownstream                 |
| -tsv                            | resource_usage dir from self implemented tools in benchmark out dir |

### Trajectory analysis
Run trajectory pseudotime analysis with 4 different methods: PAGA, Phlower, Palantir and scFates
```
trajectory_analysis/trajectory.py
```
| Flags                             |                         input                         |
|-----------------------------------|:-----------------------------------------------------:|
| trajectory_analysis/trajectory.py |                run trajectory analysis                |
| -adata                            | yaml config file containing all datasets for analysis |
| -output                           |                     output folder                     |
| -lk                               |                   labelkey in adata                   |


### Gene Regulatory Networks
Construct Gene Reguatory Networks with pyScenic
```
GRN/pyscenic.py
```
| Flags           |                                               input                                                |
|-----------------|:--------------------------------------------------------------------------------------------------:|
| GRN/pyscenic.py |                                      run trajectory analysis                                       |
| -a              |                                 h5ad file of cluster / subcluster                                  |
| -i              |                                    image (singularity / docker)                                    |
| -d              |                                        data dir for docker                                         |
| -f              | filter adata according to [Van de Sande et al.](https://www.nature.com/articles/s41596-020-0336-2) |
| -m              |                                      activate multiprocessing                                      |
| -n              |                                         Number of workers                                          |


