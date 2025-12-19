# Trophoblast Organoid Atlas

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

Is automatically done in the benchmarking step <br>

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


