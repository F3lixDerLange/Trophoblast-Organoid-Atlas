#!/bin/bash
# Job metadata
#SBATCH --job-name shib-benchmark_integration
#SBATCH --time=04:00:00
#SBATCH --no-requeue
# Working directory
#SBATCH --chdir=/nfs/home/students/f.lang/Trophoblast-Organoid-Atlas
# Environment
#SBATCH --export=ALL
#SBATCH --get-user-env
# Output and error
#SBATCH --output /nfs/home/students/f.lang/torg/logs/%j/%u_%x.out
#SBATCH --error /nfs/home/students/f.lang/torg/logs/%j/%u_%x.err
# CPU specific
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

conda activate test-env

srun python3 -m tro_org.benchmark.scIB_benchmark -input /nfs/home/students/f.lang/nf-core_out/pipline_single/Shibata_single_samplesheet_scvi/finalized/merged.h5ad -output /nfs/home/students/f.lang/torg/benchmark/shibata -bk batch -lk label -m /nfs/home/students/f.lang/nf-core_out/pipline_single

conda deactivate

