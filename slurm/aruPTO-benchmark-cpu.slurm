#!/bin/bash
# Job metadata
#SBATCH --job-name aru-benchmark_integration
#SBATCH --time=01:00:00
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

srun python3 -m tro_org.benchmark.scIB_benchmark -input /nfs/proj/tropho_org_atlas/database/Arutyunyan/Arutyunyan_PTO/Organoid_PTO_cellxgene_fixed.h5ad -output /nfs/home/students/f.lang/torg/benchmark/ArutyunyanPTO -bk sample -lk cell_annotation

conda deactivate

