#!/bin/bash
# Job metadata
#SBATCH --job-name benchmark_integration
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
# GPU specific

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE"/etc/profile.d/conda.sh
conda activate test-env

srun python3 -m tro_org.benchmark.scIB_benchmark -input /nfs/proj/tropho_org_atlas/database/Shibata/annotated/GSE241052_ari_org_annotated_fixed_normalized.h5ad -output /nfs/home/students/f.lang/torg

conda deactivate

