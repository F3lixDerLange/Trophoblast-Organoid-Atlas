#!/bin/bash
# Metadata
#SBATCH --job-name scenic_troph
#SBATCH --time=4-00:00:00
#SBATCH --no-requeue
# Working directory
#SBATCH --chdir /nfs/home/students/f.lang/Trophoblast-Organoid-Atlas
# Environment
#SBATCH --export=ALL
#SBATCH --get-user-env
# Output and error
#SBATCH --output /nfs/home/students/f.lang/torg/logs/%j/%u_%x.out
#SBATCH --error /nfs/home/students/f.lang/torg/logs/%j/%u_%x.err
# CPU specific
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export TMPDIR=$SCRATCH/tmp_$SLURM_JOB_ID
mkdir -p "$TMPDIR"
export DASK_TEMPORARY_DIRECTORY="$TMPDIR"

srun python3 -m tro_org.GRN.pyscenic -i singularity -a /nfs/home/students/f.lang/torg/GRN/data/subcluster_Trophoblast_integrated.h5ad -d /nfs/home/students/f.lang/torg/GRN/data/ -f