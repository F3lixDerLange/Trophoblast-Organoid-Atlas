#!/bin/bash
# Metadata
#SBATCH --job-name scenic_Troph
#SBATCH --time=2-00:00:00
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
#SBATCH --cpus-per-task=32
#SBATCH --mem=156G

export TMPDIR=/nfs/home/students/f.lang/torg/TMP/$SLURM_JOB_ID
mkdir -p $TMPDIR

srun python3 -m tro_org.GRN.pyscenic -i singularity -a /nfs/home/students/f.lang/torg/GRN/data/subcluster_Trophoblast-CT_integrated.h5ad -d /nfs/home/students/f.lang/torg/GRN/data/ -n 64 -f -m