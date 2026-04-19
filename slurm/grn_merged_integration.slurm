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
#SBATCH --cpus-per-task=32
#SBATCH --mem=196G

srun python3 -m tro_org.GRN.pyscenic -i singularity -a /nfs/home/students/f.lang/torg/GRN/data/merged_integration_final_integrated.h5ad -d /nfs/home/students/f.lang/torg/GRN/data/ -n 64 -f -m