#!/bin/bash
# Metadata
#SBATCH --job-name trajectory
#SBATCH --time=8:00:00
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
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

srun python3 -m tro_org.trajectory_analysis.trajectory -adata /nfs/home/students/f.lang/Trophoblast-Organoid-Atlas/tro_org/trajectory_analysis/trajectory_datasets.yaml -output /nfs/home/students/f.lang/torg/trajectory_analysis -lk label