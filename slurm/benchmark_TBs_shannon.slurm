#!/bin/bash
# Metadata
#SBATCH --job-name TBs_shannon
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
#SBATCH --mem=64G

srun python3 -m tro_org.benchmark.scIB_benchmark -input /nfs/home/students/f.lang/nf-core_out/TBs_shannon/TBs_Shannon_samplesheet_scvi/finalized/merged.h5ad -output /nfs/home/students/f.lang/torg/benchmark/TBs_shannon -bk batch -lk label -m /nfs/home/students/f.lang/nf-core_out/TBs_shannon -gtf /nfs/proj/tropho_org_atlas/database/gene_annotation/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz -u
