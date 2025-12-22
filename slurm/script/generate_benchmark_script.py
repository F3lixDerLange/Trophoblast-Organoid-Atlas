
# CPU specific
NUM_NODES = 1
NUM_TASKS = 1
CPUS_PER_TASK = 8
MEM = 64
OUT = f"/nfs/home/students/f.lang/torg/logs/%j/%u_%x.out"
ERROR = f"/nfs/home/students/f.lang/torg/logs/%j/%u_%x.err"

# GPU specific
PARTITION = "shared-gpu"
GPU_NUM = 1

# Working directory
WOR_DIR = "/nfs/home/students/f.lang/Trophoblast-Organoid-Atlas"

RUN_COMMAND = f"srun python3 -m tro_org.benchmark.scIB_benchmark"


def generate_benchmark_script(experiment_name, gpu, time, srun, script_dir):
    args = dict(
        job_name=experiment_name,
        time=time,
        chdir=WOR_DIR,
        num_nodes=NUM_NODES,
        num_tasks=NUM_TASKS,
        core_limit=CPUS_PER_TASK,
        memory_limit_gb=MEM,
        output_file=OUT,
        error_file=ERROR
    )

    gpu_args = dict(
        partition=PARTITION,
        num_gpus=GPU_NUM
    )

    if gpu:
        args.update(gpu_args) #check
        experiment_name = f"{experiment_name}_gpu"

    template = generate_sbatch_template(
        **args,
        command=srun
    )

    write_template(template, experiment_name, script_dir)

def write_template(template, experiment_name ,save_dir):

    with open(f"{save_dir}/benchmark_{experiment_name}.slurm", "w") as f:
        f.write(template)

def generate_sbatch_template(
        job_name: str,
        time: str,
        chdir: str,
        num_nodes: int,
        num_tasks: int,
        core_limit: int,
        memory_limit_gb: int,
        command: str,
        output_file: str = "stdout_%j.log",
        error_file: str = "stderr_%j.log",
        partition: str | None = None,
        num_gpus: int | None = None,
):
    gpu_lines = []
    if partition is not None and num_gpus is not None and num_gpus > 0:
        gpu_lines.append("# GPU specific")
        gpu_lines.append(f"#SBATCH --partition={partition}")
        gpu_lines.append(f"#SBATCH --gres=gpu:{num_gpus}")

    gpu_block = "\n".join(gpu_lines)

    sbatch = f"""#!/bin/bash
# Metadata
#SBATCH --job-name {job_name}
#SBATCH --time={time}
#SBATCH --no-requeue
# Working directory
#SBATCH --chdir {chdir}
# Environment
#SBATCH --export=ALL
#SBATCH --get-user-env
# Output and error
#SBATCH --output {output_file}
#SBATCH --error {error_file}
# CPU specific
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_tasks}
#SBATCH --cpus-per-task={core_limit}
#SBATCH --mem={memory_limit_gb}G
{gpu_block}

srun {command}
"""

    return sbatch

def generate_srun(input_file, output, bk, lk, m, u, gtf):
    srun = f"{RUN_COMMAND} -input {input_file} -output {output} -bk {bk} -lk {lk} -m {m} -gtf {gtf} {'-u' if u else ''}"
    return srun


def main():
    gtf = "/nfs/proj/tropho_org_atlas/database/gene_annotation/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
    datasets = [["/Shibata_single_samplesheet_scvi/finalized/merged.h5ad", "batch", "label", "/nfs/home/students/f.lang/nf-core_out/pipline_single"],
                ["/annotated_samplesheet_scvi/finalized/merged.h5ad", "batch", "label", "/Users/felixlang/Downloads/ann_integration_in_order"]
                ]
    out_dir = "/nfs/home/students/f.lang/torg/benchmark/"
    sbatch_time = "8:00:00"
    gpu = False
    monitor_usage = True
    save_script_dir = "slurm/test"
    for dataset in datasets:
        experiment_name = dataset[3].split("/")[-1]
        print(f"Generating {experiment_name}")
        srun = generate_srun(f"{dataset[3]}{dataset[0]}",
                             f"{out_dir}{experiment_name}",
                             dataset[1],
                             dataset[2],
                             dataset[3],
                             monitor_usage,
                             gtf)
        print(f"Srun: {srun}")
        generate_benchmark_script(experiment_name, gpu, sbatch_time, srun, save_script_dir)

if __name__ == '__main__':
    main()






