#!/bin/bash

sample_sheet=$1
output_dir=$2
filename="${sample_sheet##*/}"
methods=("scvi" "harmony" "bbknn" "combat")

for tool in "${methods[@]}"; do
  outpath="$output_dir/${filename%.*}_${tool}" # important: last part of the name is the tool. needed for benchmark
  report_path="$output_dir/report"
  mkdir -p "$outpath"
  mkdir -p "$report_path"
  echo "Starting nf-core/scdownstream with method $tool..."
  nextflow run nf-core/scdownstream \
    -r dev \
    --input "$sample_sheet" \
    --outdir "$outpath" \
    --integration_methods "$tool" \
    -profile daisybio,apptainer \
    -c /nfs/home/students/f.lang/Trophoblast-Organoid-Atlas/nf-core/custom_config/mem_profile.config \
    -with-report "${report_path}/${tool}_report.html" \
    -with-trace "${report_path}/${tool}_trace.txt" \
    -with-timeline "${report_path}/${tool}_timeline.html" \
    2> "${report_path}/${tool}_time.txt" \
    --prep_cellxgene
    # &    # use & if all tools should run in parallel (has maybe an impact on resource usage)
done

