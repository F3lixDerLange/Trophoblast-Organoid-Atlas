#!/bin/bash

sample_sheet=$1
output_dir=$2
monitor=false
bg=false # start in parallel
filename="${sample_sheet##*/}"
methods=("scvi" "harmony" "bbknn" "combat")

if [[ "${3:-}" == "--monitor" ]]; then
  monitor=true
fi

if [[ "${4:-}" == "--bg" ]]; then
  bg=true
fi

for tool in "${methods[@]}"; do
  outpath="$output_dir/${filename%.*}_${tool}" # important: last part of the name is the tool. needed for benchmark
  report_path="$output_dir/report"
  mkdir -p "$outpath"
  mkdir -p "$report_path"
  echo "Starting nf-core/scdonstream with method $tool..."
  cmd=(
    nextflow run nf-core/scdownstream
    -r dev
    --input "$sample_sheet"
    --outdir "$outpath"
    --integration_methods "$tool"
    --integration_hvgs 0
    --prep_cellxgene
    -profile daisybio,apptainer
    )

  if $monitor; then
    cmd+=(
    -c /nfs/home/students/f.lang/Trophoblast-Organoid-Atlas/nf-core/custom_config/mem_profile.config
    -with-report "${report_path}/${tool}_report.html"
    -with-trace "${report_path}/${tool}_trace.txt"
    -with-timeline "${report_path}/${tool}_timeline.html"
    )
  fi

  if $bg; then # use & if all tools should run in parallel (has maybe an impact on resource usage)
    echo "${cmd[@]}" 2> "${report_path}/${tool}_time.txt" &
  else
    echo "${cmd[@]}" 2> "${report_path}/${tool}_time.txt"
  fi

done



