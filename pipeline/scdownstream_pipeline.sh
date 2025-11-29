#!/bin/bash

sample_sheet=$1
output_dir=$2
filename="${sample_sheet##*/}"
methods=("scvi" "harmony" "bbknn" "combat")

for tool in "${methods[@]}"; do
  outpath="$output_dir/${filename%.*}_${tool}"
  mkdir -p "$outpath"
  echo "Starting nf-core/scdonstream with method $tool..."
  nextflow run nf-core/scdownstream \
    -r dev \
    --input "$sample_sheet" \
    --outdir "$outpath" \
    --integration_methods "$tool" \
    --integration_hvgs 0 \
    -profile daisybio,apptainer \
    -with-report "$outpath/report/${tool}_report.html" \
    -with-trace "$outpath/report/${tool}_trace.txt" \
    -with-timeline "$outpath/report/${tool}_timeline.html" \
    2> "$outpath/report/${tool}_time.txt" \
    &
done



