#!/bin/bash
# bash scripts/run_R_pipeline.sh     (from the repo root; after the upstream jobs have produced Salmon_Results_vM25/)
set -e
export REPO_DIR=$(pwd)
for s in 03_tximport 04_DESeq2 05_enrichment_ORA_GSEA 06_figures 07_supplementary_tables; do
  echo "================ $s"; Rscript scripts/R/$s.R 2>&1 | tee "results/log_$s.txt"
done
