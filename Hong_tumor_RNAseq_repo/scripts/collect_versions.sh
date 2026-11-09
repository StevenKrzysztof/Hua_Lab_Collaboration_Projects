#!/bin/bash
# bash scripts/collect_versions.sh   -> writes config/software_versions.txt (run on the cluster, from the repo root)
OUT=config/software_versions.txt
source ~/.bashrc
source /ix1/hzhang/SteveWang/TEtranscripts_venv/bin/activate 2>/dev/null
conda activate /vast/bnacev/ziw183/conda_envs/nacevlab
v() { printf "%-16s %s\n" "$1" "$($2 2>&1 | grep -m1 -E '[0-9]+\.[0-9]+' | tr -s ' ')"; }
{
  echo "Software versions — collected $(date '+%Y-%m-%d') on $(hostname)"
  echo "=================================================================="
  echo "conda env: /vast/bnacev/ziw183/conda_envs/nacevlab ; TE tools: /ix1/hzhang/SteveWang/TEtranscripts_venv"; echo
  v cutadapt  "cutadapt --version"
  v STAR      "STAR --version"
  v samtools  "samtools --version"
  v salmon    "salmon --version"
  v TElocal   "TElocal --version"
  v TEcount   "TEcount --version"
  v R         "R --version"
  echo; echo "R / Bioconductor packages"; echo "-------------------------"
  Rscript -e 'p <- c("tximport","DESeq2","clusterProfiler","enrichplot","org.Mm.eg.db","msigdbr","fgsea","AnnotationDbi",
                    "data.table","dplyr","tidyr","ggplot2","ggrepel","pheatmap","VennDiagram","matrixStats","rtracklayer","ChIPseeker");
              for (x in p) cat(sprintf("%-16s %s\n", x, tryCatch(as.character(packageVersion(x)), error = function(e) "not installed")));
              cat(sprintf("%-16s %s\n", "Bioconductor", tryCatch(as.character(BiocManager::version()), error = function(e) "unknown")))'
  echo; echo "MSigDB release used by msigdbr:"
  Rscript -e 'suppressMessages(library(msigdbr)); x <- tryCatch(msigdbr_version(), error = function(e) NA); cat(if (is.na(x[1])) "see msigdbr package version above" else x, "\n")' 2>/dev/null
} > "$OUT"
cat "$OUT"
