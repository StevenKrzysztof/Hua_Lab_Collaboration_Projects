# =============================================================================
# 00_config.R — libraries, paths, samples, thresholds, colours, helper functions.
# Sourced by every numbered R script. Paths can be overridden with environment variables:
#   TUMOR_RNASEQ_BASE (folder containing Salmon_Results_vM25/), TUMOR_RNASEQ_OUT, TX2GENE_VM25, REPO_DIR
# =============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages({
  library(data.table); library(matrixStats)
  library(tximport); library(DESeq2)
  library(clusterProfiler); library(enrichplot); library(org.Mm.eg.db); library(msigdbr)
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
  library(VennDiagram); library(grid); library(scales)
})
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
count <- dplyr::count; select <- dplyr::select; filter <- dplyr::filter
first <- dplyr::first; last <- dplyr::last; between <- dplyr::between
set.seed(1)

# ---- INPUT / OUTPUT LOCATIONS (edit these two lines for a new system) ----
BASE       <- Sys.getenv("TUMOR_RNASEQ_BASE", "/ix1/hzhang/SteveWang/Hong_data/tumor_RNA_Seq")
SALMON_DIR <- file.path(BASE, "Salmon_Results_vM25")                       # {sample}/quant.sf
TX2GENE    <- Sys.getenv("TX2GENE_VM25", "/ix1/bnacev/ziw183/ref/grcm38/vM25/tx2gene_vM25.tsv")   # from scripts/upstream/00_build_refs_vM25.sh
OUTDIR     <- Sys.getenv("TUMOR_RNASEQ_OUT", file.path(BASE, "DESeq2_vM25_PD1_vs_PBS"))
RDS        <- file.path(OUTDIR, "rds"); dir.create(RDS, recursive = TRUE, showWarnings = FALSE)
for (d in c("", "QC", "DE", "plots", "enrichment/ORA", "enrichment/GSEA", "supplementary"))
  dir.create(file.path(OUTDIR, d), recursive = TRUE, showWarnings = FALSE)

REPO <- Sys.getenv("REPO_DIR", getwd())   # repository root (scripts set this)
samples <- read.csv(file.path(REPO, "metadata", "sample_metadata.csv"), stringsAsFactors = FALSE)
samples$group <- paste(samples$genotype, samples$treatment, sep = "_")
coldata <- data.frame(row.names = samples$sample,
                      genotype  = factor(samples$genotype,  levels = c("NT", "sg7")),
                      treatment = factor(samples$treatment, levels = c("PBS", "PD1")),
                      group     = factor(samples$group, levels = c("NT_PBS", "NT_PD1", "sg7_PBS", "sg7_PD1")))

PADJ_CUT <- 0.05; LFC_CUT <- 1
MIN_COUNT <- 10; MIN_SAMPLES <- 2          # >= 10 counts in >= 2 of 3 samples of at least one group
gs <- read.csv(file.path(REPO, "metadata", "gene_sets.csv"), stringsAsFactors = FALSE)
GENES_OF_INTEREST <- gs$gene[gs$set == "genes_of_interest"]
MODULES <- split(gs$gene[gs$set != "genes_of_interest"], gs$set[gs$set != "genes_of_interest"])
MODULES <- MODULES[unique(gs$set[gs$set != "genes_of_interest"])]

grp_cols  <- c(NT_PBS = "#9ECAE1", NT_PD1 = "#2171B5", sg7_PBS = "#FDAE6B", sg7_PD1 = "#D94801")
geno_cols <- c(NT = "#2171B5", sg7 = "#D94801"); trt_cols <- c(PBS = "grey80", PD1 = "grey30")
hm_cols   <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
ann_colors <- list(genotype = geno_cols, treatment = trt_cols, group = grp_cols)

# ---- helper functions ----
strip_ver <- function(x) sub("\\..*$", "", x)
save_plot <- function(p, path, w = 8, h = 6) { ggsave(path, p, width = w, height = h); print(p) }
zrow      <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
nrow0     <- function(x) if (is.null(x)) 0L else nrow(as.data.frame(x))

safe_venn <- function(a, b, labels, fill, path, main = NULL) {
  if (length(a) + length(b) == 0) { message("Venn skipped (empty): ", basename(path)); return(invisible(NULL)) }
  v <- draw.pairwise.venn(length(a), length(b), length(intersect(a, b)), category = labels,
                          fill = fill, alpha = 0.5, cex = 1.6, cat.cex = 1.2, ind = FALSE)
  draw <- function() { grid.newpage(); if (!is.null(main)) grid.text(main, y = 0.95, gp = gpar(fontsize = 13)); grid.draw(v) }
  draw(); pdf(path, 7, 6); draw(); dev.off()
}

# msigdbr changed arguments/columns in v10 — handle both
get_hallmark <- function() {
  m <- tryCatch(msigdbr(species = "Mus musculus", collection = "H"),
                error = function(e) msigdbr(species = "Mus musculus", category = "H"))
  ens <- intersect(c("ensembl_gene", "db_ensembl_gene"), names(m))[1]
  ent <- intersect(c("entrez_gene", "ncbi_gene", "db_ncbi_gene"), names(m))[1]
  list(ens = unique(data.frame(term = m$gs_name, gene = m[[ens]])),
       ent = unique(data.frame(term = m$gs_name, gene = as.character(m[[ent]]))))
}
