# =============================================================================
# 03_tximport.R — Salmon quant.sf -> gene-level counts with tximport (GENCODE vM25 tx2gene)
# Input : Salmon_Results_vM25/<sample>/quant.sf, tx2gene_vM25.tsv
# Output: rds/01_txi.RData
# Run   : Rscript scripts/R/03_tximport.R      (from the repository root)
# =============================================================================
REPO_DIR <- Sys.getenv("REPO_DIR", getwd()); Sys.setenv(REPO_DIR = REPO_DIR)
source(file.path(REPO_DIR, "scripts", "R", "00_config.R"))
tx2 <- fread(TX2GENE, data.table = FALSE)
gene_ann <- tx2 %>% distinct(gene_id, gene_name, gene_type) %>% mutate(ensembl = strip_ver(gene_id))
sym <- setNames(gene_ann$gene_name, gene_ann$gene_id)

files <- file.path(SALMON_DIR, samples$sample, "quant.sf"); names(files) <- samples$sample
stopifnot(all(file.exists(files)))
q1 <- fread(files[1], select = "Name")$Name
miss <- mean(!(q1 %in% tx2$transcript_id))
cat(sprintf("Transcripts in quant.sf missing from tx2gene: %.3f%%\n", 100 * miss))
if (miss > 0.01) stop("More than 1% of transcripts are missing from tx2gene — index and annotation do not match.")

txi <- tximport(files, type = "salmon", tx2gene = tx2[, c("transcript_id", "gene_id")])
cat("Genes:", nrow(txi$counts), "\n")
for (g in c("Gzmb", "Prf1", "Cbx5")) {
  id <- gene_ann$gene_id[gene_ann$gene_name == g][1]
  cat(sprintf("  %-5s %s  counts: %s\n", g, id, paste(round(txi$counts[id, ]), collapse = " ")))
}

save(txi, tx2, gene_ann, sym, files, file = file.path(RDS, "01_txi.RData"))
cat("Saved", file.path(RDS, "01_txi.RData"), "\n")
