# =============================================================================
# 07_supplementary_tables.R — raw counts, TPM, normalized counts, VST, sample sheet, parameters + sessionInfo
# Input : rds/02_deseq.RData
# Output: supplementary/*
# Run   : Rscript scripts/R/07_supplementary_tables.R      (from the repository root)
# =============================================================================
REPO_DIR <- Sys.getenv("REPO_DIR", getwd()); Sys.setenv(REPO_DIR = REPO_DIR)
source(file.path(REPO_DIR, "scripts", "R", "00_config.R"))
load(file.path(RDS, "02_deseq.RData"))
ann_out <- function(m) data.frame(gene_id = rownames(m), gene_name = sym[rownames(m)],
                                  gene_type = gene_ann$gene_type[match(rownames(m), gene_ann$gene_id)], m, check.names = FALSE)
sup <- file.path(OUTDIR, "supplementary")
write.csv(samples, file.path(sup, "Table_S_sample_sheet.csv"), row.names = FALSE)
write.csv(ann_out(round(txi$counts)), file.path(sup, "Table_S_raw_counts_all_genes.csv"), row.names = FALSE)
write.csv(ann_out(round(txi$abundance, 3)), file.path(sup, "Table_S_TPM_all_genes.csv"), row.names = FALSE)
write.csv(ann_out(round(counts(dds, normalized = TRUE), 2)), file.path(sup, "Table_S_normalized_counts_filtered.csv"), row.names = FALSE)
write.csv(ann_out(round(vmat, 4)), file.path(sup, "Table_S_VST_filtered.csv"), row.names = FALSE)
writeLines(c(paste("Salmon dir:", SALMON_DIR), paste("tx2gene:", TX2GENE),
             paste("Filter: >=", MIN_COUNT, "counts in >=", MIN_SAMPLES, "samples of >= 1 group"),
             paste("Design: ~ group;  thresholds padj <", PADJ_CUT, "and |log2FC| >", LFC_CUT),
             capture.output(sessionInfo())), file.path(sup, "analysis_parameters_and_sessionInfo.txt"))
cat("Done. Outputs in", OUTDIR, "\n")