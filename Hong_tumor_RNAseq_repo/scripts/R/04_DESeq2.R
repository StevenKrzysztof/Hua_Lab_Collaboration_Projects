# =============================================================================
# 04_DESeq2.R — DESeq2 (~ group), QC, contrasts PD1 vs PBS within NT and sg7 + interaction, NT-vs-sg7 comparison
# Input : rds/01_txi.RData
# Output: DE/*.csv, QC/*.pdf, plots/Venn*, Quadrant*, rds/02_deseq.RData
# Run   : Rscript scripts/R/04_DESeq2.R      (from the repository root)
# =============================================================================
REPO_DIR <- Sys.getenv("REPO_DIR", getwd()); Sys.setenv(REPO_DIR = REPO_DIR)
source(file.path(REPO_DIR, "scripts", "R", "00_config.R"))
load(file.path(RDS, "01_txi.RData"))
dds <- DESeqDataSetFromTximport(txi, coldata, ~ group)
keep <- rep(FALSE, nrow(dds))
for (g in levels(coldata$group)) {
  cols <- which(coldata$group == g)
  keep <- keep | rowSums(counts(dds)[, cols, drop = FALSE] >= MIN_COUNT) >= MIN_SAMPLES
}
dds <- dds[keep, ]
cat("Genes after filter:", nrow(dds), "\n")
dds <- DESeq(dds)
print(resultsNames(dds))

# interaction = (sg7_PD1 - sg7_PBS) - (NT_PD1 - NT_PBS), in terms of the ~group coefficients
rn <- resultsNames(dds)
need <- c("group_NT_PD1_vs_NT_PBS", "group_sg7_PBS_vs_NT_PBS", "group_sg7_PD1_vs_NT_PBS")
stopifnot(all(need %in% rn))
int_vec <- setNames(rep(0, length(rn)), rn)
int_vec["group_sg7_PD1_vs_NT_PBS"] <- 1; int_vec["group_sg7_PBS_vs_NT_PBS"] <- -1; int_vec["group_NT_PD1_vs_NT_PBS"] <- -1

CONTRASTS <- list(
  NT_PD1_vs_PBS  = list(type = "name", value = c("group", "NT_PD1",  "NT_PBS")),
  sg7_PD1_vs_PBS = list(type = "name", value = c("group", "sg7_PD1", "sg7_PBS")),
  Interaction    = list(type = "num",  value = unname(int_vec)))
MAIN <- c("NT_PD1_vs_PBS", "sg7_PD1_vs_PBS")

vsd <- vst(dds, blind = FALSE)
vmat <- assay(vsd)

# ---- QC ----
pca  <- prcomp(t(vmat[head(order(rowVars(vmat), decreasing = TRUE), 500), ]))
pv   <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pdat <- data.frame(pca$x[, 1:3], coldata, sample = rownames(coldata))
for (ax in list(c("PC1", "PC2"), c("PC2", "PC3"))) {
  p <- ggplot(pdat, aes(.data[[ax[1]]], .data[[ax[2]]], color = genotype, shape = treatment, label = sample)) +
    geom_point(size = 5) + geom_text_repel(size = 3, color = "black") +
    scale_color_manual(values = geno_cols) + theme_bw(base_size = 13) +
    labs(title = paste("PCA (top 500 variable genes):", ax[1], "vs", ax[2]),
         x = paste0(ax[1], " ", pv[as.integer(sub("PC", "", ax[1]))], "%"),
         y = paste0(ax[2], " ", pv[as.integer(sub("PC", "", ax[2]))], "%"))
  save_plot(p, file.path(OUTDIR, "QC", paste0("PCA_", ax[1], "_", ax[2], ".pdf")), 8, 6)
}
sd_mat <- as.matrix(dist(t(vmat)))
ann_col <- data.frame(genotype = coldata$genotype, treatment = coldata$treatment, row.names = rownames(coldata))
ann_colors <- list(genotype = geno_cols, treatment = trt_cols, group = grp_cols)
pheatmap(sd_mat, annotation_col = ann_col, annotation_colors = ann_colors, main = "Sample distances (VST)",
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
pheatmap(sd_mat, annotation_col = ann_col, annotation_colors = ann_colors, main = "Sample distances (VST)",
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         filename = file.path(OUTDIR, "QC", "sample_distance_heatmap.pdf"), width = 8, height = 7)
cat("Note: ND_10R_2 and ND_11L_2 came from the 'org' delivery folder — check they sit with the other tumors here.\n")

# ---- differential expression ----
get_res <- function(ct) {
  r <- if (ct$type == "name") results(dds, contrast = ct$value, alpha = PADJ_CUT)
       else results(dds, contrast = ct$value, alpha = PADJ_CUT)
  as.data.frame(r) %>% rownames_to_column("gene_id") %>%
    left_join(gene_ann, by = "gene_id") %>%
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
           sig = !is.na(padj) & padj < PADJ_CUT & abs(log2FoldChange) > LFC_CUT,
           Direction = case_when(sig & log2FoldChange > 0 ~ "Up", sig & log2FoldChange < 0 ~ "Down", TRUE ~ "NS")) %>%
    arrange(padj)
}
ent_map <- suppressMessages(mapIds(org.Mm.eg.db, keys = unique(gene_ann$ensembl), keytype = "ENSEMBL",
                                   column = "ENTREZID", multiVals = "first"))
RES <- lapply(CONTRASTS, function(ct) { r <- get_res(ct); r$entrez <- unname(ent_map[r$ensembl]); r })
for (n in names(RES)) {
  write.csv(RES[[n]], file.path(OUTDIR, "DE", paste0("DE_", n, "_all.csv")), row.names = FALSE)
  write.csv(RES[[n]] %>% filter(sig), file.path(OUTDIR, "DE", paste0("DE_", n, "_sig.csv")), row.names = FALSE)
}
de_summary <- bind_rows(lapply(names(RES), function(n) { r <- RES[[n]]
  tibble(contrast = n, tested = sum(!is.na(r$padj)), up = sum(r$Direction == "Up"), down = sum(r$Direction == "Down")) }))
write.csv(de_summary, file.path(OUTDIR, "DE", "DE_summary.csv"), row.names = FALSE)
knitr::kable(de_summary)

# ---- NT vs sg7 PD1 response ----
cmp <- RES$NT_PD1_vs_PBS %>% select(gene_id, gene_name, LFC_NT = log2FoldChange, padj_NT = padj, sig_NT = sig) %>%
  inner_join(RES$sg7_PD1_vs_PBS %>% select(gene_id, LFC_sg7 = log2FoldChange, padj_sg7 = padj, sig_sg7 = sig), by = "gene_id") %>%
  left_join(RES$Interaction %>% select(gene_id, LFC_int = log2FoldChange, padj_int = padj), by = "gene_id") %>%
  mutate(category = case_when(sig_NT & sig_sg7 & sign(LFC_NT) == sign(LFC_sg7) ~ "Shared",
                              sig_NT & sig_sg7 ~ "Opposite",
                              sig_NT ~ "NT only", sig_sg7 ~ "sg7 only", TRUE ~ "Neither"))
write.csv(cmp, file.path(OUTDIR, "DE", "Compare_NT_vs_sg7_PD1_response.csv"), row.names = FALSE)
print(table(cmp$category))

for (dirn in c("Up", "Down")) {
  a <- RES$NT_PD1_vs_PBS$gene_id[RES$NT_PD1_vs_PBS$Direction == dirn]
  b <- RES$sg7_PD1_vs_PBS$gene_id[RES$sg7_PD1_vs_PBS$Direction == dirn]
  safe_venn(a, b, c("NT", "sg7"), c("#2171B5", "#D94801"),
            file.path(OUTDIR, "plots", paste0("Venn_PD1_", dirn, "_NT_vs_sg7.pdf")), main = paste("PD1", dirn, "genes"))
}

d <- cmp %>% filter(category != "Neither")
lab <- bind_rows(d %>% arrange(pmin(padj_NT, padj_sg7, na.rm = TRUE)) %>% head(20),
                 d %>% filter(gene_name %in% GENES_OF_INTEREST),
                 d %>% filter(!is.na(padj_int), padj_int < PADJ_CUT) %>% arrange(padj_int) %>% head(10)) %>%
  distinct(gene_id, .keep_all = TRUE)
lim <- max(abs(c(d$LFC_NT, d$LFC_sg7)), na.rm = TRUE)
p <- ggplot(d, aes(LFC_NT, LFC_sg7, color = category)) +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
  geom_abline(linetype = "dashed", color = "grey50") + geom_point(alpha = 0.7, size = 1.4) +
  geom_text_repel(data = lab, aes(label = gene_name), color = "black", size = 3, max.overlaps = 40) +
  scale_color_manual(values = c(Shared = "purple", `NT only` = "#2171B5", `sg7 only` = "#D94801", Opposite = "darkgreen")) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) + theme_bw(base_size = 13) +
  labs(title = "PD1 response: NT vs sg7 (genes significant in either)",
       subtitle = "on the diagonal = same response; off-diagonal = genotype-dependent",
       x = "log2FC PD1 vs PBS in NT", y = "log2FC PD1 vs PBS in sg7")
save_plot(p, file.path(OUTDIR, "plots", "Quadrant_PD1_response_NT_vs_sg7.pdf"), 8, 8)

save(dds, vsd, vmat, CONTRASTS, MAIN, RES, ent_map, de_summary, cmp, gene_ann, sym, txi, file = file.path(RDS, "02_deseq.RData"))
cat("Saved", file.path(RDS, "02_deseq.RData"), "\n")
