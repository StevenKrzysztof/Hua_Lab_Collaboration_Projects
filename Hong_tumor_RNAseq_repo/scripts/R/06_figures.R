# =============================================================================
# 06_figures.R — volcano plots (shared axes), heatmaps (VST z-scores), immune module heatmaps and scores
# Input : rds/02_deseq.RData
# Output: plots/*.pdf, plots/immune_module_scores.csv
# Run   : Rscript scripts/R/06_figures.R      (from the repository root)
# =============================================================================
REPO_DIR <- Sys.getenv("REPO_DIR", getwd()); Sys.setenv(REPO_DIR = REPO_DIR)
source(file.path(REPO_DIR, "scripts", "R", "00_config.R"))
load(file.path(RDS, "02_deseq.RData"))

# ---- volcano plots ----
xl <- quantile(abs(unlist(lapply(RES[MAIN], function(r) r$log2FoldChange))), 0.999, na.rm = TRUE)
ally <- unlist(lapply(RES[MAIN], function(r) -log10(r$padj))); ally <- ally[is.finite(ally)]
yl <- quantile(ally, 0.999, na.rm = TRUE)
xl <- max(ceiling(xl), 2); yl <- max(ceiling(yl), 3)

plot_volcano <- function(r, title, path, xlim = NULL, ylim = NULL, n_label = 15) {
  d <- r %>% filter(!is.na(padj)) %>% mutate(y = pmin(-log10(pmax(padj, 1e-300)), 300))
  if (!is.null(xlim)) d <- d %>% mutate(out = abs(log2FoldChange) > xlim | y > ylim,
                                        log2FoldChange = pmax(pmin(log2FoldChange, xlim), -xlim), y = pmin(y, ylim))
  else d$out <- FALSE
  lab <- bind_rows(d %>% filter(sig) %>% arrange(padj) %>% head(n_label),
                   d %>% filter(gene_name %in% GENES_OF_INTEREST)) %>% distinct(gene_id, .keep_all = TRUE)
  n_up <- sum(d$Direction == "Up"); n_dn <- sum(d$Direction == "Down")
  p <- ggplot(d, aes(log2FoldChange, y, color = Direction)) +
    geom_point(aes(shape = out), alpha = 0.6, size = 1.2) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_color_manual(values = c(NS = "grey75", Up = "#D7301F", Down = "#2171B5")) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = lab, aes(label = gene_name), color = "black", size = 3,
                    max.overlaps = 40, min.segment.length = 0, box.padding = 0.4) +
    theme_bw(base_size = 13) +
    labs(title = title, subtitle = paste0("Up: ", n_up, "   Down: ", n_dn,
                                          "   (padj < ", PADJ_CUT, ", |log2FC| > ", LFC_CUT, ")"),
         x = "log2 fold change", y = "-log10 adjusted p")
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = c(-xlim, xlim), ylim = c(0, ylim))
  save_plot(p, path, 8, 7)
}
for (n in MAIN) plot_volcano(RES[[n]], gsub("_", " ", n), file.path(OUTDIR, "plots", paste0("Volcano_", n, ".pdf")), xl, yl)
plot_volcano(RES$Interaction, "Interaction: sg7 PD1 response minus NT PD1 response",
             file.path(OUTDIR, "plots", "Volcano_Interaction.pdf"))

# ---- heatmaps ----
ord <- order(coldata$group)
gaps <- cumsum(table(coldata$group[ord]))[-nlevels(coldata$group)]
ann_c <- data.frame(group = coldata$group, genotype = coldata$genotype, treatment = coldata$treatment,
                    row.names = rownames(coldata))[ord, ]
sym <- setNames(gene_ann$gene_name, gene_ann$gene_id)

plot_hm <- function(ids, title, path, max_named = 80, row_ann = NULL) {
  ids <- intersect(ids, rownames(vmat))
  if (length(ids) < 2) { message("Heatmap skipped (<2 genes): ", title); return(invisible(NULL)) }
  m <- zrow(vmat[ids, ord, drop = FALSE]); rownames(m) <- make.unique(sym[rownames(m)])
  if (!is.null(row_ann)) rownames(row_ann) <- make.unique(sym[ids])
  args <- list(mat = m, color = hm_cols, breaks = seq(-2, 2, length.out = 256), cluster_cols = FALSE,
               cluster_rows = is.null(row_ann), gaps_col = gaps, annotation_col = ann_c, annotation_row = row_ann,
               annotation_colors = ann_colors, show_rownames = length(ids) <= max_named, fontsize_row = 7,
               main = paste0(title, " (n = ", length(ids), ")"))
  do.call(pheatmap, args)
  do.call(pheatmap, c(args, list(filename = path, width = 8, height = max(5, min(16, 0.13 * length(ids) + 3)))))
}

best <- cmp %>% filter(category != "Neither") %>% mutate(b = pmin(padj_NT, padj_sg7, na.rm = TRUE)) %>% arrange(b)
plot_hm(head(best$gene_id, 50), "Top 50 PD1-responsive genes (NT or sg7)", file.path(OUTDIR, "plots", "Heatmap_top50_PD1_genes.pdf"))
plot_hm(best$gene_id, "All PD1-responsive genes (NT or sg7)", file.path(OUTDIR, "plots", "Heatmap_all_PD1_genes.pdf"))
for (cat_i in c("Shared", "NT only", "sg7 only"))
  plot_hm(cmp$gene_id[cmp$category == cat_i], paste("PD1 genes —", cat_i),
          file.path(OUTDIR, "plots", paste0("Heatmap_PD1_", gsub(" ", "_", cat_i), ".pdf")))
int_sig <- RES$Interaction %>% filter(sig) %>% pull(gene_id)
plot_hm(int_sig, "Interaction genes (sg7 changes the PD1 response)", file.path(OUTDIR, "plots", "Heatmap_interaction_genes.pdf"))

# immune modules, rows grouped by module
mod_df <- bind_rows(lapply(names(MODULES), function(m) data.frame(module = m, gene_name = MODULES[[m]])))
mod_df <- mod_df %>% left_join(gene_ann %>% select(gene_id, gene_name), by = "gene_name") %>%
  filter(gene_id %in% rownames(vmat)) %>% distinct(gene_id, .keep_all = TRUE)
cat("Module genes found:", nrow(mod_df), "of", sum(lengths(MODULES)),
    "| missing:", paste(setdiff(unlist(MODULES), mod_df$gene_name), collapse = ", "), "\n")
if (nrow(mod_df) >= 2) {
  ra <- data.frame(module = mod_df$module)
  ann_colors$module <- setNames(brewer.pal(max(3, length(MODULES)), "Set2")[seq_along(MODULES)], names(MODULES))
  plot_hm(mod_df$gene_id, "Immune modules", file.path(OUTDIR, "plots", "Heatmap_immune_modules.pdf"),
          max_named = 200, row_ann = ra)
  # group-mean version: one column per condition
  gm <- sapply(levels(coldata$group), function(g) rowMeans(vmat[mod_df$gene_id, coldata$group == g, drop = FALSE]))
  gz <- zrow(gm); rownames(gz) <- make.unique(mod_df$gene_name)
  rownames(ra) <- rownames(gz)
  pheatmap(gz, color = hm_cols, breaks = seq(-1.5, 1.5, length.out = 256), cluster_cols = FALSE, cluster_rows = FALSE,
           annotation_row = ra, annotation_colors = ann_colors, main = "Immune modules — group means (z)",
           gaps_row = cumsum(table(factor(mod_df$module, levels = names(MODULES))))[-length(MODULES)],
           filename = file.path(OUTDIR, "plots", "Heatmap_immune_modules_group_means.pdf"), width = 6, height = 10)
}

# ---- immune module scores (custom signatures) ----
# Score = mean of per-gene z-scores (VST, z-scored across all 12 tumors) over the module genes present.
mod_scores <- sapply(names(MODULES), function(m) {
  ids <- gene_ann$gene_id[gene_ann$gene_name %in% MODULES[[m]]]; ids <- intersect(ids, rownames(vmat))
  if (length(ids) < 2) return(rep(NA_real_, ncol(vmat)))
  colMeans(zrow(vmat[ids, , drop = FALSE]))
})
rownames(mod_scores) <- colnames(vmat)
ms <- data.frame(sample = rownames(mod_scores), coldata[rownames(mod_scores), ], mod_scores, check.names = FALSE)
write.csv(ms, file.path(OUTDIR, "plots", "immune_module_scores.csv"), row.names = FALSE)
msl <- ms %>% pivot_longer(all_of(names(MODULES)), names_to = "module", values_to = "score")
p <- ggplot(msl, aes(group, score, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) + facet_wrap(~ module, scales = "free_y") +
  scale_fill_manual(values = grp_cols) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune module scores (mean z-scored VST)", x = NULL, y = "module score")
save_plot(p, file.path(OUTDIR, "plots", "immune_module_scores.pdf"), 11, 7)
