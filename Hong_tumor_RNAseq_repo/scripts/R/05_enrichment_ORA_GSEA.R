# =============================================================================
# 05_enrichment_ORA_GSEA.R — ORA (GO BP, Hallmark, KEGG) and GSEA (Hallmark, GO BP, KEGG) for every contrast
# Input : rds/02_deseq.RData
# Output: enrichment/ORA/*, enrichment/GSEA/*, rds/03_enrichment.RData
# Run   : Rscript scripts/R/05_enrichment_ORA_GSEA.R      (from the repository root)
# =============================================================================
REPO_DIR <- Sys.getenv("REPO_DIR", getwd()); Sys.setenv(REPO_DIR = REPO_DIR)
source(file.path(REPO_DIR, "scripts", "R", "00_config.R"))
load(file.path(RDS, "02_deseq.RData"))

# ---- ORA ----
HALL <- get_hallmark()
universe_ens <- unique(strip_ver(rownames(dds)))
universe_ent <- unique(na.omit(unname(ent_map[universe_ens])))

run_ora <- function(ens, label) {
  ens <- unique(ens); out <- list()
  if (length(ens) < 5) { message("ORA skipped (<5 genes): ", label); return(out) }
  ent <- unique(na.omit(unname(ent_map[ens])))
  out$GO_BP <- tryCatch(enrichGO(ens, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                                 universe = universe_ens, pAdjustMethod = "BH", readable = TRUE), error = function(e) NULL)
  out$Hallmark <- tryCatch(enricher(ens, TERM2GENE = HALL$ens, universe = universe_ens, pAdjustMethod = "BH"),
                           error = function(e) NULL)
  out$KEGG <- tryCatch(setReadable(enrichKEGG(ent, organism = "mmu", universe = universe_ent, pAdjustMethod = "BH"),
                                   OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
                       error = function(e) { message("KEGG ORA failed (needs internet?): ", conditionMessage(e)); NULL })
  for (db in names(out)) {
    x <- out[[db]]; if (nrow0(x) == 0) next
    df <- as.data.frame(x)
    if (db == "Hallmark") df$geneSymbol <- vapply(strsplit(df$geneID, "/"), function(g)
      paste(sym[gene_ann$gene_id[match(g, gene_ann$ensembl)]], collapse = "/"), character(1))
    write.csv(df, file.path(OUTDIR, "enrichment/ORA", paste0("ORA_", db, "_", label, ".csv")), row.names = FALSE)
    if (sum(df$p.adjust < 0.05) > 0) {
      p <- dotplot(x, showCategory = 15, font.size = 8) + ggtitle(paste("ORA", db, "—", gsub("_", " ", label)))
      save_plot(p, file.path(OUTDIR, "enrichment/ORA", paste0("ORA_", db, "_", label, ".pdf")), 9, 7)
    }
  }
  out
}
gid <- function(df, filt) strip_ver(df$gene_id[filt])
ORA_SETS <- list(
  NT_PD1_up    = gid(RES$NT_PD1_vs_PBS,  RES$NT_PD1_vs_PBS$Direction  == "Up"),
  NT_PD1_down  = gid(RES$NT_PD1_vs_PBS,  RES$NT_PD1_vs_PBS$Direction  == "Down"),
  sg7_PD1_up   = gid(RES$sg7_PD1_vs_PBS, RES$sg7_PD1_vs_PBS$Direction == "Up"),
  sg7_PD1_down = gid(RES$sg7_PD1_vs_PBS, RES$sg7_PD1_vs_PBS$Direction == "Down"),
  Shared_up    = gid(cmp, cmp$category == "Shared" & cmp$LFC_NT > 0),
  Shared_down  = gid(cmp, cmp$category == "Shared" & cmp$LFC_NT < 0),
  NTonly_up    = gid(cmp, cmp$category == "NT only"  & cmp$LFC_NT  > 0),
  NTonly_down  = gid(cmp, cmp$category == "NT only"  & cmp$LFC_NT  < 0),
  sg7only_up   = gid(cmp, cmp$category == "sg7 only" & cmp$LFC_sg7 > 0),
  sg7only_down = gid(cmp, cmp$category == "sg7 only" & cmp$LFC_sg7 < 0))
print(sapply(ORA_SETS, length))
ORA <- lapply(names(ORA_SETS), function(n) run_ora(ORA_SETS[[n]], n)); names(ORA) <- names(ORA_SETS)
ora_summary <- bind_rows(lapply(names(ORA), function(n) bind_rows(lapply(names(ORA[[n]]), function(db) {
  df <- if (nrow0(ORA[[n]][[db]])) as.data.frame(ORA[[n]][[db]]) else data.frame(p.adjust = numeric())
  tibble(set = n, database = db, n_genes = length(ORA_SETS[[n]]), sig_terms = sum(df$p.adjust < 0.05))
}))))
write.csv(ora_summary, file.path(OUTDIR, "enrichment/ORA", "ORA_summary.csv"), row.names = FALSE)
print(ora_summary, n = Inf)

# ---- GSEA ----
ranked <- function(r, id = c("ensembl", "entrez")) {
  id <- match.arg(id)
  d <- r %>% filter(!is.na(stat)) %>% mutate(key = if (id == "ensembl") ensembl else entrez) %>%
    filter(!is.na(key)) %>% arrange(desc(abs(stat))) %>% distinct(key, .keep_all = TRUE)
  sort(setNames(d$stat, d$key), decreasing = TRUE)
}
# call with eps/seed (newer clusterProfiler); retry without them on older versions
run_gsea <- function(fun, args, label) {
  res <- tryCatch(do.call(fun, args), error = function(e) e)
  if (inherits(res, "error") && grepl("unused argument", conditionMessage(res)))
    res <- tryCatch(do.call(fun, args[setdiff(names(args), c("eps", "seed"))]), error = function(e) e)
  if (inherits(res, "error")) { message("GSEA failed (", label, "): ", conditionMessage(res)); return(NULL) }
  res
}
GSEA_RES <- list()
for (n in names(RES)) {
  rk_e <- ranked(RES[[n]], "ensembl"); rk_n <- ranked(RES[[n]], "entrez")
  set.seed(1)
  common <- list(pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500,
                 eps = 0, seed = TRUE, verbose = FALSE)
  GSEA_RES[[n]] <- list(
    Hallmark = run_gsea(GSEA,    c(list(geneList = rk_e, TERM2GENE = HALL$ens), common), paste(n, "Hallmark")),
    GO_BP    = run_gsea(gseGO,   c(list(geneList = rk_e, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP"), common),
                        paste(n, "GO_BP")),
    KEGG     = run_gsea(gseKEGG, c(list(geneList = rk_n, organism = "mmu"), common), paste(n, "KEGG (needs internet)")))
  cat(n, ": Hallmark", nrow0(GSEA_RES[[n]]$Hallmark), "| GO", nrow0(GSEA_RES[[n]]$GO_BP),
      "| KEGG", nrow0(GSEA_RES[[n]]$KEGG), "terms tested\n")
}

clean_terms <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))
core_symbols <- function(core, db) vapply(strsplit(core, "/"), function(g) {
  if (db == "KEGG") paste(na.omit(suppressMessages(mapIds(org.Mm.eg.db, g, "SYMBOL", "ENTREZID"))), collapse = "/")
  else paste(na.omit(sym[gene_ann$gene_id[match(g, gene_ann$ensembl)]]), collapse = "/")
}, character(1))

plot_gsea_lollipop <- function(df, title, path, top_n = 20) {
  df <- df %>% filter(p.adjust < 0.05)
  if (!nrow(df)) { message("No significant terms: ", title); return(invisible(NULL)) }
  df <- df %>% arrange(desc(abs(NES))) %>% head(top_n) %>%
    mutate(Description = factor(Description, levels = Description[order(NES)]),
           Direction = ifelse(NES > 0, "Up in PD1", "Down in PD1"))
  p <- ggplot(df, aes(NES, Description)) +
    geom_segment(aes(x = 0, xend = NES, yend = Description), color = "grey60") +
    geom_point(aes(color = Direction, size = -log10(p.adjust))) +
    scale_color_manual(values = c(`Up in PD1` = "#D7301F", `Down in PD1` = "#2171B5")) +
    geom_vline(xintercept = 0, linetype = "dashed") + theme_minimal(base_size = 11) +
    labs(title = title, x = "NES", y = NULL, size = "-log10 padj")
  save_plot(p, path, 9, 0.3 * nrow(df) + 2)
}

for (n in names(GSEA_RES)) for (db in names(GSEA_RES[[n]])) {
  x <- GSEA_RES[[n]][[db]]; if (nrow0(x) == 0) next
  df <- as.data.frame(x)
  if (db == "Hallmark") df$Description <- clean_terms(df$Description)
  df$core_genes_symbol <- core_symbols(df$core_enrichment, db)
  write.csv(df, file.path(OUTDIR, "enrichment/GSEA", paste0("GSEA_", db, "_", n, ".csv")), row.names = FALSE)
  plot_gsea_lollipop(df, paste("GSEA", db, "—", gsub("_", " ", n)),
                     file.path(OUTDIR, "enrichment/GSEA", paste0("GSEA_", db, "_", n, "_lollipop.pdf")))
}
# running-score plots for the top Hallmark sets
for (n in MAIN) {
  x <- GSEA_RES[[n]]$Hallmark; if (nrow0(x) == 0) next
  top <- as.data.frame(x) %>% filter(p.adjust < 0.05) %>% arrange(desc(abs(NES))) %>% head(4)
  if (nrow(top)) {
    p <- gseaplot2(x, geneSetID = top$ID, title = paste("Hallmark —", gsub("_", " ", n)), pvalue_table = FALSE)
    pdf(file.path(OUTDIR, "enrichment/GSEA", paste0("GSEA_Hallmark_", n, "_running.pdf")), 9, 7); print(p); dev.off(); print(p)
  }
}

# NES comparison heatmap: NT vs sg7 (and interaction), * = padj < 0.05
nes_heatmap <- function(db, path, top_n = 40) {
  cols <- c(MAIN, "Interaction")
  dfs <- lapply(cols, function(n) { x <- GSEA_RES[[n]][[db]]
    if (nrow0(x) == 0) return(NULL)
    as.data.frame(x) %>% transmute(ID, Description, NES, padj = p.adjust, contrast = n) })
  all <- bind_rows(dfs); if (!nrow(all)) return(invisible(NULL))
  keep <- all %>% filter(contrast %in% MAIN, padj < 0.05) %>% group_by(ID) %>%
    summarise(m = max(abs(NES)), .groups = "drop") %>% arrange(desc(m)) %>% head(top_n) %>% pull(ID)
  if (length(keep) < 2) { message("NES heatmap skipped for ", db); return(invisible(NULL)) }
  w  <- all %>% filter(ID %in% keep)
  nm <- w %>% distinct(ID, Description)
  mat <- w %>% select(ID, contrast, NES) %>% pivot_wider(names_from = contrast, values_from = NES) %>%
    column_to_rownames("ID") %>% as.matrix()
  st  <- w %>% mutate(s = ifelse(padj < 0.05, "*", "")) %>% select(ID, contrast, s) %>%
    pivot_wider(names_from = contrast, values_from = s, values_fill = "") %>% column_to_rownames("ID") %>% as.matrix()
  cc <- intersect(cols, colnames(mat)); mat <- mat[, cc, drop = FALSE]; st <- st[rownames(mat), cc, drop = FALSE]
  mat[is.na(mat)] <- 0
  rownames(mat) <- clean_terms(nm$Description[match(rownames(mat), nm$ID)])
  args <- list(mat = mat, display_numbers = st, color = colorRampPalette(c("#2171B5", "white", "#D7301F"))(101),
               breaks = seq(-3, 3, length.out = 102), cluster_cols = FALSE, cluster_rows = TRUE,
               fontsize_row = 7, number_color = "black", main = paste("GSEA NES —", db, "(* padj < 0.05)"))
  do.call(pheatmap, args)
  do.call(pheatmap, c(args, list(filename = path, width = 7, height = max(5, 0.18 * nrow(mat) + 2))))
}
for (db in c("Hallmark", "GO_BP", "KEGG"))
  nes_heatmap(db, file.path(OUTDIR, "enrichment/GSEA", paste0("GSEA_NES_compare_", db, ".pdf")))

save(ORA, ORA_SETS, GSEA_RES, HALL, file = file.path(RDS, "03_enrichment.RData"))
cat("Saved", file.path(RDS, "03_enrichment.RData"), "\n")
