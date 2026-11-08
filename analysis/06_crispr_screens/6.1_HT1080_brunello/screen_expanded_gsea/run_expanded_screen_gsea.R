# Preranked GSEA of HT1080 and A549 screen rankings across Hallmark,
# Reactome, GO Biological Process and GO Cellular Component (MSigDB
# 2026.1.Hs, frozen for the Chapter 3 RNA analyses). Extends the Hallmark-only
# analysis in research/ht1080_hallmark_screen_gsea_20261021 and reuses its
# symbol mapping and common-essential list, so rankings are identical.
#
# Statistic: -normZ (positive NES = sensitizer, guides relatively depleted in
# LOY). Primary FDR: BH across all eligible sets of the four collections,
# separately per ranking; within-collection BH is saved as secondary.
#
#   Rscript run_expanded_screen_gsea.R

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(ggplot2)
})
set.seed(20261021)

proj <- "A:/Dropbox/LOY_paralogs/projects"
thesis <- file.path(proj, "SEPT_11_2026_THESIS")
skill_dir <- "A:/Dropbox/LOY_paralogs/.claude/skills/thesis-figure-style"
hallmark_data <- file.path(thesis, "figures/chapter4/source_data/ht1080_hallmark_screen_gsea_20261021")
out_data <- file.path(thesis, "figures/chapter4/source_data/screen_expanded_gsea_20261021")
out_fig <- file.path(thesis, "figures/chapter4/review")
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)

src <- list(
  msigdb = file.path(proj, "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/resources/msigdb_all_collections.rds"),
  metadata = file.path(proj, "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/resources/msigdb_all_set_metadata.csv"),
  mapping = file.path(hallmark_data, "symbol_mapping_and_scores.csv"),
  essentials = file.path(hallmark_data, "common_essentials_used.csv")
)
stopifnot(all(file.exists(unlist(src))))

# ---- gene sets ------------------------------------------------------------
msig <- readRDS(src$msigdb)
meta <- fread(src$metadata)
collections <- list(
  Hallmark = msig$H,
  Reactome = msig$C2[meta[gs_subcollection == "CP:REACTOME", gs_name]],
  `GO:BP` = msig$C5[meta[gs_subcollection == "GO:BP", gs_name]],
  `GO:CC` = msig$C5[meta[gs_subcollection == "GO:CC", gs_name]]
)
stopifnot(!anyNA(names(unlist(lapply(collections, names)))))
set_collection <- rbindlist(lapply(names(collections), function(k)
  data.table(pathway = names(collections[[k]]), collection = k)))
pathways <- unlist(unname(collections), recursive = FALSE)
stopifnot(!anyDuplicated(names(pathways)))

# ---- rankings (identical to the Hallmark analysis) ------------------------
mapping <- fread(src$mapping)
ess <- unique(na.omit(fread(src$essentials)$symbol))
make_ranks <- function(lab, drop = character()) {
  keep <- mapping[ranking == lab & status %in% c("current", "alias_updated") &
                    !(msigdb_symbol %in% drop)]
  sort(setNames(-keep$normZ, keep$msigdb_symbol), decreasing = TRUE)
}
rankings <- list(
  `HT1080 drugZ` = make_ranks("drugZ"),
  `HT1080 IsogenicZ` = make_ranks("IsogenicZ"),
  `HT1080 drugZ, no CE` = make_ranks("drugZ", ess),
  `HT1080 IsogenicZ, no CE` = make_ranks("IsogenicZ", ess),
  `A549 drugZ` = make_ranks("A549 drugZ")
)
stopifnot(length(rankings$`HT1080 drugZ`) == 18960L,
          length(rankings$`A549 drugZ`) == 18960L)

# ---- GSEA -----------------------------------------------------------------
res <- rbindlist(lapply(names(rankings), function(nm) {
  f <- as.data.table(fgsea(pathways = pathways, stats = rankings[[nm]],
                           minSize = 15, maxSize = 500, gseaParam = 1,
                           eps = 0, nPermSimple = 10000, nproc = 4))
  f[, `:=`(ranking = nm, n_ranked_genes = length(rankings[[nm]]))]
  f[, leadingEdge := vapply(leadingEdge, paste, "", collapse = ";")]
  f
}))
res <- merge(res, set_collection, by = "pathway")
setnames(res, "padj", "FDR_global")
res[, FDR_within_collection := p.adjust(pval, "BH"), by = .(ranking, collection)]
setcolorder(res, c("ranking", "collection", "pathway", "NES", "pval",
                   "FDR_global", "FDR_within_collection", "size"))
setorder(res, ranking, pval)
fwrite(res, file.path(out_data, "expanded_fgsea_all.csv"))

# Wide summary of HT1080 drugZ-significant sets with sensitivity columns.
wide <- dcast(res, collection + pathway + size ~ ranking,
              value.var = c("NES", "FDR_global"))
hits <- res[ranking == "HT1080 drugZ" & FDR_global <= 0.05, pathway]
fwrite(wide[pathway %in% hits][order(`NES_HT1080 drugZ`)],
       file.path(out_data, "ht1080_drugz_significant_sets_wide.csv"))
counts <- res[FDR_global <= 0.05, .(n = .N, sensitizer = sum(NES > 0),
                                    resistor = sum(NES < 0)),
              by = .(ranking, collection)]
fwrite(counts, file.path(out_data, "significant_set_counts.csv"))

# ---- draft figure: HT1080 drugZ only (author decision 21 October 2026; A549
# has no set at FDR <= 0.05); top sets per direction at FDR <= 0.05 ---------
source(file.path(skill_dir, "scripts/thesis_style.R"))
register_thesis_fonts(skill_dir)
plot_rankings <- c(`HT1080 drugZ` = "HT1080")
n_per_dir <- 12
top <- res[ranking == "HT1080 drugZ" & FDR_global <= 0.05]
top <- rbind(top[NES > 0][order(-NES)][seq_len(min(.N, n_per_dir))],
             top[NES < 0][order(NES)][seq_len(min(.N, n_per_dir))])
plot_dt <- res[ranking %in% names(plot_rankings) & pathway %in% top$pathway]
pretty_name <- function(p) {
  p <- sub("^(HALLMARK|REACTOME|GOBP|GOCC)_", "", p)
  p <- gsub("_", " ", p)
  ifelse(nchar(p) > 48, paste0(substr(p, 1, 46), "\u2026"), p)
}
plot_dt[, label := paste0(pretty_name(pathway), "  [", collection, "]")]
row_order <- plot_dt[ranking == "HT1080 drugZ"][order(NES), label]
plot_dt[, label := factor(label, levels = row_order)]
plot_dt[, ranking := factor(ranking, levels = names(plot_rankings))]
plot_dt[, sig := FDR_global <= 0.05]
# Printed statistics, e.g. "NES −3.40   FDR 2.7 × 10⁻³¹".
fmt_fdr <- function(q) {
  sup <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³",
           "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷",
           "8" = "⁸", "9" = "⁹", "-" = "⁻")
  vapply(q, function(x) {
    if (x >= 0.001) return(formatC(x, format = "fg", digits = 2, flag = "#"))
    e <- floor(log10(x)); m <- x / 10^e
    paste0(formatC(m, format = "f", digits = 1), " × 10",
           paste(sup[strsplit(as.character(e), "")[[1]]], collapse = ""))
  }, "")
}
plot_dt[, stat := paste0("NES ", sub("-", "−", sprintf("%.2f", NES)),
                         "   FDR ", fmt_fdr(FDR_global))]
fwrite(plot_dt, file.path(out_data, "figure_values.csv"))

lim <- ceiling(max(abs(plot_dt$NES)) * 2) / 2
p <- ggplot(plot_dt, aes(ranking, label)) +
  geom_point(data = plot_dt[sig == FALSE], colour = "#C8C8C8", size = 0.8) +
  geom_point(data = plot_dt[sig == TRUE], aes(fill = NES, size = abs(NES)),
             shape = 21, colour = "#333333", stroke = 0.25) +
  geom_text(aes(x = 1.3, label = stat), hjust = 0, size = 6.5 / .pt,
            family = "Open Sans", colour = "#222222") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim),
                       breaks = c(-2, 0, 2),
                       labels = c("\u22122\nresistor", "0", "2\nsensitizer"),
                       name = "NES",
                       guide = guide_colourbar(direction = "horizontal",
                                               title.position = "top",
                                               barwidth = unit(70, "pt"),
                                               barheight = unit(7, "pt"),
                                               order = 1)) +
  scale_size_continuous(range = c(1.5, 5), limits = c(1, lim), name = "|NES|",
                        guide = guide_legend(order = 2)) +
  scale_x_discrete(limits = names(plot_rankings), labels = plot_rankings,
                   expand = expansion(add = c(0.35, 3.4))) +
  labs(x = NULL, y = NULL) +
  theme_thesis(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.2),
        legend.title = element_text(face = "plain"),
        plot.margin = margin(10, 16, 2, 2, "pt"))

w <- 6.8; h <- 0.8 + 0.16 * uniqueN(plot_dt$pathway)
stem <- file.path(out_fig, "ch04_supp_screen_expanded_gsea_draft")
ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = cairo_pdf)
ggsave(paste0(stem, ".svg"), p, width = w, height = h, device = svglite::svglite)
ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 600, device = ragg::agg_png)

# ---- provenance -----------------------------------------------------------
jsonlite::write_json(list(
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  script = "SEPT_11_2026_THESIS/research/screen_expanded_gsea_20261021/run_expanded_screen_gsea.R",
  sources = lapply(src, function(f) list(path = f, md5 = unname(tools::md5sum(f)))),
  fgsea_version = as.character(packageVersion("fgsea")),
  collections = vapply(collections, length, 1L),
  settings = list(stat = "-normZ", minSize = 15, maxSize = 500, gseaParam = 1,
                  eps = 0, nPermSimple = 10000, seed = 20261021,
                  fdr_primary = "BH across all eligible sets of the four collections per ranking",
                  fdr_secondary = "BH within collection per ranking"),
  ranked_genes = vapply(rankings, length, 1L),
  tested_sets = res[, .N, by = ranking],
  significant = counts
), file.path(out_data, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
message("done")
