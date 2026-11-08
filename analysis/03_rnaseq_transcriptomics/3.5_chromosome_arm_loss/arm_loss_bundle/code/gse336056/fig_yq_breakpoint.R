#!/usr/bin/env Rscript
# ============================================================================
# GSE336056: where does the sgYq-induced chrY loss start?
# Thesis-style figure (thesis-figure-style skill via brain_atlas 00_fig_theme.R):
#   a  per-sample % of cells that are Yq-null while retaining Yp (Wilson CI),
#      Fisher P sgYq vs sgNTC per day
#   b  per-cell breakpoint map: for every expressed MSY gene (hg38 order),
#      detection in Yq-null sgYq cells relative to depth-matched Yq-intact sgYq
#      cells. Filled = genes used to define Yq-null (trivially 0 on Yq), open =
#      HELD-OUT genes, which are the real test of a segmental loss.
#   c  pseudobulk log2FC (sgYq / same-lane sgNTC) per MSY gene vs coordinate,
#      per day, bootstrap CI over cells.
# Cells >= 10k UMI throughout (the D10 lane is ~22k median, D5 ~40k).
# Outputs: analysis/GSE336056/figures/main_gse336056_yq_breakpoint.{pdf,svg,png}
#          + source_data/ per panel.
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
source(file.path(ROOT, "brain_atlas/scripts/figures/00_fig_theme.R"))
FIG_DIR <- file.path(ROOT, "analysis/GSE336056/figures"); SRC_DIR <- file.path(FIG_DIR, "source_data")
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
STAGE <- file.path(ROOT, "work/GSE336056_staged")
set.seed(1)
MIN_UMI <- 10000; CEN <- c(10.32, 10.54)          # hg38 chrY centromere (Mb)

# ---- MSY gene coordinates (GENCODE v49, PAR excluded) ---------------------------
source(file.path(ROOT, "scripts/resolve_hg38_msy_genes.R"), local = TRUE)
first <- readRDS(list.files(STAGE, "\\.rds$", full.names = TRUE)[1])
msy <- as.data.table(resolve_hg38_msy_genes(dataset_genes = rownames(first),
         gtf_path = file.path(ROOT, "metadata/references/gencode.v49.annotation.gtf.gz")))
msy <- unique(msy[present_in_dataset == TRUE, .(gene = matched_feature, gene_type, start, end)], by = "gene")
msy[, `:=`(mb = start / 1e6, arm = fifelse(start / 1e6 < CEN[1], "Yp", "Yq"))]
setorder(msy, start)
# selection panel (same 14 genes as yq_arm_check.R)
SEL <- c("RPS4Y1","ZFY","PCDH11Y","AMELY","TBL1Y","TSPY1",
         "USP9Y","DDX3Y","UTY","TMSB4Y","NLGN4Y","KDM5D","EIF1AY","RPS4Y2")
msy[, selection := gene %in% SEL]
message(sprintf("MSY genes in dataset: %d (%d Yp, %d Yq)", nrow(msy), sum(msy$arm=="Yp"), sum(msy$arm=="Yq")))

# ---- per-cell table ----------------------------------------------------------------
cells <- list(); Y <- list()
for (f in list.files(STAGE, "\\.rds$", full.names = TRUE)) {
  sid <- sub("\\.rds$", "", basename(f)); m <- readRDS(f)
  n <- Matrix::colSums(m); keep <- n >= MIN_UMI; m <- m[, keep, drop = FALSE]; n <- n[keep]
  g <- intersect(msy$gene, rownames(m)); ym <- m[g, , drop = FALSE]
  p <- tstrsplit(sid, "_")
  yp_sel <- Matrix::colSums(ym[intersect(g, SEL[SEL %in% msy[arm=="Yp", gene]]), , drop = FALSE])
  yq_sel <- Matrix::colSums(ym[intersect(g, SEL[SEL %in% msy[arm=="Yq", gene]]), , drop = FALSE])
  cells[[sid]] <- data.table(cell = colnames(m), sample = sid, arm_grp = p[[1]], day = p[[2]], hto = p[[3]],
                             nUMI = n, yp_sel = yp_sel, yq_sel = yq_sel)
  Y[[sid]] <- ym
}
D <- rbindlist(cells); Ymat <- do.call(cbind, Y); stopifnot(identical(colnames(Ymat), D$cell))
D[, day := factor(day, c("D5", "D10"))]
D[, arm_grp := factor(arm_grp, c("sgNTC", "sg3p", "sgYq"))]
D[, status := fifelse(yp_sel >= 5 & yq_sel == 0, "Yq_null",
              fifelse(yp_sel >= 5 & yq_sel >= 2, "Yq_intact", "ambiguous"))]
PAL_ARM <- c(sgNTC = "grey55", sg3p = NPG[["lavender"]], sgYq = NPG[["teal"]])
PAL_Y   <- c(Yp = PAL[["Y_intact"]], Yq = PAL[["LOY"]])          # retained arm / lost arm

# ---- a: % Yq-null (of Yp-retained cells) per sample -----------------------------
A <- D[yp_sel >= 5, .(n = .N, k = sum(status == "Yq_null")), by = .(sample, arm_grp, day, hto)]
A[, pct := 100 * k / n]; A[, c("lo", "hi") := wilson(k, n)]
fish <- D[yp_sel >= 5 & arm_grp %in% c("sgYq", "sgNTC"), {
  t <- table(factor(arm_grp, c("sgYq","sgNTC")), factor(status == "Yq_null", c(TRUE, FALSE)))
  ft <- fisher.test(t); .(OR = unname(ft$estimate), p = ft$p.value,
                          pct_yq = 100 * t["sgYq","TRUE"]/sum(t["sgYq",]), pct_ntc = 100 * t["sgNTC","TRUE"]/sum(t["sgNTC",]))
}, by = day]
setorder(fish, day); print(fish)
for (i in seq_len(nrow(fish))) record_num(sprintf("a_fisher_OR_%s", fish$day[i]), fish$OR[i], note = sprintf("P=%s", fmt_p(fish$p[i])))
fish[, lab := sprintf('atop("OR = %.1f", %s)', OR, sapply(p, p_expr))]
save_source(A, "main_gse336056_yq_breakpoint", "a")
A[, x := paste0(arm_grp, "\n", hto)]
pa <- ggplot(A, aes(x = arm_grp, y = pct, colour = arm_grp, group = sample)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = LINE_W * 1.5, position = position_dodge(0.5)) +
  geom_point(size = 1.8, position = position_dodge(0.5)) +
  facet_grid(. ~ day, labeller = labeller(day = c(D5 = "Day 5", D10 = "Day 10"))) +
  scale_colour_manual(values = PAL_ARM, guide = "none") +
  scale_y_continuous(limits = c(0, max(A$hi) * 1.35), expand = expansion(0)) +
  geom_text(data = fish, inherit.aes = FALSE, parse = TRUE, hjust = 0, vjust = 1, size = SMALL_PT / .pt,
            aes(x = 0.6, y = max(A$hi) * 1.33, label = lab)) +
  scale_x_discrete(labels = c(sgNTC = "NTC", sg3p = "3p", sgYq = "Yq"), expand = expansion(add = 0.6)) +
  labs(x = "sgRNA", y = "Yq-null cells, Yp retained (%)")

# ---- b: per-cell breakpoint map (held-out genes) --------------------------------
yq <- D[arm_grp == "sgYq"]
nul <- yq[status == "Yq_null"]; int <- yq[status == "Yq_intact"]
# depth-match: for each Yq-null cell take the nearest-depth unused Yq-intact cell of the same sample
matched <- character()
for (s in unique(nul$sample)) {
  pool <- int[sample == s]; pool <- pool[order(nUMI)]; used <- logical(nrow(pool))
  for (u in nul[sample == s, nUMI]) {
    j <- which.min(abs(pool$nUMI - u) + ifelse(used, 1e12, 0)); used[j] <- TRUE; matched <- c(matched, pool$cell[j])
  }
}
det <- function(cellset) Matrix::rowMeans(Ymat[, cellset, drop = FALSE] > 0)
B <- data.table(gene = rownames(Ymat), det_null = det(nul$cell), det_int = det(matched),
                k_null = Matrix::rowSums(Ymat[, nul$cell, drop = FALSE] > 0))
B <- merge(B, msy, by = "gene")[det_int >= 0.05]
B[, ratio := det_null / det_int]
B[, c("lo", "hi") := { w <- wilson(k_null, nrow(nul)); .(w$lo / 100 / det_int, w$hi / 100 / det_int) }]
setorder(B, start)
brk <- c(B[arm == "Yp" & ratio > 0.5, max(end)] / 1e6, B[arm == "Yq" & ratio < 0.5, min(start)] / 1e6)
message(sprintf("b: %d Yq-null vs %d matched Yq-intact sgYq cells; %d genes; break interval %.2f-%.2f Mb",
                nrow(nul), length(matched), nrow(B), brk[1], brk[2]))
record_num("b_n_yq_null", nrow(nul)); record_num("b_break_lo_mb", brk[1]); record_num("b_break_hi_mb", brk[2])
held <- B[selection == FALSE & arm == "Yq"]
record_num("b_heldout_yq_median_ratio", median(held$ratio), n = nrow(held))
save_source(B[, .(gene, gene_type, arm, selection, start, end, mb, det_null, det_int, ratio, lo, hi)],
            "main_gse336056_yq_breakpoint", "b")
pb <- ggplot(B, aes(x = mb, y = ratio, colour = arm)) +
  annotate("rect", xmin = CEN[1], xmax = CEN[2], ymin = -Inf, ymax = Inf, fill = "grey85", colour = NA) +
  annotate("rect", xmin = brk[1], xmax = brk[2], ymin = -Inf, ymax = Inf, fill = PAL[["LOY"]], alpha = 0.08, colour = NA) +
  geom_hline(yintercept = c(0, 1), linewidth = LINE_W, colour = "grey70", linetype = c("solid", "22")) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = LINE_W) +
  geom_point(aes(shape = selection), size = 1.9, stroke = 0.6) +
  ggrepel::geom_text_repel(aes(label = gene), size = 6.5 / .pt, family = FONT, max.overlaps = 30,
                           segment.size = 0.2, segment.colour = "grey60", box.padding = 0.25, seed = 1,
                           show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = NULL,
                     labels = c(`TRUE` = "selection gene", `FALSE` = "held-out gene"),
                     breaks = c("TRUE", "FALSE")) +
  scale_colour_manual(values = PAL_Y, name = NULL, labels = c(Yp = "Yp (retained)", Yq = "Yq (lost)")) +
  scale_x_continuous(breaks = seq(0, 25, 5), limits = c(2, 23)) +
  scale_y_continuous(limits = c(-0.05, max(1.6, max(B$hi) * 1.05)), breaks = seq(0, 1.5, 0.5)) +
  annotate("text", x = CEN[2] + 0.3, y = 1.58, hjust = 0, vjust = 1, size = SMALL_PT / .pt, lineheight = 0.95,
           label = sprintf("break interval %.1f\u2013%.1f Mb\n(centromere %.1f\u2013%.1f Mb)\n%d Yq-null vs %d depth-matched cells",
                           brk[1], brk[2], CEN[1], CEN[2], nrow(nul), length(matched))) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)), shape = guide_legend(order = 2)) +
  labs(x = "chrY position (Mb, hg38)", y = "Detection in Yq-null cells\n(relative to matched Yq-intact)") +
  inside(0.60, 0.66) + theme(legend.spacing.y = unit(0, "pt"), legend.key.height = unit(8, "pt"))

# ---- c: per-cell Yp vs Yq expression, sgNTC vs sgYq ------------------------------
Cc <- D[arm_grp %in% c("sgNTC", "sgYq"), .(cell, arm_grp, day, nUMI, status,
        yp_per10k = 1e4 * yp_sel / nUMI, yq_per10k = 1e4 * yq_sel / nUMI)]
Cc[, status := factor(status, c("Yq_intact", "Yq_null", "ambiguous"))]
cnt <- Cc[, .(n = .N, k = sum(status == "Yq_null")), by = arm_grp][, lab := sprintf("%d / %s cells Yq-null (%.1f%%)", k, format(n, big.mark = ","), 100 * k / n)]
save_source(Cc, "main_gse336056_yq_breakpoint", "c")
PAL_ST <- c(Yq_intact = PAL[["Y_intact"]], Yq_null = PAL[["LOY"]], ambiguous = "#BDBDBD")
pc <- ggplot(Cc[order(status == "Yq_null")], aes(x = yp_per10k, y = yq_per10k, colour = status)) +
  rast(geom_point(size = 0.35, alpha = 0.55, shape = 16)) +
  facet_grid(. ~ arm_grp, labeller = labeller(arm_grp = c(sgNTC = "sgNTC (non-targeting)", sgYq = "sgYq (chrYq break)"))) +
  geom_text(data = cnt, inherit.aes = FALSE, aes(x = 0.15, y = Inf, label = lab), hjust = 0, vjust = 1.4, size = SMALL_PT / .pt) +
  scale_colour_manual(values = PAL_ST, name = NULL,
                      labels = c(Yq_intact = "Yq intact", Yq_null = "Yq-null, Yp retained", ambiguous = "ambiguous (Yq = 1 or Yp < 5)")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Yp UMI per 10k (RPS4Y1, ZFY, TBL1Y, PCDH11Y, AMELY, TSPY1)", y = "Yq UMI per 10k\n(USP9Y \u2026 EIF1AY, 8 genes)") +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- assemble ---------------------------------------------------------------------
p <- ((pa | pb) + plot_layout(widths = c(0.36, 0.64))) / pc + plot_layout(heights = c(0.55, 0.45)) + tags()
save_fig(p, "main_gse336056_yq_breakpoint", h_mm = 150)
message("done")
