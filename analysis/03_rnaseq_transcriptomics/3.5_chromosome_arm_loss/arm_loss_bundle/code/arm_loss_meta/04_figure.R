#!/usr/bin/env Rscript
# ============================================================================
# Arm-loss meta-analysis figure (thesis style).
#   a  NES heatmap: 14 engineered arm-loss contrasts x watch-list pathways,
#      rows grouped by lost segment (ordered by protein-coding genes lost);
#      dot = FDR < 0.05.
#   b  profile agreement: Spearman rho of Hallmark+GO:BP NES profiles, arm loss vs
#      arm loss, and arm loss vs p53 induction in the monosomic backgrounds
#      (the diploid ip53 contrast is excluded: induction is inert there).
#   c  p53 pathway NES vs protein-coding genes lost; open points = contrasts,
#      filled = segment mean; statistic at the SEGMENT level (n = 7).
#   d  RPTEC 3p-loss clones, day 18 vs day 90 (GSE336576): the acute -> chronic flip.
# Inputs: results/gsea/*_GSEA.tsv.gz, results/meta_*.tsv (03_meta.R).
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
source("/global/home/hpc5173/projects/LOY/scRNA/brain_atlas/scripts/figures/00_fig_theme.R")
ROOT    <- "/global/home/hpc5173/projects/LOY/RNA/arm_loss_meta"
RES     <- file.path(ROOT, "results")
FIG_DIR <- file.path(ROOT, "figures"); SRC_DIR <- file.path(FIG_DIR, "source_data")
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
STEM <- "main_arm_loss_meta"

G <- rbindlist(lapply(list.files(file.path(RES, "gsea"), "\\.tsv\\.gz$", full.names = TRUE), zc), fill = TRUE)
G[, leadingEdge := NULL]
SC <- fread(file.path(RES, "meta_scaling.tsv"))            # contrast, cis, n_lost, ...
SC[, seg := sub("^(arm|chrom):", "", cis)]
SC[seg == "chr10,chr18", seg := "10+18"][seg == "chr13", seg := "13"][seg == "chrY", seg := "Y*"]

LAB <- c(ARPE19_chrY = "ARPE-19",
         A2058_1q = "A2058", A2780_1q_vsWT = "A2780 (vs WT)", A2780_1q_vsTK = "A2780 (vs TK)",
         AGS_1q_vsTK = "AGS (vs TK)", MCF10A_1q = "MCF10A", MCF10Ahr_1q = "MCF10A-HRAS",
         A2058_7p = "A2058", A2058_8q = "A2058", HCT116_8q = "HCT116",
         RPTEC_3p_early = "RPTEC day 18", RPTEC_3p_late = "RPTEC day 90",
         RPE1_mono1018 = "RPE1", RPE1_mono13 = "RPE1", RPE1_mono19p = "RPE1")

# ---- a: heatmap ------------------------------------------------------------------
PW <- c(HALLMARK_P53_PATHWAY = "p53 pathway", HALLMARK_E2F_TARGETS = "E2F targets",
        HALLMARK_G2M_CHECKPOINT = "G2M checkpoint", HALLMARK_MYC_TARGETS_V1 = "MYC targets V1",
        HALLMARK_OXIDATIVE_PHOSPHORYLATION = "OXPHOS", RIBOSOMAL_PROTEIN = "Ribosomal proteins",
        HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-α response",
        HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-γ response",
        HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFα–NF-κB",
        HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
        HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT")
A <- merge(G[type == "arm_loss" & pathway %in% names(PW), .(contrast, pathway, NES, padj)],
           SC[, .(contrast, seg, n_lost)], by = "contrast")
A[, pw := factor(PW[pathway], PW)]
A[, seg_lab := sprintf("%s (%d)", seg, n_lost)]
seg_order <- unique(SC[order(n_lost, seg), sprintf("%s (%d)", seg, n_lost)])
A[, seg_lab := factor(seg_lab, rev(seg_order))]
A[, row := factor(LAB[contrast], unique(LAB))]
A[, sig := padj < 0.05]
save_source(A[, .(contrast, segment = seg, n_lost, pathway, NES, padj)], STEM, "a")
lim <- max(abs(A$NES))
pa <- ggplot(A, aes(pw, row, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_point(data = A[sig == TRUE], size = 0.7, colour = "#222222", show.legend = FALSE) +
  facet_grid(seg_lab ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = NPG[["blue"]], mid = "white", high = NPG[["vermilion"]],
                       limits = c(-lim, lim), name = "NES") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = SMALL_PT),
        axis.text.y = element_text(size = SMALL_PT), axis.line = element_blank(), axis.ticks = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = SMALL_PT),
        strip.background = element_blank(), panel.spacing.y = unit(1.5, "pt"),
        legend.position = "right", legend.key.height = unit(16, "pt"))

# ---- b: profile agreement --------------------------------------------------------
aa <- fread(file.path(RES, "meta_arm_vs_arm_rho.tsv"))[, .(cmp = "arm vs\narm", rho)]
ap <- fread(file.path(RES, "meta_arm_vs_p53_rho.tsv"))[p53 != "ip53_diploid", .(cmp = "arm vs\np53", rho)]
B <- rbind(aa, ap); B[, cmp := factor(cmp, unique(cmp))]
save_source(B, STEM, "b")
med <- B[, .(m = median(rho), n = .N), by = cmp]
for (k in seq_len(nrow(med))) record_num(sprintf("b_median_rho_%d", k), med$m[k], n = med$n[k], note = gsub("\n", " ", as.character(med$cmp[k])))
pb <- ggplot(B, aes(cmp, rho)) +
  geom_hline(yintercept = 0, linewidth = LINE_W, colour = "grey70", linetype = "22") +
  geom_point(position = position_jitter(width = 0.18, height = 0, seed = 1), size = 0.9, shape = 16, colour = "grey55") +
  geom_errorbar(data = med, aes(y = m, ymin = m, ymax = m), width = 0.5, linewidth = 0.6, colour = "#222222") +
  geom_text(data = med, aes(y = 0.95, label = sprintf("%.2f\nn=%d", m, n)), size = 6.5 / .pt, lineheight = 0.9) +
  scale_y_continuous(limits = c(-0.6, 1.05), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = NULL, y = expression("NES profile "*rho))

# ---- c: p53 vs genes lost ---------------------------------------------------------
Cc <- merge(G[type == "arm_loss" & pathway == "HALLMARK_P53_PATHWAY", .(contrast, NES)],
            SC[, .(contrast, seg, n_lost)], by = "contrast")
Cm <- Cc[, .(NES = mean(NES), k = .N), by = .(seg, n_lost)]
# chrY is hemizygous (1 copy -> 0), not a halving like every other segment here,
# so it is plotted for reference but excluded from the dose-response fit.
Cm[, hemi := seg == "Y*"]
ct <- cor.test(Cm[hemi == FALSE, n_lost], Cm[hemi == FALSE, NES], method = "spearman", exact = TRUE)
record_num("c_segment_rho", ct$estimate, n = nrow(Cm[hemi == FALSE]),
           note = sprintf("P = %.3f; chrY excluded (hemizygous)", ct$p.value))
save_source(Cc, STEM, "c")
pc <- ggplot(Cc, aes(n_lost, NES)) +
  geom_hline(yintercept = 0, linewidth = LINE_W, colour = "grey70", linetype = "22") +
  geom_point(shape = 1, size = 1.5, colour = "grey45", stroke = 0.4) +
  geom_point(data = Cm[hemi == FALSE], size = 2.2, colour = NPG[["blue"]]) +
  geom_point(data = Cm[hemi == TRUE], size = 2.2, shape = 17, colour = "grey45") +
  ggrepel::geom_text_repel(data = Cm, aes(label = seg), size = 6.5 / .pt, family = FONT, seed = 1,
                           box.padding = 0.35, min.segment.length = 0.2, segment.size = 0.2, segment.colour = "grey60") +
  labs(title = sprintf("segment level: ρ = %+.2f, P = %.2f (n = %d)\n1q: p53 up in 6/6 backgrounds; ▲ Y* hemizygous, not fitted",
                       ct$estimate, ct$p.value, nrow(Cm[hemi == FALSE]))) +
  scale_y_continuous(limits = c(-2.6, 2.6)) +
  scale_x_continuous(breaks = c(300, 600, 900)) +
  labs(x = "Protein-coding genes lost", y = "p53 pathway NES") +
  theme(plot.title = element_text(size = 6.5, face = "plain", lineheight = 0.95))

# ---- d: RPTEC day 18 -> day 90 ---------------------------------------------------------
DP <- c(HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-α response", HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-γ response",
        HALLMARK_E2F_TARGETS = "E2F targets", HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
        HALLMARK_MYC_TARGETS_V1 = "MYC targets V1", HALLMARK_P53_PATHWAY = "p53 pathway",
        RIBOSOMAL_PROTEIN = "Ribosomal proteins")
D <- G[contrast %in% c("RPTEC_3p_early", "RPTEC_3p_late") & pathway %in% names(DP), .(contrast, pathway, NES, padj)]
D[, day := fifelse(contrast == "RPTEC_3p_early", "day 18", "day 90")]
D[, pw := factor(DP[pathway], rev(DP))]
save_source(D, STEM, "d")
Dw <- dcast(D, pw ~ day, value.var = "NES")
pd <- ggplot(D, aes(NES, pw)) +
  geom_vline(xintercept = 0, linewidth = LINE_W, colour = "grey70") +
  geom_segment(data = Dw, aes(x = `day 18`, xend = `day 90`, y = pw, yend = pw),
               linewidth = 0.4, colour = "grey60",
               arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  geom_point(aes(colour = day, shape = padj < 0.05), size = 1.8) +
  scale_colour_manual(values = c(`day 18` = NPG[["teal"]], `day 90` = NPG[["brown"]]), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
  labs(x = "NES vs unedited", y = NULL, title = "RPTEC 3p-loss clones") +
  theme(plot.title = element_text(size = SMALL_PT, face = "plain"), legend.position = "bottom",
        legend.box = "vertical", legend.spacing.y = unit(0, "pt"))

bottom <- (pb | pc | pd) + plot_layout(widths = c(0.7, 1, 1.15))
p <- (pa / free(bottom)) + plot_layout(heights = c(1.35, 1)) + tags()
save_fig(p, STEM, h_mm = 190)
message("done")
