#!/usr/bin/env Rscript
# GSE336056: 3p-loss vs Yq-null pathway response, NES against NES (thesis style).
#   a  Hallmark (49 sets): NES in sg3p 3p-loss vs sgYq Yq-null; shared + interferon labelled.
#   b  GO:BP (3,178 sets): same axes, opposite-direction terms labelled.
# Class = FDR < 0.05 in both arms (same / opposite sign) or in one arm only.
# Inputs: de/GSE336056_3p_vs_Yq_pathway_overlap.tsv (compare_3p_yq_pathways.R).
suppressPackageStartupMessages({ library(data.table) })
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
source(file.path(ROOT, "brain_atlas/scripts/figures/00_fig_theme.R"))
FIG_DIR <- file.path(ROOT, "analysis/GSE336056/figures"); SRC_DIR <- file.path(FIG_DIR, "source_data")
G <- fread(file.path(ROOT, "analysis/GSE336056/de/GSE336056_3p_vs_Yq_pathway_overlap.tsv"))
STEM <- "main_gse336056_3p_vs_yq_nes"

CLS <- c(shared = "both, same direction", opposite = "both, opposite",
         `3p_only` = "3p-loss only", Yq_only = "Yq-null only", neither = "neither")
PAL_C <- setNames(c(NPG[["blue"]], NPG[["vermilion"]], NPG[["teal"]], NPG[["salmon"]], "#C8C8C8"), CLS)
G[, cls := factor(CLS[class], CLS)]
rho <- G[, .(rho = cor(NES_3p, NES_Yq, method = "spearman"), n = .N), by = collection]
for (k in rho$collection) record_num(paste0("rho_", k), rho[collection == k, rho], n = rho[collection == k, n])

nes_axes <- function(lim) list(
  geom_abline(slope = 1, intercept = 0, linewidth = LINE_W, colour = "grey70", linetype = "22"),
  geom_hline(yintercept = 0, linewidth = LINE_W, colour = "grey85"),
  geom_vline(xintercept = 0, linewidth = LINE_W, colour = "grey85"),
  scale_colour_manual(values = PAL_C, name = "FDR < 0.05 in", drop = FALSE,
                      guide = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1))),
  coord_equal(xlim = lim, ylim = lim),
  labs(x = "NES, 3p-loss vs intact\n(649 vs 8,934 cells)", y = "NES, Yq-null vs intact\n(194 vs 6,407 cells)"))
stat_lab <- function(k) {
  r <- rho[collection == k]; x <- G[collection == k]
  sprintf("Spearman ρ = %.2f (%s sets)\n%d shared, %d opposite", r$rho, format(r$n, big.mark = ","),
          sum(x$class == "shared"), sum(x$class == "opposite"))
}

# ---- a: Hallmark ------------------------------------------------------------------------
H <- G[collection == "H"]
NAMES <- c(E2F_TARGETS = "E2F targets", G2M_CHECKPOINT = "G2M checkpoint", MYC_TARGETS_V1 = "MYC V1", MYC_TARGETS_V2 = "MYC V2",
  MITOTIC_SPINDLE = "Mitotic spindle", P53_PATHWAY = "p53 pathway", EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  TNFA_SIGNALING_VIA_NFKB = "TNFα–NF-κB", APOPTOSIS = "Apoptosis", COMPLEMENT = "Complement",
  INTERFERON_ALPHA_RESPONSE = "IFN-α response", INTERFERON_GAMMA_RESPONSE = "IFN-γ response",
  OXIDATIVE_PHOSPHORYLATION = "OXPHOS", MTORC1_SIGNALING = "mTORC1", UNFOLDED_PROTEIN_RESPONSE = "UPR", HYPOXIA = "Hypoxia")
H[, short := sub("^HALLMARK_", "", pathway)]
H[, lab := fifelse(short %in% names(NAMES), NAMES[short], "")]
save_source(H[, .(pathway, NES_3p, padj_3p, NES_Yq, padj_Yq, class)], STEM, "a")
lim_h <- range(c(H$NES_3p, H$NES_Yq)) + c(-0.3, 0.3)
HU <- H[lab != "" & NES_3p > 0][order(-NES_Yq)]
HU[, ly := seq(lim_h[2] - 0.15, 0.45, length.out = .N)]
# negative side: fixed positions (x, y, hjust); segment starts at the text edge nearest the point
POS <- data.table(lab = c("E2F targets", "G2M checkpoint", "MYC V1", "MYC V2", "Mitotic spindle", "mTORC1", "UPR"),
                  lx  = c(-3.25, -3.25, -1.90, -1.85, -3.25, -0.55, -0.55),
                  ly  = c(-3.35, -2.35, -3.35, -1.95, -1.40, -1.35, -2.05),
                  sx  = c(-3.00, -2.90, -1.95, -1.90, -2.40, -0.60, -0.60),   # segment start (text edge)
                  sy  = c(-3.25, -2.45, -3.30, -1.95, -1.48, -1.35, -2.05),
                  hj  = 0)
HD <- merge(H[lab != "" & NES_3p < 0], POS, by = "lab")
pa <- ggplot(H, aes(NES_3p, NES_Yq, colour = cls)) + nes_axes(lim_h) +
  geom_point(size = 1.6, shape = 16, show.legend = TRUE) +
  geom_segment(data = HD, aes(x = sx, y = sy, xend = NES_3p, yend = NES_Yq), linewidth = 0.2, colour = "grey60") +
  geom_text(data = HD, aes(x = lx, y = ly, label = lab, hjust = hj), size = 6.5 / .pt, family = FONT, show.legend = FALSE) +
  geom_segment(data = HU, aes(x = 1.30, y = ly, xend = NES_3p, yend = NES_Yq), linewidth = 0.2, colour = "grey60") +
  geom_text(data = HU, aes(x = 1.25, y = ly, label = lab), hjust = 1, size = 6.5 / .pt, family = FONT, show.legend = FALSE) +
  annotate("text", x = lim_h[2], y = lim_h[1], hjust = 1, vjust = 0, size = SMALL_PT / .pt, lineheight = 0.95, label = stat_lab("H")) +
  ggtitle("Hallmark")

# ---- b: GO:BP ---------------------------------------------------------------------------
B <- G[collection == "GOBP"]
B[, lab := fifelse(class == "opposite", gsub("type i ", "type I ", gsub(" t cell", " T cell", tolower(gsub("_", " ", sub("^GOBP_", "", pathway))))), "")]
save_source(B[, .(pathway, NES_3p, padj_3p, NES_Yq, padj_Yq, class)], STEM, "b")
lim_b <- range(c(B$NES_3p, B$NES_Yq)) + c(-0.3, 0.3)
pb <- ggplot(B, aes(NES_3p, NES_Yq, colour = cls)) + nes_axes(lim_b) +
  rast(geom_point(data = B[class == "neither"], size = 0.5, alpha = 0.6, shape = 16, show.legend = FALSE)) +
  geom_point(data = B[class %in% c("3p_only", "Yq_only")], size = 0.7, alpha = 0.8, shape = 16, show.legend = FALSE) +
  geom_point(data = B[class %in% c("shared", "opposite")], size = 0.9, shape = 16, show.legend = FALSE) +
  ggrepel::geom_label_repel(data = B[lab != ""], aes(label = lab), size = 6.5 / .pt, family = FONT, fill = alpha("white", 0.85),
                           label.size = NA, label.padding = unit(0.1, "lines"),
                           segment.size = 0.2, segment.colour = "grey60", box.padding = 0.5, min.segment.length = 0,
                           nudge_y = 0.4, seed = 1, show.legend = FALSE) +
  annotate("text", x = lim_b[2], y = lim_b[1], hjust = 1, vjust = 0, size = SMALL_PT / .pt, lineheight = 0.95, label = stat_lab("GOBP")) +
  ggtitle("GO biological process")

p <- (pa | pb) + plot_layout(guides = "collect") + tags() & theme(legend.position = "bottom")
save_fig(p, STEM, h_mm = 105)
message("done")
