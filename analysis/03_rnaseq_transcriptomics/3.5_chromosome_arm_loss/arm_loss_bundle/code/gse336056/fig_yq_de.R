#!/usr/bin/env Rscript
# GSE336056 sgYq arm: Yq-null vs Yq-intact DE figure (thesis style).
#   a  volcano, MAST pooled fit (autosomal + chrX); p53 targets and proliferation
#      genes labelled; chrY genes have no finite MAST logFC (complete separation)
#      and are listed in the caption instead.
#   b  Hallmark GSEA (padj < 0.05), NES bars, ranked by MAST logFC.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
source(file.path(ROOT, "brain_atlas/scripts/figures/00_fig_theme.R"))
FIG_DIR <- file.path(ROOT, "analysis/GSE336056/figures"); SRC_DIR <- file.path(FIG_DIR, "source_data")
DE <- file.path(ROOT, "analysis/GSE336056/de")
res  <- fread(file.path(DE, "GSE336056_sgYq_YqNull_vs_intact_MAST.tsv"))
gsea <- fread(file.path(DE, "GSE336056_sgYq_YqNull_vs_intact_GSEA.tsv"))

# ---- a: volcano ----------------------------------------------------------------------
V <- res[gene_class %in% c("autosomal", "chrX") & is.finite(logFC_pooled)]
V[, nlp := -log10(pmax(p_pooled, 1e-300))]
V[, sig := fifelse(padj_pooled < 0.05 & logFC_pooled > 0, "up", fifelse(padj_pooled < 0.05 & logFC_pooled < 0, "down", "ns"))]
P53   <- c("CDKN1A","MDM2","BAX","RPS27L","ZMAT3","PHLDA3","DDB2","RRM2B","GDF15","MIR34AHG","TRIM22","SESN1","TP53I3","FDXR","AEN","TNFRSF10B","BBC3","PLK3","SPATA18","CCNG1")
PROLIF <- c("HMGB1","TUBA1B","SMC3","PTMA","RAN","HNRNPD","HNRNPF","LRR1","TCF19","HSP90AA1","PSMB1","NUP54","PRELID1")
V[, lab := fifelse(gene %in% c(P53, PROLIF) & padj_pooled < 0.05, gene, "")]
V[, cls := fifelse(gene %in% P53 & padj_pooled < 0.05, "p53 target", fifelse(sig == "ns", "ns", "other significant"))]
n_up <- sum(V$sig == "up"); n_dn <- sum(V$sig == "down")
record_num("a_n_up", n_up); record_num("a_n_down", n_dn)
save_source(V[, .(gene, gene_class, logFC_pooled, p_pooled, padj_pooled, logFC_D5, logFC_D10, sig, cls)], "main_gse336056_yq_de", "a")
PAL_V <- c(`p53 target` = NPG[["red"]], `other significant` = "#222222", ns = "#C8C8C8")
pa <- ggplot(V, aes(x = logFC_pooled, y = nlp, colour = cls)) +
  rast(geom_point(data = V[cls == "ns"], size = 0.5, alpha = 0.6, shape = 16)) +
  geom_point(data = V[cls != "ns"], size = 1.2, shape = 16) +
  ggrepel::geom_text_repel(data = V[lab != ""], aes(label = lab), size = 6.5 / .pt, family = FONT, max.overlaps = 40,
                           segment.size = 0.2, segment.colour = "grey60", box.padding = 0.3, seed = 1, show.legend = FALSE) +
  geom_hline(yintercept = -log10(V[sig != "ns", max(p_pooled)]), linewidth = LINE_W, colour = "grey70", linetype = "22") +
  geom_vline(xintercept = 0, linewidth = LINE_W, colour = "grey70") +
  scale_colour_manual(values = PAL_V, name = NULL, breaks = c("p53 target", "other significant", "ns"),
                      labels = c("p53 target (FDR < 0.05)", "other FDR < 0.05", "not significant")) +
  annotate("text", x = min(V$logFC_pooled), y = max(V$nlp), hjust = 0, vjust = 1, size = SMALL_PT / .pt, lineheight = 0.95,
           label = sprintf("194 Yq-null vs 6,407 intact\n%d up, %d down (FDR < 0.05)\nMAST ~ status + CDR + library", n_up, n_dn)) +
  labs(x = expression("MAST log"[2]~"FC, Yq-null vs Yq-intact"), y = expression(-log[10]~italic(P))) +
  inside(0.80, 0.74)

# ---- b: Hallmark GSEA ---------------------------------------------------------------------
G <- gsea[collection == "H" & padj < 0.05]
NAMES <- c(HALLMARK_E2F_TARGETS = "E2F targets", HALLMARK_MYC_TARGETS_V1 = "MYC targets V1", HALLMARK_MYC_TARGETS_V2 = "MYC targets V2",
  HALLMARK_G2M_CHECKPOINT = "G2M checkpoint", HALLMARK_P53_PATHWAY = "p53 pathway",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "Epithelial\u2013mesenchymal transition",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF\u03b1 signalling via NF-\u03baB", HALLMARK_COMPLEMENT = "Complement",
  HALLMARK_MITOTIC_SPINDLE = "Mitotic spindle", HALLMARK_UNFOLDED_PROTEIN_RESPONSE = "Unfolded protein response",
  HALLMARK_MTORC1_SIGNALING = "mTORC1 signalling", HALLMARK_COAGULATION = "Coagulation", HALLMARK_APOPTOSIS = "Apoptosis",
  HALLMARK_HYPOXIA = "Hypoxia", HALLMARK_ANGIOGENESIS = "Angiogenesis", HALLMARK_KRAS_SIGNALING_UP = "KRAS signalling up",
  HALLMARK_KRAS_SIGNALING_DN = "KRAS signalling down", HALLMARK_HEME_METABOLISM = "Heme metabolism",
  HALLMARK_ALLOGRAFT_REJECTION = "Allograft rejection", HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Oxidative phosphorylation", HALLMARK_IL2_STAT5_SIGNALING = "IL2\u2013STAT5 signalling",
  HALLMARK_DNA_REPAIR = "DNA repair", HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon-\u03b3 response",
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon-\u03b1 response", HALLMARK_GLYCOLYSIS = "Glycolysis",
  HALLMARK_XENOBIOTIC_METABOLISM = "Xenobiotic metabolism", HALLMARK_UV_RESPONSE_UP = "UV response up", HALLMARK_UV_RESPONSE_DN = "UV response down")
G[, name := ifelse(pathway %in% names(NAMES), NAMES[pathway], tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", pathway)))))]
G[, dir := fifelse(NES > 0, "up in Yq-null", "down in Yq-null")]
setorder(G, NES); G[, name := factor(name, name)]
G[, plab := sapply(padj, function(p) if (p < 1e-3) sprintf("%.0e", p) else sprintf("%.3f", p))]
save_source(G[, .(pathway, name, NES, pval, padj, size, dir)], "main_gse336056_yq_de", "b")
PAL_D <- c(`up in Yq-null` = NPG[["vermilion"]], `down in Yq-null` = NPG[["blue"]])
pb <- ggplot(G, aes(x = NES, y = name, fill = dir)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = plab, x = ifelse(NES > 0, -0.08, 0.08), hjust = ifelse(NES > 0, 1, 0)), size = 6 / .pt, colour = "#444444") +
  geom_vline(xintercept = 0, linewidth = LINE_W, colour = "grey40") +
  scale_fill_manual(values = PAL_D, name = NULL, breaks = names(PAL_D)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = "NES (Hallmark; bar labels = FDR)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme(axis.text.y = element_text(size = SMALL_PT), legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())

p <- (pa | pb) + plot_layout(widths = c(0.5, 0.5)) + tags()
save_fig(p, "main_gse336056_yq_de", h_mm = 105)
message("done")
