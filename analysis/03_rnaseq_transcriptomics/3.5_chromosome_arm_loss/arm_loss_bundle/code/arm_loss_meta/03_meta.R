#!/usr/bin/env Rscript
# ============================================================================
# Arm-loss meta-analysis: is there a conserved transcriptional program across
# independent engineered chromosome-arm / whole-chromosome losses, and how much
# of it is just p53 activation?
#
# Inputs: results/gsea/*_GSEA.tsv.gz, results/de/*_DE.tsv.gz, results/contrast_qc.tsv
#
# Framing follows the project's signed-meta convention: rank by DIRECTION and
# consistency across datasets, not by per-dataset FDR. A pathway that is FDR-
# significant in one contrast and absent in twelve is weaker evidence than one
# that leans the same way in twelve at moderate size.
#
# Three questions:
#   Q1 conserved  which pathways move the same way across independent arm losses
#                 (sign test over contrasts; direction, not per-contrast FDR)
#   Q2 p53        how much of each arm-loss profile is reproduced by the pure
#                 p53-induction contrasts (+dox vs -dox at FIXED karyotype).
#                 A high correlation means the program is p53, not dosage.
#   Q3 scaling    does effect size grow with the number of genes lost
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "/global/home/hpc5173/projects/LOY/RNA/arm_loss_meta"
OUT  <- file.path(ROOT, "results")
zf <- function(f) fread(cmd = paste("zcat", shQuote(f)))

QC  <- fread(file.path(OUT, "contrast_qc.tsv"))
ann <- zf(file.path(OUT, "gene_annotation.tsv.gz"))
G   <- rbindlist(lapply(list.files(file.path(OUT, "gsea"), "\\.tsv\\.gz$", full.names = TRUE), zf),
                 fill = TRUE)
G[, leadingEdge := NULL]
ARM <- G[type == "arm_loss"]; P53 <- G[type == "p53_induction"]
arms <- sort(unique(ARM$contrast)); n_arm <- length(arms)
cat(sprintf("arm-loss contrasts: %d (%s)\n", n_arm, paste(arms, collapse = ", ")))
cat(sprintf("p53-induction contrasts: %d (%s)\n\n", uniqueN(P53$contrast), paste(sort(unique(P53$contrast)), collapse = ", ")))

# ---- Q1: conserved direction across independent arm losses -----------------------
M <- ARM[!is.na(NES), .(
  n        = .N,
  n_up     = sum(NES > 0),
  n_dn     = sum(NES < 0),
  n_sig    = sum(padj < 0.05),
  n_sig_up = sum(padj < 0.05 & NES > 0),
  n_sig_dn = sum(padj < 0.05 & NES < 0),
  mean_NES = mean(NES),
  min_NES  = min(NES), max_NES = max(NES)
), by = .(collection, pathway)]
M <- M[n >= ceiling(0.75 * n_arm)]                      # present in most contrasts
M[, concord := pmax(n_up, n_dn) / n]
M[, sign_p := mapply(function(k, nn) binom.test(k, nn, 0.5)$p.value, pmax(n_up, n_dn), n)]
M[, sign_fdr := p.adjust(sign_p, "BH"), by = collection]
M[, dir := fifelse(mean_NES > 0, "up", "down")]
M[, abs_NES := abs(mean_NES)]
setorder(M, sign_p, -abs_NES)
fwrite(M, file.path(OUT, "meta_pathway_conservation.tsv"), sep = "\t")

show <- function(d, k = 25) print(d[seq_len(min(k, .N)),
  .(pathway, n, up = n_up, dn = n_dn, sig = n_sig, NES = round(mean_NES, 2),
    rng = sprintf("%+.1f..%+.1f", min_NES, max_NES), fdr = signif(sign_fdr, 2))], nrows = 40)
cat("=== Q1 Hallmark: most consistently DOWN across arm losses ===\n")
show(M[collection == "H" & dir == "down"][order(-concord, mean_NES)])
cat("\n=== Q1 Hallmark: most consistently UP across arm losses ===\n")
show(M[collection == "H" & dir == "up"][order(-concord, -mean_NES)])
cat("\n=== Q1 Hallmark: LEAST consistent (direction splits across contrasts) ===\n")
show(M[collection == "H"][order(concord)], 12)
cat("\n=== Q1 custom ribosomal-protein set, per contrast ===\n")
print(ARM[collection == "CUSTOM", .(contrast, NES = round(NES, 2), padj = signif(padj, 2), size)][order(NES)])
cat("\n=== Q1 interferon + inflammatory, per contrast ===\n")
IFN <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
         "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE")
print(dcast(ARM[pathway %in% IFN], contrast ~ pathway, value.var = "NES")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], nrows = 30)
cat("\n=== Q1 p53 / proliferation, per contrast ===\n")
PRO <- c("HALLMARK_P53_PATHWAY", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
         "HALLMARK_MYC_TARGETS_V1")
print(dcast(ARM[pathway %in% PRO], contrast ~ pathway, value.var = "NES")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], nrows = 30)

# ---- Q2: how much of each arm-loss profile is pure p53? --------------------------
if (nrow(P53)) {
  p53w <- dcast(P53[collection %in% c("H", "GOBP")], pathway ~ contrast, value.var = "NES")
  armw <- dcast(ARM[collection %in% c("H", "GOBP")], pathway ~ contrast, value.var = "NES")
  W <- merge(p53w, armw, by = "pathway")
  pc <- setdiff(names(p53w), "pathway"); ac <- setdiff(names(armw), "pathway")
  cr <- CJ(arm = ac, p53 = pc)
  cr[, rho := mapply(function(a, p) {
    ok <- complete.cases(W[[a]], W[[p]]); cor(W[[a]][ok], W[[p]][ok], method = "spearman")
  }, arm, p53)]
  cat("\n=== Q2 Spearman correlation of NES profiles: arm loss vs pure p53 induction ===\n")
  print(dcast(cr, arm ~ p53, value.var = "rho")[
    , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], nrows = 30)
  # arm-vs-arm agreement, for scale: how alike are two independent arm losses?
  aa <- CJ(a1 = ac, a2 = ac)[a1 < a2]
  aa[, rho := mapply(function(x, y) {
    ok <- complete.cases(W[[x]], W[[y]]); cor(W[[x]][ok], W[[y]][ok], method = "spearman")
  }, a1, a2)]
  cat(sprintf("\narm-vs-arm NES agreement: median rho = %.2f (range %.2f..%.2f, %d pairs)\n",
              median(aa$rho), min(aa$rho), max(aa$rho), nrow(aa)))
  cat(sprintf("arm-vs-p53  NES agreement: median rho = %.2f (range %.2f..%.2f)\n",
              median(cr$rho), min(cr$rho), max(cr$rho)))
  fwrite(cr, file.path(OUT, "meta_arm_vs_p53_rho.tsv"), sep = "\t")
  fwrite(aa, file.path(OUT, "meta_arm_vs_arm_rho.tsv"), sep = "\t")
}

# ---- Q3: does effect size scale with how much was lost? --------------------------
CIS <- fread(text = '
contrast|cis
A2058_1q|arm:1q
A2058_7p|arm:7p
A2058_8q|arm:8q
AGS_1q_vsTK|arm:1q
HCT116_8q|arm:8q
MCF10A_1q|arm:1q
MCF10Ahr_1q|arm:1q
A2780_1q_vsWT|arm:1q
A2780_1q_vsTK|arm:1q
RPTEC_3p_early|arm:3p
RPTEC_3p_late|arm:3p
RPE1_mono1018|chrom:chr10,chr18
RPE1_mono13|chrom:chr13
RPE1_mono19p|arm:19p
ARPE19_chrY|chrom:chrY
', sep = "|", header = TRUE)
# chrY loss in a male is HEMIZYGOUS: 1 copy -> 0, i.e. complete loss of the
# segment. Every other contrast here is 2 -> 1 (or 3 -> 2 in the aneuploid cancer
# lines), a halving. Their cis log2FC sits at -0.35..-0.84; ARPE19_chrY sits at
# -6.26 because the denominator goes to zero. These are different kinds of dosage
# change, so chrY is kept for the pathway questions (Q1/Q2) but EXCLUDED from the
# size-scaling fit (Q3), where the x-axis assumes a comparable per-gene dosage step.
CIS[, hemizygous := contrast == "ARPE19_chrY"]
CIS[, n_lost := vapply(cis, function(s) {
  p <- strsplit(s, ":")[[1]]; v <- strsplit(p[2], ",")[[1]]
  if (p[1] == "arm") ann[arm %in% v & biotype == "protein_coding", .N]
  else ann[chrom %in% v & biotype == "protein_coding", .N]
}, integer(1))]
S <- merge(CIS, QC[, .(contrast, n_padj05, median_lfc_cis, n_case, n_ctrl)], by = "contrast")
S <- merge(S, ARM[pathway == "HALLMARK_P53_PATHWAY", .(contrast, p53_NES = NES)], by = "contrast", all.x = TRUE)
S <- merge(S, ARM[pathway == "HALLMARK_E2F_TARGETS", .(contrast, e2f_NES = NES)], by = "contrast", all.x = TRUE)
setorder(S, n_lost)
cat("\n=== Q3 effect size vs protein-coding genes lost ===\n")
print(S[, .(contrast, n_lost, cis_lfc = median_lfc_cis, n_DE = n_padj05,
            p53_NES = round(p53_NES, 2), e2f_NES = round(e2f_NES, 2))], nrows = 30)
cat("  (ARPE19_chrY excluded from the fits below: hemizygous 1->0, not 2->1)\n")
for (v in c("n_padj05", "p53_NES", "e2f_NES")) {
  x <- S[!is.na(get(v)) & hemizygous == FALSE]
  if (nrow(x) >= 4) cat(sprintf("  %-10s vs genes lost: Spearman rho = %+.2f (n=%d)\n",
                                v, cor(x$n_lost, x[[v]], method = "spearman"), nrow(x)))
}
fwrite(S, file.path(OUT, "meta_scaling.tsv"), sep = "\t")
cat("\ndone\n")

# ---- Q3b: scaling at the SEGMENT level ------------------------------------------
# The contrast-level correlation is pseudo-replicated: 5 of 14 contrasts are the
# same 1q segment, 2 are 8q, 2 are 3p. The unit for a dose-response in "amount
# lost" is the distinct segment, so collapse to it before correlating.
WATCH <- c(p53 = "HALLMARK_P53_PATHWAY", e2f = "HALLMARK_E2F_TARGETS",
           myc = "HALLMARK_MYC_TARGETS_V1", ifna = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
           rp = "RIBOSOMAL_PROTEIN")
A <- merge(ARM[pathway %in% WATCH, .(contrast, pathway, NES)], CIS[, .(contrast, cis, n_lost)],
           by = "contrast")
A[, seg := sub("^(arm|chrom):", "", cis)]
A <- A[seg != "chrY"]                      # hemizygous; see note above
SEG <- dcast(A[, .(NES = mean(NES), n_bg = .N), by = .(seg, n_lost, pathway)],
             seg + n_lost ~ pathway, value.var = "NES")
setorder(SEG, n_lost)
cat("\n=== Q3b segment-level means (distinct segments, n =", nrow(SEG), ") ===\n")
print(SEG[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], nrows = 20)
for (nm in names(WATCH)) {
  p <- WATCH[[nm]]; if (!p %in% names(SEG)) next
  x <- SEG[!is.na(get(p))]
  ct <- suppressWarnings(cor.test(x$n_lost, x[[p]], method = "spearman"))
  cat(sprintf("  %-5s vs genes lost (segment level): rho = %+.2f, P = %.3f (n=%d)\n",
              nm, ct$estimate, ct$p.value, nrow(x)))
}
cat("\n=== within-segment consistency across cell backgrounds ===\n")
CON <- A[, .(n_bg = .N, n_up = sum(NES > 0), mean = round(mean(NES), 2),
             rng = sprintf("%+.1f..%+.1f", min(NES), max(NES))), by = .(pathway, seg, n_lost)]
print(CON[n_bg >= 2][order(pathway, n_lost)], nrows = 40)
fwrite(SEG, file.path(OUT, "meta_scaling_segment.tsv"), sep = "\t")
fwrite(CON, file.path(OUT, "meta_within_segment.tsv"), sep = "\t")
