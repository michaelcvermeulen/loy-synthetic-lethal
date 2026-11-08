# Matched region sets, motif presence and motif tests for the HT1080
# differential-region heatmap figure.
#
# Targets: primary FDR < 0.05 lost (663) and gained (95) regions.
# Controls: tested regions with primary FDR >= 0.20 and no N bases, sampled without
# replacement within the SEPT_08 matching strata (5-percentage-point GC bins x
# baseline logCPM deciles; motifs/peak_covariates.csv), K controls per target.
# Lost controls are drawn first and gained controls from the remainder, so the two
# control sets are disjoint. Strata without enough controls are subsampled on the
# target side (reported in the accounting table), as in the original motif analysis.
# Motif presence: chromVAR motifmatchr matches (JASPAR2024, site P < 5e-5) for the
# seven representative HOMER motifs. Test: two-sided Cochran-Mantel-Haenszel test
# stratified by matching stratum, target vs its own controls; BH across all
# 7 motifs x 2 directions.
# Run: Rscript --vanilla <absolute path to this script>
suppressPackageStartupMessages({library(Matrix); library(data.table)})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
fig <- dirname(dirname(normalizePath(script, winslash = "/")))
data_dir <- file.path(fig, "source_data/ht1080_atac_accessibility")
run <- file.path(dirname(dirname(dirname(fig))), "ATACseq/SEPT_08_2026_REPROCESS")
K <- 3L
SEED <- 20260922L

da <- fread(file.path(data_dir, "differential_peaks.csv"))
cov <- fread(file.path(run, "motifs/peak_covariates.csv"))
acc <- fread(file.path(run, "motifs/TF_extension/chromVAR/peak_accounting.csv"))
mm <- readRDS(file.path(run, "motifs/TF_extension/chromVAR/motif_matches.rds"))
chromvar_peaks <- acc[included == TRUE, peak_id]
stopifnot(nrow(mm) == length(chromvar_peaks))
reps <- data.table(family = c("TFAP2", "EBF", "ZIC", "TEAD", "RUNX", "GATA", "SOX6"),
                   motif_id = c("MA0524.3", "MA0154.5", "MA0696.1", "MA0090.4", "MA0002.3", "MA1104.3", "MA0515.1"))

d <- merge(da[, .(peak_id, chrom, start, end, logFC, FDR, direction, annotation, nearest_symbol)],
           cov[, .(peak_id, GC, N_fraction, logCPM, stratum)], by = "peak_id")
stopifnot(nrow(d) == nrow(da))
d[, in_chromvar := peak_id %in% chromvar_peaks]
pool <- d[FDR >= 0.20 & N_fraction == 0 & in_chromvar == TRUE]

set.seed(SEED)
used <- character(0)
sets <- list(); accounting <- list()
for (dir in c("lost", "gained")) {
  tg <- d[direction == dir & N_fraction == 0 & in_chromvar == TRUE]
  for (s in sort(unique(tg$stratum))) {
    t_s <- tg[stratum == s][order(peak_id)]
    c_s <- pool[stratum == s & !(peak_id %in% used)][order(peak_id)]
    n_t <- min(nrow(t_s), floor(nrow(c_s) / K))
    keep_t <- if (n_t < nrow(t_s)) t_s[sort(sample.int(nrow(t_s), n_t))] else t_s
    pick_c <- if (n_t > 0) c_s[sort(sample.int(nrow(c_s), n_t * K))] else c_s[0]
    used <- c(used, pick_c$peak_id)
    sets[[length(sets) + 1]] <- rbind(keep_t[, group := dir], pick_c[, group := paste0(dir, "_ref")])
    accounting[[length(accounting) + 1]] <- data.table(direction = dir, stratum = s, targets = nrow(t_s),
                                                       matched_targets = n_t, controls = nrow(pick_c))
  }
}
regions <- rbindlist(sets)
accounting <- rbindlist(accounting)
regions[, centre := (start + end) %/% 2L]
row <- match(regions$peak_id, chromvar_peaks)
stopifnot(!anyNA(row))
for (i in seq_len(nrow(reps))) regions[[reps$family[i]]] <- as.logical(mm[row, reps$motif_id[i]])
setcolorder(regions, c("group", "peak_id", "chrom", "start", "end", "centre", "logFC", "logCPM", "GC", "FDR",
                       "stratum", "annotation", "nearest_symbol"))
regions[, c("N_fraction", "in_chromvar", "direction") := NULL]

balance <- regions[, .(n = .N, mean_logCPM = mean(logCPM), mean_GC = mean(GC)), by = group]
smd <- function(a, b) (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
balance[, SMD_vs_targets := NA_real_]
for (dir in c("lost", "gained")) {
  t_ <- regions[group == dir]; c_ <- regions[group == paste0(dir, "_ref")]
  balance[group == paste0(dir, "_ref"), `:=`(SMD_logCPM = smd(c_$logCPM, t_$logCPM), SMD_GC = smd(c_$GC, t_$GC))]
}
balance[, SMD_vs_targets := NULL]

tests <- rbindlist(lapply(c("lost", "gained"), function(dir) rbindlist(lapply(reps$family, function(fam) {
  x <- regions[group %in% c(dir, paste0(dir, "_ref"))]
  x[, is_target := group == dir]
  # Strata with a single row carry no information for CMH; they are dropped by construction (K >= 1).
  tab <- table(factor(x$is_target, c(TRUE, FALSE)), factor(x[[fam]], c(TRUE, FALSE)), x$stratum)
  tab <- tab[, , apply(tab, 3, sum) > 1, drop = FALSE]
  mh <- mantelhaen.test(tab, alternative = "two.sided", correct = TRUE)
  data.table(direction = dir, family = fam, motif_id = reps[family == fam, motif_id],
             target_n = sum(x$is_target), control_n = sum(!x$is_target),
             target_pct = 100 * mean(x[is_target == TRUE][[fam]]),
             control_pct = 100 * mean(x[is_target == FALSE][[fam]]),
             MH_odds_ratio = unname(mh$estimate), CI_low = mh$conf.int[1], CI_high = mh$conf.int[2],
             CMH_P = mh$p.value, strata = dim(tab)[3])
}))))
tests[, BH_FDR := p.adjust(CMH_P, method = "BH")]

fwrite(regions, file.path(data_dir, "region_heatmap_regions.csv"))
fwrite(accounting, file.path(data_dir, "region_heatmap_matching_accounting.csv"))
fwrite(balance, file.path(data_dir, "region_heatmap_matching_balance.csv"))
fwrite(tests, file.path(data_dir, "region_heatmap_motif_CMH.csv"))
fwrite(reps, file.path(data_dir, "region_heatmap_motif_representatives.csv"))
writeLines(c(sprintf("K = %d controls per target; seed %d; controls primary FDR >= 0.20, no N, in chromVAR universe.", K, SEED),
             "Strata: SEPT_08 5-pt GC bins x baseline logCPM deciles (motifs/peak_covariates.csv).",
             "Test: two-sided CMH (continuity corrected) stratified by matching stratum; BH across 14 tests."),
           file.path(data_dir, "region_heatmap_matching_METHODS.txt"))
print(accounting[, .(targets = sum(targets), matched = sum(matched_targets), controls = sum(controls)), by = direction])
print(balance)
print(tests[, .(direction, family, target_pct = round(target_pct, 1), control_pct = round(control_pct, 1),
                OR = round(MH_odds_ratio, 2), BH_FDR = signif(BH_FDR, 3))])
