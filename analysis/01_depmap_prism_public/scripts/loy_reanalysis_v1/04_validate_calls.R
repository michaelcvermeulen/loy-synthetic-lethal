# Stage 04: internal validation of the Y-status calls. Descriptive only: nothing
# here feeds back into a threshold. The comparison with earlier call tables is
# in 20_compare_prior.R, the only stage that reads earlier outputs.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))

core <- params$y_genes_core
dna_par <- params$dna_call
calls <- read_result("y_status_calls")
males <- calls[Sex == "Male"]

# Gene-only DNA class from one assay, so WGS and WES can be compared on the same lines.
gene_call <- function(dt, prefix) {
  v <- as.matrix(dt[, paste0(prefix, core), with = FALSE])
  n_available <- rowSums(!is.na(v))
  low <- rowSums(v <= dna_par$loss_ratio_max, na.rm = TRUE)
  high <- rowSums(v >= dna_par$retained_ratio_min, na.rm = TRUE)
  out <- fifelse(low == n_available, "loss", fifelse(high == n_available, "retained", fifelse(low > 0 & high > 0, "partial", "mosaic")))
  out[n_available < 5] <- NA_character_
  out
}

# --- RNA against DNA ---------------------------------------------------------
rna_dna <- males[!is.na(rna_call) & !is.na(dna_call), .N, by = .(rna_call, dna_call, dna_source)][order(dna_source, rna_call, dna_call)]
write_result(rna_dna, "validation_rna_vs_dna")
definite <- males[rna_call %in% c("depleted", "expressed") & dna_call %in% c("loss", "retained")]
agreement <- definite[, .(
  n = .N,
  both_negative = sum(rna_call == "depleted" & dna_call == "loss"),
  both_positive = sum(rna_call == "expressed" & dna_call == "retained"),
  rna_depleted_dna_retained = sum(rna_call == "depleted" & dna_call == "retained"),
  rna_expressed_dna_loss = sum(rna_call == "expressed" & dna_call == "loss")
), by = dna_source]
agreement[, agreement := (both_negative + both_positive) / n]
write_result(agreement, "validation_rna_dna_agreement")
print(agreement)

# --- WGS against WES ---------------------------------------------------------
both_assays <- males[has_wgs == TRUE & has_wes == TRUE]
both_assays[, wgs_class := gene_call(both_assays, "wgs_")]
both_assays[, wes_class := gene_call(both_assays, "wes_")]
wgs_wes <- both_assays[, .N, by = .(wgs_class, wes_class)][order(wgs_class, wes_class)]
write_result(wgs_wes, "validation_wgs_vs_wes")
print(wgs_wes)

# --- females as the no-Y reference -------------------------------------------
females <- calls[Sex == "Female"]
females[, dna_class := fifelse(has_wgs, gene_call(females, "wgs_"), gene_call(females, "wes_"))]
female_check <- females[, .(
  n = .N, with_rna = sum(has_rna), rna_mean_le_1 = sum(rna_y_mean <= params$rna_call$depleted_mean_max, na.rm = TRUE),
  with_dna = sum(!is.na(dna_class)), dna_loss = sum(dna_class %in% "loss"), dna_not_loss = sum(!is.na(dna_class) & dna_class != "loss")
)]
write_result(female_check, "validation_female_reference")
write_result(females[!is.na(dna_class) & dna_class != "loss", .(ModelID, StrippedCellLineName, OncotreeLineage, dna_class, dna_y_ratio, rna_y_mean)], "validation_female_with_y_signal")
print(female_check)

# --- unknown sex, for the record ---------------------------------------------
unknown <- calls[!Sex %in% c("Male", "Female") & has_rna == TRUE, .(
  n = .N, y_expressing = sum(rna_y_mean >= params$rna_call$expressed_mean_min), y_absent = sum(rna_y_mean <= params$rna_call$depleted_mean_max)
)]
write_result(unknown, "validation_unknown_sex")

# --- stability and flags -----------------------------------------------------
stability <- males[, .(n = .N, unstable = sum(!call_stable)), by = .(y_tier)][order(y_tier)]
write_result(stability, "validation_call_stability")
flags <- males[!is.na(y_call), .(
  n = .N, xist_high = sum(xist_high %in% TRUE), x_two_copy = sum(x_two_copy %in% TRUE), possible_female = sum(possible_female)
), by = y_call]
write_result(flags, "validation_sex_flags")
print(flags)
write_result(males[possible_female == TRUE, .(ModelID, StrippedCellLineName, OncotreeLineage, y_tier, rna_y_mean, rna_XIST, x_ratio, dna_y_ratio)], "validation_possible_female_lines")

# --- partial deletions -------------------------------------------------------
partial <- males[dna_call == "partial", .N, by = .(dna_partial_type, rna_call)][order(dna_partial_type, rna_call)]
write_result(partial, "validation_partial_deletions")

by_lineage <- males[!is.na(y_call), .(y_negative = sum(y_call == "Y_negative"), y_positive = sum(y_call == "Y_positive")), by = OncotreeLineage][order(-(y_negative + y_positive))]
write_result(by_lineage, "validation_calls_by_lineage")

wgs_depleted <- males[rna_call == "depleted" & dna_source == "WGS" & dna_call %in% c("loss", "retained")]
jsonlite::write_json(list(
  rna_depleted_wgs_definite = nrow(wgs_depleted),
  rna_depleted_wgs_loss = sum(wgs_depleted$dna_call == "loss"),
  females_with_rna = female_check$with_rna, females_rna_absent = female_check$rna_mean_le_1,
  females_with_dna = female_check$with_dna, females_dna_loss = female_check$dna_loss,
  males_unstable = sum(!males$call_stable), possible_female = sum(males$possible_female)
), file.path(derived_dir, "call_validation.json"), auto_unbox = TRUE, pretty = TRUE)
log_msg("RNA-depleted males with a definite WGS call: %d, of which DNA loss: %d", nrow(wgs_depleted), sum(wgs_depleted$dna_call == "loss"))
