# Stage 03: call Y status for every male model from the Stage 02 evidence.
# Thresholds come from parameters.json and from the female (no-Y) reference in
# the evidence table. No outcome table (CRISPR, drugs) is read here.
# Output: y_status_calls.csv (one row per model), call_thresholds.json, call_summary.csv.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))

core <- params$y_genes_core
rna_par <- params$rna_call
dna_par <- params$dna_call
rna_cols <- paste0("rna_", core)
yp_genes <- c("ZFY", "RPS4Y1")
yq_genes <- setdiff(core, yp_genes)

evidence <- fread(file.path(derived_dir, "y_evidence.csv"))
assert_unique(evidence$ModelID, "y_evidence")

# --- female reference --------------------------------------------------------
females <- evidence[Sex == "Female" & has_rna == TRUE]
males <- evidence[Sex == "Male" & has_rna == TRUE]
female_gene_max <- vapply(rna_cols, function(col) max(females[[col]]), numeric(1))
female_mean_max <- max(females$rna_y_mean)
male_counts <- c(
  le_depleted = sum(males$rna_y_mean <= rna_par$depleted_mean_max),
  between = sum(males$rna_y_mean > rna_par$depleted_mean_max & males$rna_y_mean < rna_par$expressed_mean_min),
  ge_expressed = sum(males$rna_y_mean >= rna_par$expressed_mean_min)
)
expected_counts <- c(rna_par$expected_male_le_depleted, rna_par$expected_male_between, rna_par$expected_male_ge_expressed)
log_msg("female reference: n=%d, mean max %.3f; male counts %s", nrow(females), female_mean_max, paste(male_counts, collapse = "/"))
# The thresholds were chosen on these numbers; a different input would need a new decision.
stopifnot(
  abs(female_mean_max - rna_par$expected_female_mean_max) < 0.005,
  female_mean_max <= rna_par$depleted_mean_max,
  all(abs(male_counts - expected_counts) <= rna_par$tolerance_count)
)

# --- call functions ----------------------------------------------------------
call_rna <- function(dt, mean_depleted, mean_expressed, gene_max, breadth = TRUE) {
  n_on <- rowSums(sweep(as.matrix(dt[, rna_cols, with = FALSE]), 2, gene_max, ">"))
  mean_y <- dt$rna_y_mean
  depleted <- mean_y <= mean_depleted
  expressed <- mean_y >= mean_expressed
  if (breadth) {
    depleted <- depleted & n_on <= rna_par$depleted_max_genes_on
    expressed <- expressed & n_on >= rna_par$expressed_min_genes_on
  }
  out <- fifelse(depleted, "depleted", fifelse(expressed, "expressed", "intermediate"))
  out[!dt$has_rna] <- NA_character_
  list(call = out, n_on = fifelse(dt$has_rna, as.integer(n_on), NA_integer_))
}

# WGS is used when present, otherwise WES. Gene ratios are relative to the
# sample's own mean, so one Y copy in a diploid male is about 0.5.
dna_values <- function(dt) {
  source <- fifelse(dt$has_wgs, "WGS", fifelse(dt$has_wes, "WES", NA_character_))
  values <- matrix(NA_real_, nrow(dt), length(core), dimnames = list(NULL, core))
  for (g in core) {
    values[, g] <- fifelse(source == "WGS", dt[[paste0("wgs_", g)]], dt[[paste0("wes_", g)]])
  }
  list(source = source, values = values)
}

call_dna <- function(dt, loss_max, retained_min, floor_min) {
  dv <- dna_values(dt)
  v <- dv$values
  n_available <- rowSums(!is.na(v))
  low <- v <= loss_max
  high <- v >= retained_min
  all_low <- rowSums(low, na.rm = TRUE) == n_available
  all_high <- rowSums(high, na.rm = TRUE) == n_available
  any_low <- rowSums(low, na.rm = TRUE) > 0
  any_high <- rowSums(high, na.rm = TRUE) > 0
  # For WGS, whole loss also needs the segment evidence when it exists.
  floor_ok <- dv$source != "WGS" | is.na(dt$wgs_y_floor_fraction) | dt$wgs_y_floor_fraction >= floor_min
  out <- fifelse(
    all_low & floor_ok, "loss",
    fifelse(all_high, "retained", fifelse(any_low & any_high, "partial", "mosaic"))
  )
  out[is.na(dv$source) | n_available < 5] <- NA_character_
  arm <- function(genes, test) rowSums(test[, genes, drop = FALSE], na.rm = TRUE) == rowSums(!is.na(v[, genes, drop = FALSE]))
  partial_type <- fifelse(
    arm(yq_genes, low) & arm(yp_genes, high), "Yq_loss",
    fifelse(arm(yp_genes, low) & arm(yq_genes, high), "Yp_loss", "other")
  )
  partial_type[is.na(out) | out != "partial"] <- NA_character_
  list(call = out, source = dv$source, partial_type = partial_type, ratio = apply(v, 1, median, na.rm = TRUE))
}

# RNA leads; DNA fills in where RNA is missing and vetoes a clear contradiction.
integrate_calls <- function(rna, dna, male) {
  rna_as <- fifelse(rna == "depleted", "Y_negative", fifelse(rna == "expressed", "Y_positive", NA_character_))
  dna_as <- fifelse(dna == "loss", "Y_negative", fifelse(dna == "retained", "Y_positive", NA_character_))
  rna_definite <- !is.na(rna_as)
  dna_definite <- !is.na(dna_as)
  dna_ambiguous <- !is.na(dna) & dna %in% c("partial", "mosaic")
  tier <- rep(NA_character_, length(rna))
  tier[rna_definite & dna_definite & rna_as == dna_as] <- "concordant"
  tier[rna_definite & dna_definite & rna_as != dna_as] <- "discordant"
  tier[rna_definite & dna_ambiguous] <- "rna_led_dna_ambiguous"
  tier[rna_definite & is.na(dna)] <- "rna_only"
  tier[is.na(rna) & dna_definite] <- "dna_only"
  tier[!is.na(rna) & rna == "intermediate"] <- "rna_intermediate"
  tier[is.na(rna) & dna_ambiguous] <- "dna_ambiguous_no_rna"
  tier[is.na(rna) & is.na(dna)] <- "no_data"
  call <- fifelse(
    tier %in% c("concordant", "rna_led_dna_ambiguous", "rna_only"), rna_as,
    fifelse(tier == "dna_only", dna_as, NA_character_)
  )
  tier[!male] <- "not_male"
  call[!male] <- NA_character_
  list(call = call, tier = tier)
}

make_calls <- function(dt, scale = 1, breadth = TRUE) {
  rna <- call_rna(dt, rna_par$depleted_mean_max * scale, rna_par$expressed_mean_min * scale, female_gene_max * scale, breadth)
  dna <- call_dna(dt, dna_par$loss_ratio_max * scale, dna_par$retained_ratio_min * scale, dna_par$wgs_floor_fraction_min)
  integrated <- integrate_calls(rna$call, dna$call, dt$Sex == "Male")
  list(rna = rna, dna = dna, integrated = integrated)
}

# --- primary calls -----------------------------------------------------------
primary <- make_calls(evidence)
calls <- copy(evidence)
calls[, rna_n_genes_on := primary$rna$n_on]
calls[, rna_call := primary$rna$call]
for (g in core) calls[, (paste0("rna_on_", g)) := fifelse(has_rna, get(paste0("rna_", g)) > female_gene_max[[paste0("rna_", g)]], NA)]
calls[, dna_source := primary$dna$source]
calls[, dna_y_ratio := primary$dna$ratio]
calls[, dna_call := primary$dna$call]
calls[, dna_partial_type := primary$dna$partial_type]
calls[, y_call := primary$integrated$call]
calls[, y_tier := primary$integrated$tier]

# Descriptive copy estimate; ploidy is missing for many profiles, so it never enters the call.
calls[, y_copies := dna_y_ratio * Ploidy]
calls[, control_one_copy := y_copies >= dna_par$one_copy_min & y_copies <= dna_par$one_copy_max]
calls[, control_gain := y_copies >= dna_par$gain_copies_min]
calls[, control_relative_gain := dna_y_ratio > dna_par$relative_gain_ratio_min]

# Continuous score: 0 at the female median, 1 at the median of expressed males.
female_median <- median(females$rna_y_mean)
expressed_median <- median(calls[Sex == "Male" & rna_call == "expressed", rna_y_mean])
calls[, rna_y_score := (rna_y_mean - female_median) / (expressed_median - female_median)]

# Possible sex mis-annotation: a "male" without Y that also looks two-X and XIST-high.
calls[, x_ratio := fifelse(!is.na(wgs_x_ratio), wgs_x_ratio, wes_x_ratio)]
calls[, xist_high := rna_XIST >= params$sex_flags$xist_log2tpm_min]
calls[, x_two_copy := x_ratio >= params$sex_flags$x_ratio_two_copy_min]
calls[, possible_female := Sex == "Male" & y_call %in% "Y_negative" & xist_high %in% TRUE & x_two_copy %in% TRUE]

# --- alternative exposures for sensitivity analyses --------------------------
male <- calls$Sex == "Male"
lenient <- make_calls(evidence, breadth = FALSE)
calls[, rna_call_lenient := lenient$rna$call]
calls[, y_call_lenient := lenient$integrated$call]
calls[, y_call_concordant := fifelse(y_tier == "concordant", y_call, NA_character_)]
calls[, y_call_rna := fifelse(male & rna_call == "depleted", "Y_negative", fifelse(male & rna_call == "expressed", "Y_positive", NA_character_))]
calls[, y_call_dna := fifelse(male & dna_call == "loss", "Y_negative", fifelse(male & dna_call == "retained", "Y_positive", NA_character_))]
calls[, y_call_one_copy := fifelse(y_call == "Y_negative" | (y_call == "Y_positive" & control_one_copy %in% TRUE), y_call, NA_character_)]
calls[, y_call_no_relative_gain := fifelse(y_call == "Y_positive" & control_relative_gain %in% TRUE, NA_character_, y_call)]
calls[, y_call_no_xist := fifelse(xist_high %in% TRUE, NA_character_, y_call)]

# Stability: the same call after shifting every threshold down or up by the stated fraction.
shift <- params$stability_shift
low_shift <- make_calls(evidence, scale = 1 - shift)$integrated$call
high_shift <- make_calls(evidence, scale = 1 + shift)$integrated$call
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
calls[, call_stable := same(y_call, low_shift) & same(y_call, high_shift)]
calls[, y_call_stable := fifelse(call_stable, y_call, NA_character_)]

# --- outputs -----------------------------------------------------------------
assert_unique(calls$ModelID, "y_status_calls")
write_result(calls, "y_status_calls")

summary_tab <- calls[Sex == "Male", .N, by = .(y_tier, y_call)][order(y_tier, y_call)]
write_result(summary_tab, "call_summary")
print(summary_tab)
totals <- calls[Sex == "Male", .(
  males = .N, y_negative = sum(y_call %in% "Y_negative"), y_positive = sum(y_call %in% "Y_positive"),
  set_aside = sum(is.na(y_call)), unstable = sum(!call_stable), possible_female = sum(possible_female)
)]
print(totals)
print(calls[Sex == "Male" & has_rna == TRUE, .N, by = rna_call])
print(calls[Sex == "Male" & !is.na(dna_call), .N, by = .(dna_source, dna_call)][order(dna_source, dna_call)])

jsonlite::write_json(list(
  female_reference_n = nrow(females),
  female_mean_max = female_mean_max,
  female_gene_max = as.list(setNames(female_gene_max, core)),
  female_median = female_median,
  expressed_male_median = expressed_median,
  male_rna_counts = as.list(male_counts),
  rna_call = rna_par, dna_call = dna_par, sex_flags = params$sex_flags,
  totals = as.list(totals)
), file.path(derived_dir, "call_thresholds.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("calls written: %d Y-negative, %d Y-positive, %d set aside", totals$y_negative, totals$y_positive, totals$set_aside)
