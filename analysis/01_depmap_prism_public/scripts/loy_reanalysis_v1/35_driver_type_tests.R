# Stage 35: driver status inside cancer types, the tests.
# For every split frozen by stage 34 (a driver gene in a cancer type with five altered
# lines per Y group), one model on the lines of that cancer type:
#   rank-normalised score ~ library + status + Y among altered + Y among unaltered.
# It gives the Y effect among altered lines and, where all four cells hold five lines,
# the difference between the two Y effects (specificity). Ranks are the screen's, over
# all lines; the test is limma's moderated t, as in 12_context_genes.R. Raw-scale fits
# give the effect sizes. Each test is also given by plain least squares on the same ranks,
# where a gene keeps its own variance (limma's test borrows nearly all of it from the other
# genes at these sizes); those p-values describe and make no call. Controls that do not
# involve Y run first and stop the stage if they fail. Y status shuffled within cancer type and status calibrates the procedure.
# The Sanger organoids are fitted the same way without weights, and splits present in
# both sources are combined as in 28_combine_sources.R.
# LOY_EXTENSION_TRIAL=1 runs a few hundred genes and few shuffles on Y status shuffled
# within each split's cells, so a trial never shows a real result (use with
# LOY_EXTENSION_OUTPUT).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("organoid", dir = ext_dir)
assert_extension_freeze("driver_type")

dp <- ext_params$driver_type
sp <- ext_params$shared
op <- ext_params$organoid
fdr_cut <- sp$fdr
gate <- sp$min_cell
covariates <- dp$covariates
trial <- Sys.getenv("LOY_EXTENSION_TRIAL") == "1"
null_runs <- if (trial) dp$trial_null_runs else dp$null_runs
# A trial fits its "observed" run on labels shuffled like a null run, under a seed of its own.
observed_shuffle <- if (trial) 900000 else 0
lambda_of <- function(p) median(qchisq(p[is.finite(p)], 1, lower.tail = FALSE)) / qchisq(0.5, 1)
rank_rows <- function(mat_t) {
  out <- t(apply(mat_t, 1, rank_normalise))
  dimnames(out) <- dimnames(mat_t)
  out
}
unit_of <- function(design, coef) sqrt(solve(crossprod(design))[coef, coef])

splits <- read_ext("driver_type_splits")[fitted == TRUE]
split_lines <- read_ext("driver_type_split_lines")
status <- read_ext("driver_type_status")
families <- read_ext("driver_type_genes")
pairs <- read_ext("driver_type_pairs")
organoid_splits <- read_ext("driver_type_organoid_splits")
organoid_lines <- read_ext("driver_type_organoid_split_lines")

# --- cell lines: the screen's rows ------------------------------------------------------------
cohort <- read_result("cohort")
effect <- load_matrix("crispr_effect")
frame <- analysis_frame(cohort, effect)
frame[, MSI := as.numeric(MSIScore >= params$context$msi_score_min)]
pan_hits <- read_result("gene_pan_hits")
named <- pan_hits$feature
cell_genes <- families[source == "Cell lines"]
selective <- cell_genes[selective == TRUE, feature]
universe <- cell_genes$feature
if (trial) {
  universe <- union(head(selective, sp$trial_genes), cell_genes[focused == TRUE, feature])
  selective <- intersect(selective, universe)
  log_msg("trial run: %d genes, %d shuffles, Y status shuffled within each split's cells", length(universe), null_runs)
}
stopifnot(all(named %in% universe))
mat_raw <- t(effect[frame$ModelID, universe])
mat_rank <- rank_rows(mat_raw)
log_msg("driver-by-type tests: %d lines, %d genes (%d selective), %d splits (%d with a specificity test)", nrow(frame), length(universe), length(selective), nrow(splits), splits[, sum(label == split_labels[1])])

# The lines of one split as frozen, with Y status shuffled within status for shuffle > 0.
split_rows <- function(id, shuffle = 0, lines_table = split_lines, source_frame = frame, ids = splits$split_id) {
  lines <- lines_table[split_id == id]
  rows <- source_frame[match(lines$ModelID, ModelID)]
  stopifnot(nrow(lines) > 0, !anyNA(rows$ModelID), identical(as.numeric(rows$exposed), as.numeric(lines$exposed)))
  if (shuffle > 0) {
    set.seed(dp$seed + 1000 * shuffle + match(id, ids))
    rows[, exposed := permute_within(exposed, lines$ctx)]
  }
  split_frame(rows, rep(TRUE, nrow(rows)), lines$ctx)
}

# --- controls that do not involve Y status ----------------------------------------------------
controls <- as.data.table(dp$controls[c("gene", "status_gene", "level", "type", "direction")])
control_fit <- function(i) {
  rows <- frame[get(dp$levels[[controls$level[i]]]) %in% controls$type[i]]
  value <- status[[controls$status_gene[i]]][match(rows$ModelID, status$ModelID)]
  feature <- colnames(effect)[gene_symbol(colnames(effect)) == controls$gene[i]]
  stopifnot(length(feature) == 1, nrow(rows) > 0)
  nuisance <- full_rank(cbind(intercept = 1, as.matrix(rows[, covariates, with = FALSE])))
  fit <- fit_ols(effect[rows$ModelID, feature], value, nuisance)
  data.table(
    control = sprintf("%s by %s status inside %s", controls$gene[i], controls$status_gene[i], controls$type[i]),
    altered = sum(value %in% 1), unaltered = sum(value %in% 0), estimate = fit$estimate, se = fit$se, p = fit$p, expected_direction = controls$direction[i]
  )
}
control_table <- rbindlist(lapply(seq_len(nrow(controls)), control_fit))
control_table[, pass := sign(estimate) == expected_direction & p < dp$controls$p_max]
write_ext(control_table, "driver_type_controls")
print(control_table)
if (!all(control_table$pass)) stop("a control failed: ", paste(control_table[pass == FALSE, control], collapse = "; "))

# --- cell lines: the fits ------------------------------------------------------------------------
# limma warns about partially missing rows in every fit, as it did in the screen.
scan_split <- function(id, shuffle = 0, with_raw = TRUE) {
  sframe <- split_rows(id, shuffle)
  ranked <- suppressWarnings(fit_split(sframe, mat_rank, covariates, gate))
  out <- ranked[, .(
    split_id = id, feature, n, n_neg_in, n_pos_in, n_neg_out, n_pos_out,
    inside_rank = inside_estimate, inside_p, inside_df, outside_rank = outside_estimate, outside_p, interaction_rank = interaction_estimate, interaction_p
  )]
  out[, c("inside_p_ols", "interaction_p_ols") := fit_split_ols(sframe, mat_rank, covariates, gate)[, .(inside_p_ols, interaction_p_ols)]]
  if (with_raw) {
    raw <- suppressWarnings(fit_split(sframe, mat_raw, covariates, gate, trend = TRUE))
    out[, c("inside", "inside_se", "inside_p_raw", "outside", "outside_se", "interaction", "interaction_se", "interaction_p_raw") :=
      raw[, .(inside_estimate, inside_se, inside_p, outside_estimate, outside_se, interaction_estimate, interaction_se, interaction_p)]]
  }
  out
}
call_hits <- function(dt) {
  dt[, inside_fdr := bh(inside_p)]
  dt[, interaction_fdr := bh(interaction_p)]
  dt[, inside_fdr_split := bh(inside_p), by = split_id]
  # A call needs the Y effect among altered lines and its difference from the unaltered lines.
  dt[, specific := !is.na(interaction_fdr) & interaction_fdr < fdr_cut & inside_fdr < fdr_cut]
  dt[, inside_only := inside_fdr < fdr_cut & !specific]
  # Least squares with each gene's own variance: a description beside the rule, never a call.
  dt[, inside_fdr_ols := bh(inside_p_ols)]
  dt[]
}
# The four paired genes inside every split: a family named in advance.
call_named <- function(dt, features = named) {
  out <- dt[feature %in% features]
  out[, named_fdr := bh(inside_p)]
  out[, shown := named_fdr < fdr_cut & inside_rank < 0]
  out[]
}
tail_row <- function(dt, run, keep = selective, named_features = named) {
  screen <- call_hits(dt[feature %in% keep])
  shown <- call_named(dt, named_features)
  rate <- function(p, cut) if (any(is.finite(p))) mean(p[is.finite(p)] < cut) else NA_real_
  data.table(
    run = run, inside_tests = sum(is.finite(screen$inside_p)), interaction_tests = sum(is.finite(screen$interaction_p)),
    inside_below_0.05 = rate(screen$inside_p, 0.05), inside_below_0.001 = rate(screen$inside_p, 0.001), inside_below_1e5 = rate(screen$inside_p, 1e-5),
    interaction_below_0.05 = rate(screen$interaction_p, 0.05), interaction_below_0.001 = rate(screen$interaction_p, 0.001), interaction_below_1e5 = rate(screen$interaction_p, 1e-5),
    inside_lambda = lambda_of(screen$inside_p), interaction_lambda = lambda_of(screen$interaction_p),
    inside_min_p = min(screen$inside_p, na.rm = TRUE), interaction_min_p = min(screen$interaction_p, na.rm = TRUE),
    inside_fdr_hits = sum(screen$inside_fdr < fdr_cut, na.rm = TRUE), split_alone_hits = sum(screen$inside_fdr_split < fdr_cut, na.rm = TRUE),
    specific = sum(screen$specific), inside_only = sum(screen$inside_only, na.rm = TRUE),
    named_tests = nrow(shown), named_shown = sum(shown$shown, na.rm = TRUE), named_below_0.05 = rate(shown$inside_p, 0.05),
    ols_below_0.05 = rate(screen$inside_p_ols, 0.05), ols_below_1e5 = rate(screen$inside_p_ols, 1e-5), ols_lambda = lambda_of(screen$inside_p_ols),
    ols_fdr_hits = sum(screen$inside_fdr_ols < fdr_cut, na.rm = TRUE), ols_min_p = min(screen$inside_p_ols, na.rm = TRUE)
  )
}

results <- rbindlist(parallel::mclapply(splits$split_id, scan_split, shuffle = observed_shuffle, with_raw = TRUE, mc.cores = n_cores))
stopifnot(nrow(results) == nrow(splits) * length(universe), results[, all(is.finite(inside_p))])
# No specificity result for a split without four cells of five lines.
with_four <- splits[label == split_labels[1], split_id]
stopifnot(results[!split_id %in% with_four, all(is.na(interaction_p))], results[split_id %in% with_four, all(is.finite(interaction_p))])

# The same alteration in the other cancer types of that level, as a description beside each split.
elsewhere_fit <- function(id) {
  row <- splits[split_id == id]
  altered <- status[[row$gene]][match(frame$ModelID, status$ModelID)] %in% 1
  other <- frame[altered & !(get(dp$levels[[row$level]]) %in% row$type)]
  if (observed_shuffle > 0) {
    set.seed(dp$seed + observed_shuffle + match(id, splits$split_id))
    other[, exposed := permute_within(exposed, OncotreeLineage)]
  }
  if (min(sum(other$exposed == 1), sum(other$exposed == 0)) < gate) return(NULL)
  design <- build_design(other, tested = "exposed", covariates = covariates)
  if (!design$estimable || nrow(design$design) - ncol(design$design) < sp$min_residual_df) return(NULL)
  lines <- other$ModelID[design$rows]
  ranked <- suppressWarnings(limma_scan(mat_rank[, lines, drop = FALSE], design$design, "exposed", trend = FALSE))
  raw <- suppressWarnings(limma_scan(mat_raw[, lines, drop = FALSE], design$design, "exposed", trend = TRUE))
  data.table(split_id = id, feature = ranked$feature, elsewhere_rank = ranked$estimate, elsewhere_p = ranked$p, elsewhere = raw$estimate, elsewhere_lines = length(lines))
}
elsewhere <- rbindlist(parallel::mclapply(splits$split_id, elsewhere_fit, mc.cores = n_cores))
if (nrow(elsewhere) == 0) elsewhere <- data.table(split_id = character(), feature = character(), elsewhere_rank = numeric(), elsewhere_p = numeric(), elsewhere = numeric(), elsewhere_lines = integer())
results <- merge(results, elsewhere, by = c("split_id", "feature"), all.x = TRUE, sort = FALSE)
results <- merge(results, splits[, .(split_id, level, type, gene, label, shared_group, altered_set, kept_inside, low_information)], by = "split_id", sort = FALSE)
results[, symbol := gene_symbol(feature)]
results[, selective := feature %in% selective]

screen <- call_hits(results[selective == TRUE])
screen[, kind := fifelse(specific, "differs by driver status inside the cancer type", fifelse(inside_only %in% TRUE, "inside only", NA_character_))]
paired <- call_named(results)
paired <- merge(paired, pan_hits[, .(feature, pan_cancer_estimate = estimate, pan_cancer_p = p)], by = "feature", sort = FALSE)
# Exactly tied tests keep the order of the split table, so that the row order, and with it
# the seed of each follow-up below, does not depend on the processor.
split_position <- function(id) match(id, splits$split_id)
paired <- paired[tied_order(inside_p, split_position(split_id), first = symbol)]
write_ext(paired, "driver_type_paired_genes")
print(paired[, .(splits = .N, shown = sum(shown), median_estimate = median(inside), smallest_p = min(inside_p)), by = symbol])

# --- the same procedure under shuffled Y status -------------------------------------------------
run_null <- function(b) tail_row(rbindlist(lapply(splits$split_id, scan_split, shuffle = b, with_raw = FALSE)), sprintf("shuffled_%03d", b))
calibration <- rbind(tail_row(results, "observed"), rbindlist(parallel::mclapply(seq_len(null_runs), run_null, mc.cores = n_cores, mc.preschedule = FALSE)))
calibration[, source := "Cell lines"]
null_rows <- calibration[run != "observed"]
stopifnot(nrow(null_rows) == null_runs)
print(null_rows[, .(runs = .N, inside_below_0.05 = mean(inside_below_0.05), inside_below_0.001 = mean(inside_below_0.001), inside_lambda_mean = mean(inside_lambda), inside_lambda_min = min(inside_lambda), inside_lambda_max = max(inside_lambda),
  interaction_below_0.05 = mean(interaction_below_0.05), runs_with_inside_hit = sum(inside_fdr_hits > 0), runs_with_call = sum(specific > 0), runs_with_named = sum(named_shown > 0))])

saveRDS(screen, file.path(ext_out_dir, "driver_type_screen.rds"))
fwrite(screen, file.path(ext_out_dir, "driver_type_screen.csv.gz"))
hits <- screen[!is.na(kind)]
hits <- hits[tied_order(inside_p, split_position(split_id))]
write_ext(hits, "driver_type_hits")
split_results <- screen[, .(
  tests = .N, inside_min_p = min(inside_p), inside_lambda = lambda_of(inside_p), inside_fdr_hits = sum(inside_fdr < fdr_cut), split_alone_hits = sum(inside_fdr_split < fdr_cut),
  interaction_min_p = if (all(is.na(interaction_p))) NA_real_ else min(interaction_p, na.rm = TRUE), specific = sum(specific), inside_only = sum(inside_only),
  crit_split_alone = bh_next_critical(inside_p, fdr_cut)
), by = .(split_id, level, type, gene, label, n_neg_in, n_pos_in, n_neg_out, n_pos_out)]
write_ext(split_results, "driver_type_split_results")
log_msg("cell lines: %d tests of the Y effect among altered lines, %d of specificity; %d calls, %d inside only; shuffled runs gave %d and %d in total",
  nrow(screen), sum(is.finite(screen$interaction_p)), sum(screen$specific), sum(screen$inside_only), sum(null_rows$specific), sum(null_rows$inside_only))

# --- Sanger organoids ----------------------------------------------------------------------------------
organoids <- read_ext("organoid_cohort", dir = ext_dir)[in_primary == TRUE]
organoid_effect <- readRDS(file.path(ext_dir, "organoid_effect.rds"))
organoid_frame <- analysis_frame(organoids, organoid_effect)
organoid_ann <- read_ext("organoid_gene_annotation", dir = ext_dir)
organoid_genes <- families[source == "Sanger organoids"]
organoid_selective <- organoid_genes[selective == TRUE, symbol]
if (trial) organoid_genes <- organoid_genes[symbol %in% union(head(organoid_selective, sp$trial_genes), sp$x_partners)]
# Ranks over the organoids of the two tumour types, as in the stage-26 screen by cancer type.
type_frame <- organoid_frame[OncotreeLineage %in% op$tumour_types]
organoid_raw <- t(organoid_effect[type_frame$ModelID, organoid_genes$symbol])
organoid_rank <- rank_rows(organoid_raw)
organoid_rows <- function(id, shuffle = 0) split_rows(id, shuffle, lines_table = organoid_lines, source_frame = organoid_frame, ids = organoid_splits$split_id)
scan_organoid <- function(id, shuffle = 0, with_raw = TRUE) {
  sframe <- organoid_rows(id, shuffle)
  ranked <- fit_split(sframe, organoid_rank, dp$organoid_covariates, gate, lineage = FALSE)
  out <- ranked[, .(
    split_id = id, feature, n, n_neg_in, n_pos_in, n_neg_out, n_pos_out,
    inside_rank = inside_estimate, inside_p, inside_df, outside_rank = outside_estimate, outside_p, interaction_rank = interaction_estimate, interaction_p
  )]
  out[, c("inside_p_ols", "interaction_p_ols") := fit_split_ols(sframe, organoid_rank, dp$organoid_covariates, gate, lineage = FALSE)[, .(inside_p_ols, interaction_p_ols)]]
  if (with_raw) {
    raw <- fit_split(sframe, organoid_raw, dp$organoid_covariates, gate, lineage = FALSE, trend = TRUE)
    out[, c("inside", "inside_se", "interaction", "interaction_se") := raw[, .(inside_estimate, inside_se, interaction_estimate, interaction_se)]]
  }
  out
}
organoid_results <- rbindlist(lapply(organoid_splits$split_id, scan_organoid, shuffle = observed_shuffle))
# The organoid matrix names genes by symbol; EIF1AX is not in the organoid library.
organoid_named <- intersect(gene_symbol(named), organoid_genes$symbol)
organoid_calibration <- rbind(
  tail_row(organoid_results, "observed", keep = organoid_selective, named_features = organoid_named),
  rbindlist(parallel::mclapply(seq_len(null_runs), function(b) {
    tail_row(rbindlist(lapply(organoid_splits$split_id, scan_organoid, shuffle = b, with_raw = FALSE)), sprintf("shuffled_%03d", b), keep = organoid_selective, named_features = organoid_named)
  }, mc.cores = n_cores, mc.preschedule = FALSE))
)
organoid_calibration[, source := "Sanger organoids"]
write_ext(rbind(calibration, organoid_calibration), "driver_type_calibration")
organoid_screen <- call_hits(organoid_results[feature %in% organoid_selective])
setnames(organoid_screen, "feature", "symbol")
organoid_screen <- merge(organoid_screen, organoid_splits[, .(split_id, type, gene, label)], by = "split_id", sort = FALSE)
organoid_screen[, feature := organoid_ann$feature[match(symbol, organoid_ann$symbol)]]
organoid_screen[, kind := fifelse(specific, "differs by driver status inside the cancer type", fifelse(inside_only %in% TRUE, "inside only", NA_character_))]
write_ext(organoid_screen, "driver_type_organoid_screen")
log_msg("organoids: %d tests among altered organoids, %d of specificity; %d calls, %d inside only", nrow(organoid_screen), sum(is.finite(organoid_screen$interaction_p)), sum(organoid_screen$specific), sum(organoid_screen$inside_only))

# --- cell lines and organoids combined, for the splits both have ------------------------------------
combined <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  sframe <- split_rows(pairs$cell_split_id[i], observed_shuffle)
  pair_covariates <- covariates
  if (pairs$cell_line_msi_covariate[i]) {
    sframe <- sframe[!is.na(MSI)]
    pair_covariates <- c(covariates, "MSI")
  }
  cell <- suppressWarnings(fit_split(sframe, mat_rank, pair_covariates, gate))
  cell_unit <- unit_of(split_designs(sframe, pair_covariates, gate)$inside$design, "exposed_in")
  organoid_frame_i <- organoid_rows(pairs$organoid_split_id[i], observed_shuffle)
  organoid_unit <- unit_of(split_designs(organoid_frame_i, dp$organoid_covariates, gate, lineage = FALSE)$inside$design, "exposed_in")
  organoid <- organoid_screen[split_id == pairs$organoid_split_id[i] & !is.na(feature), .(feature, organoid_inside_rank = inside_rank, organoid_inside_p = inside_p)]
  both <- merge(cell[feature %in% selective, .(feature, cell_inside_rank = inside_estimate, cell_inside_p = inside_p)], organoid, by = "feature")
  z <- cbind(signed_z(both$cell_inside_p, both$cell_inside_rank), signed_z(both$organoid_inside_p, both$organoid_inside_rank))
  pooled <- stouffer(z, cbind(rep(1 / cell_unit, nrow(both)), rep(1 / organoid_unit, nrow(both))))
  both[, c("cell_z", "organoid_z", "combined_z", "combined_p") := .(z[, 1], z[, 2], pooled$z, pooled$p)]
  both[, `:=`(
    pair = pairs$pair[i], gene = pairs$gene[i], cell_unit = cell_unit, organoid_unit = organoid_unit, msi_adjusted_cell_lines = pairs$cell_line_msi_covariate[i],
    cell_n_negative = sum(sframe$exposed == 1 & sframe$ctx == 1), cell_n_positive = sum(sframe$exposed == 0 & sframe$ctx == 1),
    organoid_n_negative = sum(organoid_frame_i$exposed == 1 & organoid_frame_i$ctx == 1), organoid_n_positive = sum(organoid_frame_i$exposed == 0 & organoid_frame_i$ctx == 1),
    symbol = gene_symbol(feature)
  )]
  both
}))
combined[, combined_fdr := bh(combined_p)]
combined[, combined_hit := combined_fdr < fdr_cut]
setorder(combined, combined_p)
write_ext(combined, "driver_type_combined")
print(combined[, .(genes = .N, combined_hits = sum(combined_hit), min_combined_p = min(combined_p), cell_weight_share = round(mean((1 / cell_unit^2) / (1 / cell_unit^2 + 1 / organoid_unit^2)), 2)), by = .(pair, gene)])

# --- follow-up: every hit, and the paired genes shown inside a split ----------------------------------
# Checks on the rank scale by plain least squares, so that the observed and the shuffled
# statistics are computed alike. They can only downgrade a hit.
alternative_calls <- c(rna_only = "y_call_rna", dna_only = "y_call_dna", concordant_only = "y_call_concordant")
adjusted_covariates <- c(covariates, "WGD", "tp53_altered", "Aneuploidy")
targets <- unique(rbind(
  hits[, .(split_id, feature, role = "hit")],
  paired[shown == TRUE, .(split_id, feature, role = "paired gene, positive control")]
), by = c("split_id", "feature"))
# On shuffled labels nothing qualifies, so a trial follows up two rows to exercise the code.
if (trial && nrow(targets) == 0) targets <- paired[order(inside_p)][c(1, .N), .(split_id, feature, role = "trial: code exercise only")]
follow_one <- function(index) {
  id <- targets$split_id[index]
  target <- targets$feature[index]
  sframe <- split_rows(id, observed_shuffle)
  least_squares <- function(rows, fit_covariates, alone = FALSE) {
    design <- tryCatch(
      if (alone) build_design(rows[ctx == 1], tested = "exposed", covariates = fit_covariates) else split_designs(rows, fit_covariates, gate)$inside,
      error = function(e) NULL
    )
    if (is.null(design) || !design$estimable) return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, n = 0L))
    used <- (if (alone) rows[ctx == 1] else rows)[design$rows]
    tested <- if (alone) "exposed" else "exposed_in"
    fit_ols(mat_rank[target, used$ModelID], design$design[, tested], design$design[, colnames(design$design) != tested, drop = FALSE])
  }
  row_of <- function(fit, check, variant) data.table(check = check, variant = variant, estimate = fit$estimate, se = fit$se, p = fit$p, n = fit$n)
  design <- split_designs(sframe, covariates, gate)$inside
  used <- sframe[design$rows]
  base <- followup_checks(
    mat_rank[target, used$ModelID], design$design[, "exposed_in"], design$design[, colnames(design$design) != "exposed_in", drop = FALSE],
    strata = paste(used$ctx), n_perm = sp$followup_permutations, seed = dp$seed + index, tie_tolerance = 1e-9
  )
  alternatives <- rbindlist(lapply(names(alternative_calls), function(name) {
    calls <- sframe[[alternative_calls[[name]]]]
    rows <- copy(sframe)[!is.na(calls)]
    rows[, exposed := as.numeric(calls[!is.na(calls)] == "Y_negative")]
    # A trial shuffles these labels too: no variant of a trial may use a real Y call.
    if (observed_shuffle > 0) {
      set.seed(dp$seed + observed_shuffle + index + match(name, names(alternative_calls)))
      rows[, exposed := permute_within(exposed, ctx)]
    }
    row_of(least_squares(split_frame(rows, rep(TRUE, nrow(rows)), rows$ctx), covariates), "Y-status definition", name)
  }))
  one_per_donor <- sframe[order(ModelID)][!duplicated(PatientID)]
  out <- rbind(
    row_of(base, "primary", "one model per split"),
    data.table(check = "permutation within status", variant = sprintf("%d shuffles", sp$followup_permutations), estimate = base$estimate, se = NA_real_, p = base$permutation_p, n = base$n),
    data.table(check = "leave one line out", variant = "largest p", estimate = NA_real_, se = NA_real_, p = base$loo_max_p, n = base$n),
    data.table(check = "leave one line out", variant = "smallest |effect|", estimate = sign(base$estimate) * base$loo_min_abs, se = NA_real_, p = NA_real_, n = base$n),
    alternatives,
    row_of(least_squares(sframe, adjusted_covariates), "covariates", "plus WGD, TP53, aneuploidy"),
    row_of(least_squares(one_per_donor, covariates), "one line per donor", "first line of each donor"),
    row_of(least_squares(sframe, covariates, alone = TRUE), "altered lines fitted alone", "the form of the context screen")
  )
  out[, `:=`(
    split_id = id, feature = target, symbol = gene_symbol(target), role = targets$role[index],
    robust = base$permutation_p < sp$robust_permutation_p && base$loo_max_p < sp$robust_leave_one_out_p && base$loo_sign_stable
  )]
  out
}
followup <- if (nrow(targets) > 0) {
  rbindlist(parallel::mclapply(seq_len(nrow(targets)), follow_one, mc.cores = n_cores))
} else {
  data.table(check = character(), variant = character(), estimate = numeric(), se = numeric(), p = numeric(), n = integer(), split_id = character(), feature = character(), symbol = character(), role = character(), robust = logical())
}
if (nrow(followup) > 0) {
  # The organoid counterpart of a split, where the same gene was tested in the same cancer type.
  counterpart <- merge(targets, pairs[, .(split_id = cell_split_id, organoid_split_id)], by = "split_id")
  if (nrow(counterpart) > 0) {
    seen <- merge(counterpart, organoid_screen[, .(organoid_split_id = split_id, feature, estimate = inside_rank, p = inside_p, n)], by = c("organoid_split_id", "feature"))
    if (nrow(seen) > 0) {
      followup <- rbind(followup, seen[, .(check = "other source", variant = "Sanger organoids, same gene and cancer type", estimate, se = NA_real_, p, n, split_id, feature, symbol = gene_symbol(feature), role,
        robust = followup$robust[match(paste(split_id, feature), paste(followup$split_id, followup$feature))])])
    }
  }
  setcolorder(followup, c("split_id", "feature", "symbol", "role"))
}
write_ext(followup, "driver_type_hit_followup")
if (nrow(followup) > 0) print(unique(followup[, .(split_id, symbol, role, robust)]))

# --- summary ---------------------------------------------------------------------------------------------
family_bar <- function(p) bh_next_critical(p[is.finite(p)], fdr_cut)
summary_list <- list(
  lines = nrow(frame), genes = length(universe), selective = length(selective), splits = nrow(splits), specificity_splits = length(with_four),
  inside_tests = nrow(screen), interaction_tests = sum(is.finite(screen$interaction_p)),
  calls = sum(screen$specific), inside_only = sum(screen$inside_only), split_alone_hits = sum(screen$inside_fdr_split < fdr_cut),
  inside_lambda = lambda_of(screen$inside_p), interaction_lambda = lambda_of(screen$interaction_p),
  inside_min_p = min(screen$inside_p), interaction_min_p = min(screen$interaction_p, na.rm = TRUE),
  bars = list(family_inside = family_bar(screen$inside_p), family_interaction = family_bar(screen$interaction_p)),
  paired = list(tests = nrow(paired), shown = sum(paired$shown), by_gene = paired[, .(splits = .N, shown = sum(shown), median_estimate = median(inside), median_rank_effect = median(inside_rank)), by = symbol]),
  null = list(
    runs = null_runs, runs_with_call = sum(null_rows$specific > 0), runs_with_inside_only = sum(null_rows$inside_only > 0), runs_with_split_alone_hit = sum(null_rows$split_alone_hits > 0),
    runs_with_named_shown = sum(null_rows$named_shown > 0),
    inside_below_0.05 = mean(null_rows$inside_below_0.05), inside_below_0.001 = mean(null_rows$inside_below_0.001), inside_below_1e5 = mean(null_rows$inside_below_1e5),
    interaction_below_0.05 = mean(null_rows$interaction_below_0.05), interaction_below_0.001 = mean(null_rows$interaction_below_0.001),
    inside_lambda_mean = mean(null_rows$inside_lambda), inside_lambda_range = range(null_rows$inside_lambda),
    interaction_lambda_mean = mean(null_rows$interaction_lambda), interaction_lambda_range = range(null_rows$interaction_lambda),
    named_below_0.05 = mean(null_rows$named_below_0.05)
  ),
  organoids = list(
    splits = nrow(organoid_splits), inside_tests = nrow(organoid_screen), interaction_tests = sum(is.finite(organoid_screen$interaction_p)),
    calls = sum(organoid_screen$specific), inside_only = sum(organoid_screen$inside_only), inside_min_p = min(organoid_screen$inside_p),
    null_runs_with_call = organoid_calibration[run != "observed", sum(specific > 0)], null_runs_with_inside_only = organoid_calibration[run != "observed", sum(inside_only > 0)],
    null_inside_below_0.05 = organoid_calibration[run != "observed", mean(inside_below_0.05)]
  ),
  combined = list(pairs = nrow(pairs), tests = nrow(combined), hits = sum(combined$combined_hit), min_p = min(combined$combined_p)),
  least_squares = list(
    fdr_hits = sum(screen$inside_fdr_ols < fdr_cut, na.rm = TRUE), fdr_hits_paired = screen[inside_fdr_ols < fdr_cut & feature %in% named, .N],
    genes = sort(unique(screen[inside_fdr_ols < fdr_cut, symbol])), below_family_bar = sum(screen$inside_p_ols < family_bar(screen$inside_p), na.rm = TRUE), min_p = min(screen$inside_p_ols, na.rm = TRUE),
    lambda = lambda_of(screen$inside_p_ols), null_below_0.05 = mean(null_rows$ols_below_0.05), null_lambda_mean = mean(null_rows$ols_lambda), null_runs_with_hit = sum(null_rows$ols_fdr_hits > 0)
  ),
  raw_scale = list(below_family_bar = sum(screen$inside_p_raw < family_bar(screen$inside_p)), bh_hits = sum(bh(screen$inside_p_raw) < fdr_cut), bh_hits_not_paired = sum(bh(screen$inside_p_raw) < fdr_cut & !screen$feature %in% named)),
  distinct_tests = nrow(unique(screen[, .(feature, altered_set, p = signif(inside_p, 10))])),
  followed_up = nrow(targets), controls_pass = all(control_table$pass), trial = trial, labels_shuffled = trial,
  parameters_sha256 = ext_section_hash(c("shared", "driver_type"))
)
jsonlite::write_json(summary_list, file.path(ext_out_dir, "driver_type_tests_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
log_msg("driver-by-type tests done: %d calls, %d inside only; organoids %d and %d; %d combined hits", summary_list$calls, summary_list$inside_only, summary_list$organoids$calls, summary_list$organoids$inside_only, summary_list$combined$hits)
