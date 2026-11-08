# Stage 13 (plan Stage 5, drugs): context-specific Y-loss drug sensitivities.
# Same three quantities as for genes (interaction, effect inside, effect
# outside), fitted per dataset. The cell gate is applied per compound, on the
# lines that actually have a response.
#
# Inference uses rank-normalised responses and the classical OLS test. The
# robust (HC3) test used for the pan-cancer scan is not reliable with 5-9 lines
# in a cell: it gave 27 times the expected number of p < 0.001 for effects
# inside the smallest contexts. Effect sizes are still reported on the raw
# scale, in units of the response SD. A within-lineage-and-context shuffle of
# Y status shows how many hits the procedure yields under the null.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
assert_freeze()

fdr_cut <- params$drug_scan$fdr
min_cell <- params$context$min_cell
min_continuous <- params$context$min_group_continuous
n_null <- 5
cohort <- read_result("cohort")
context <- read_result("context_matrix")
definitions <- read_result("context_definitions")[kept == TRUE]
pan <- read_result("drug_pan_results")
consensus <- read_result("drug_pan_consensus")
secondary_screen <- readRDS(file.path(matrix_dir, "prism_secondary_screen.rds"))
datasets <- unique(pan$dataset)

# Focused family: per-dataset hits, plus every dataset entry of a compound that passes the stacked model.
consensus_hits <- consensus[!is.na(fdr) & fdr < fdr_cut, compound_id]
focus_keys <- pan[associated == TRUE | compound_id %in% consensus_hits, paste(dataset, feature)]
log_msg("focused drug family: %d dataset-compound entries", length(focus_keys))

screen_columns <- function(feature, lines) {
  screen <- secondary_screen[[feature]][match(lines, secondary_screen$depmap_id)]
  levels <- sort(unique(na.omit(screen)))
  if (length(levels) < 2) return(NULL)
  vapply(levels[-1], function(level) as.numeric(screen == level), numeric(length(screen)))
}

# One classical OLS fit per compound; `nuisance` holds everything except the tested column.
scan_features <- function(name, y, tested, nuisance, lines, features) {
  fits <- lapply(features, function(feature) {
    z <- nuisance
    if (name == "prism_secondary") {
      extra <- screen_columns(feature, lines)
      if (!is.null(extra)) z <- cbind(z, extra)
    }
    fit <- fit_ols(y[, feature], tested, z)
    list(estimate = fit$estimate, se = fit$se, p = fit$p)
  })
  rbindlist(fits)
}

loaded <- lapply(datasets, function(name) {
  mat <- load_matrix(name)
  frame <- analysis_frame(cohort, mat)
  raw <- mat[frame$ModelID, pan[dataset == name, feature], drop = FALSE]
  ranked <- apply(raw, 2, rank_normalise)
  dimnames(ranked) <- dimnames(raw)
  list(frame = frame, raw = raw, ranked = ranked)
})
names(loaded) <- datasets

fit_context <- function(name, i, shuffle = 0) {
  frame <- loaded[[name]]$frame
  id <- definitions$context_id[i]
  type <- definitions$type[i]
  value <- context[[id]][match(frame$ModelID, context$ModelID)]
  if (!context_supported(context_cells(frame$exposed, value, type), type)) return(NULL)
  ctx_frame <- with_context(frame, value, type)
  rows <- match(ctx_frame$ModelID, frame$ModelID)
  raw <- loaded[[name]]$raw[rows, , drop = FALSE]
  ranked <- loaded[[name]]$ranked[rows, , drop = FALSE]
  has <- !is.na(raw)
  # Per-compound gate on lines that have a response.
  count <- function(mask) colSums(has[mask, , drop = FALSE])
  if (type == "binary") {
    cells <- cbind(
      n_neg_in = count(ctx_frame$exposed == 1 & ctx_frame$ctx == 1), n_pos_in = count(ctx_frame$exposed == 0 & ctx_frame$ctx == 1),
      n_neg_out = count(ctx_frame$exposed == 1 & ctx_frame$ctx == 0), n_pos_out = count(ctx_frame$exposed == 0 & ctx_frame$ctx == 0)
    )
    ok <- apply(cells, 1, min) >= min_cell
  } else {
    cells <- cbind(n_neg_in = count(ctx_frame$exposed == 1), n_pos_in = count(ctx_frame$exposed == 0), n_neg_out = NA, n_pos_out = NA)
    ok <- cells[, "n_neg_in"] >= min_continuous & cells[, "n_pos_in"] >= min_continuous
  }
  keep <- colnames(raw)[ok]
  if (length(keep) == 0) return(NULL)
  if (shuffle > 0) {
    set.seed(params$seed + 1000 * shuffle + i)
    strata <- if (type == "binary") paste(ctx_frame$OncotreeLineage, ctx_frame$ctx) else ctx_frame$OncotreeLineage
    ctx_frame[, exposed := permute_within(exposed, strata)]
    ctx_frame[, interaction := exposed * ctx]
  }
  design <- build_design(ctx_frame, tested = "interaction", covariates = c("exposed", "ctx"))
  if (!design$estimable) return(NULL)
  nuisance <- design$design[, colnames(design$design) != "interaction", drop = FALSE]
  out <- data.table(feature = keep)
  out[, c("interaction_rank", "interaction_rank_se", "interaction_p") := scan_features(name, ranked, ctx_frame$interaction, nuisance, ctx_frame$ModelID, keep)]
  if (shuffle == 0) out[, c("interaction", "interaction_se", "interaction_p_raw") := scan_features(name, raw, ctx_frame$interaction, nuisance, ctx_frame$ModelID, keep)]
  if (type == "binary") {
    for (side in if (shuffle == 0) c("inside", "outside") else "inside") {
      sub <- ctx_frame$ctx == as.numeric(side == "inside")
      sub_frame <- ctx_frame[sub]
      sub_design <- build_design(sub_frame, tested = "exposed")
      if (!sub_design$estimable) next
      sub_nuisance <- sub_design$design[, colnames(sub_design$design) != "exposed", drop = FALSE]
      rank_fit <- scan_features(name, ranked[sub, , drop = FALSE], sub_frame$exposed, sub_nuisance, sub_frame$ModelID, keep)
      out[, (paste0(side, "_p")) := rank_fit$p]
      if (shuffle == 0) {
        raw_fit <- scan_features(name, raw[sub, , drop = FALSE], sub_frame$exposed, sub_nuisance, sub_frame$ModelID, keep)
        out[, (c(side, paste0(side, "_se"))) := raw_fit[, .(estimate, se)]]
      }
    }
  }
  out <- cbind(out, as.data.table(cells[ok, , drop = FALSE]))
  out[, context_id := id]
  out[, dataset := name]
  out
}
run_scan <- function(shuffle) {
  tasks <- CJ(dataset = datasets, i = seq_len(nrow(definitions)), sorted = FALSE)
  out <- rbindlist(parallel::mclapply(seq_len(nrow(tasks)), function(k) fit_context(tasks$dataset[k], tasks$i[k], shuffle), mc.cores = n_cores), fill = TRUE)
  for (col in c("inside_p", "outside_p")) if (!col %in% names(out)) out[, (col) := NA_real_]
  out <- out[is.finite(interaction_p)]
  merge(out, definitions[, .(context_id, class, label, type)], by = "context_id", sort = FALSE)
}
call_hits <- function(dt) {
  dt <- add_family_fdr(dt, "interaction_p", "interaction")
  dt <- add_family_fdr(dt, "inside_p", "inside")
  dt[, context_specific := interaction_fdr_class < fdr_cut & (type == "continuous" | inside_fdr_class < fdr_cut)]
  dt[]
}
cell_bin <- function(dt) cut(pmin(dt$n_neg_in, dt$n_pos_in, dt$n_neg_out, dt$n_pos_out, na.rm = TRUE), c(4, 6, 9, 14, 24, 49, Inf), labels = c("5-6", "7-9", "10-14", "15-24", "25-49", "50+"))

results <- run_scan(0)
for (col in c("inside", "inside_se", "outside", "outside_se")) if (!col %in% names(results)) results[, (col) := NA_real_]
results <- merge(results, pan[, .(dataset, feature, compound_id, compound_name, mechanism, targets, response_sd)], by = c("dataset", "feature"), sort = FALSE)
results[, interaction_standardized := interaction / response_sd]
results[, inside_standardized := inside / response_sd]
results[, outside_standardized := outside / response_sd]
results[, interaction_mde_standardized := power_multiplier * interaction_se / response_sd]
results[, focused := paste(dataset, feature) %in% focus_keys]
results <- call_hits(results)
results[, detected_inside_only := !(context_specific %in% TRUE) & type == "binary" & inside_fdr_class < fdr_cut]
log_msg("observed scan: %d tests over %d contexts", nrow(results), uniqueN(results$context_id))
focused <- results[focused == TRUE]
if (nrow(focused) > 0) {
  focused[, interaction_fdr_focused := bh(interaction_p)]
  focused[, inside_fdr_focused := bh(inside_p)]
  focused[, context_specific_focused := interaction_fdr_focused < fdr_cut & (type == "continuous" | inside_fdr_focused < fdr_cut)]
}

# --- the same procedure under shuffled Y status ------------------------------
tail_table <- function(dt, run) {
  dt[, .(
    run = run, tests = .N, below_0.05 = mean(interaction_p < 0.05, na.rm = TRUE), below_0.001 = mean(interaction_p < 0.001, na.rm = TRUE),
    interaction_fdr_class_hits = sum(interaction_fdr_class < fdr_cut, na.rm = TRUE), context_specific = sum(context_specific, na.rm = TRUE)
  ), keyby = .(smallest_cell = cell_bin(dt))]
}
calibration <- rbind(
  tail_table(results, "observed"),
  rbindlist(lapply(seq_len(n_null), function(b) tail_table(call_hits(run_scan(b)), paste0("shuffled_", b))))
)
write_result(calibration, "context_drug_calibration")
null_summary <- calibration[, .(interaction_fdr_class_hits = sum(interaction_fdr_class_hits), context_specific = sum(context_specific), below_0.001 = weighted.mean(below_0.001, tests)), by = run]
print(null_summary)

saveRDS(results, file.path(derived_dir, "context_drug_screen.rds"))
fwrite(results, file.path(derived_dir, "context_drug_screen.csv.gz"))
write_result(focused[order(interaction_p)], "context_drug_focused")
write_result(results[context_specific == TRUE][order(interaction_p)], "context_drug_hits")
write_result(results[detected_inside_only == TRUE][order(inside_p)], "context_drug_inside_only")
summary_tab <- results[, .(
  tests = .N, min_interaction_p = min(interaction_p),
  interaction_fdr_class_hits = sum(interaction_fdr_class < fdr_cut, na.rm = TRUE),
  interaction_fdr_global_hits = sum(interaction_fdr_global < fdr_cut, na.rm = TRUE),
  inside_fdr_class_hits = sum(inside_fdr_class < fdr_cut, na.rm = TRUE),
  context_specific = sum(context_specific, na.rm = TRUE),
  median_interaction_mde_sd = median(interaction_mde_standardized, na.rm = TRUE)
), by = .(class, context_id, label, type)]
write_result(summary_tab[order(class, -context_specific, min_interaction_p)], "context_drug_summary")
print(summary_tab[, .(contexts = .N, tests = sum(tests), interaction_hits = sum(interaction_fdr_class_hits), inside_hits = sum(inside_fdr_class_hits), context_specific = sum(context_specific), median_mde_sd = median(median_interaction_mde_sd)), by = class])
print(results[order(interaction_p)][1:15, .(dataset, compound = substr(compound_name, 1, 22), label = substr(label, 1, 22), class, interaction_standardized, interaction_p, interaction_fdr_class, inside_p, n_neg_in, n_pos_in)])
null_hits <- null_summary[run != "observed"]
jsonlite::write_json(list(
  tests = nrow(results), contexts = uniqueN(results$context_id),
  interaction_fdr_class_hits = sum(results$interaction_fdr_class < fdr_cut, na.rm = TRUE),
  interaction_fdr_global_hits = sum(results$interaction_fdr_global < fdr_cut, na.rm = TRUE),
  context_specific = sum(results$context_specific, na.rm = TRUE),
  inside_only = sum(results$detected_inside_only, na.rm = TRUE),
  focused_tests = nrow(focused),
  focused_context_specific = if (nrow(focused) > 0) sum(focused$context_specific_focused, na.rm = TRUE) else 0,
  shuffled_interaction_fdr_class_hits = null_hits$interaction_fdr_class_hits, shuffled_context_specific = null_hits$context_specific
), file.path(derived_dir, "context_drug_summary.json"), auto_unbox = TRUE, pretty = TRUE)
log_msg("drug context scan done: %d tests, %d context-specific; shuffled runs gave %s", nrow(results), sum(results$context_specific, na.rm = TRUE), paste(null_hits$context_specific, collapse = ", "))
