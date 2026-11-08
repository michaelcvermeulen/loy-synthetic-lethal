# Stage 24: what the cancer-type tests could detect.
# The context screen (12_context_genes.R) found nothing beyond chance. This stage
# asks how large a Y-loss effect confined to one cancer type would have had to be
# for that screen to find it, at three levels of the Oncotree hierarchy (lineage,
# primary disease, subtype), and whether known vulnerabilities would have been
# found at the same sample sizes.
#
# 1. Census: Y-negative and Y-positive lines per cancer type at each level.
# 2. Clean shifts: real selective genes with a fixed shift in the Y-negative lines
#    of one cancer type, each tested as the one gene added to the screen.
# 3. Benchmarks: known genotype-dependency pairs estimated in the same cohort.
# 4. Benchmarks moved into each cancer type's Y sample sizes, two ways.
# 5. Positive controls: the four established Y-loss dependencies inside each type.
# 6. Drugs: lines per cancer type and normal-theory limits (no simulation).
#
# Bars of significance. Each is the p-value at which one more test would be
# rejected by BH in the stored family, so it describes the screen as it stands:
#   candidate     p < 0.05 for one gene named in advance
#   one_type      the cancer type analysed alone, over its own selective genes
#   class_inside  the Y effect inside the type, within the family of its class
#   as_run        the screen's rule: class_inside and the interaction, the latter
#                 also within the family of its class
# Simulation parameters are in power_parameters.json, which is not part of the
# freeze; its hash is written to power_summary.json. LOY_POWER_PARAMETERS and
# LOY_POWER_OUTPUT point a trial run at another parameter file and output folder.
# With the argument "resume" the cancer-type simulations (parts 2 and 4, the long
# step) are read back from the output folder instead of being run again.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
assert_freeze()

parameter_file <- Sys.getenv("LOY_POWER_PARAMETERS", unset = file.path(script_dir, "power_parameters.json"))
out_dir <- Sys.getenv("LOY_POWER_OUTPUT", unset = derived_dir)
power_params <- jsonlite::fromJSON(parameter_file)
resume <- "resume" %in% commandArgs(trailingOnly = TRUE)
alpha <- power_params$alpha
targets <- power_params$power_targets
cancer_classes <- power_params$cancer_type_classes
base_covariates <- c("LibraryKY", "LibraryHumagneCD")
log_msg("parameters: %s; output: %s; %d cores", parameter_file, out_dir, n_cores)

cohort <- read_result("cohort")
context <- read_result("context_matrix")
definitions <- read_result("context_definitions")
gene_ann <- read_result("gene_annotation")
effect <- load_matrix("crispr_effect")
pan_hits <- read_result("gene_pan_hits")
controls <- read_result("pipeline_controls")
screen <- readRDS(file.path(derived_dir, "context_gene_screen.rds"))
setDT(screen)

# The screen's rows: selective genes plus the focus family, rank-normalised over all lines.
focus <- union(pan_hits$feature, gene_ann[in_scan == TRUE & !is.na(y_partner), feature])
universe <- union(gene_ann[selective == TRUE, feature], focus)
frame <- analysis_frame(cohort, effect)
mat_raw <- t(effect[frame$ModelID, universe])
mat_rank <- t(apply(mat_raw, 1, rank_normalise))
dimnames(mat_rank) <- dimnames(mat_raw)
selective <- universe[universe %in% gene_ann[selective == TRUE, feature]]
log_msg("power stage: %d lines, %d rows fitted, %d selective genes", nrow(frame), nrow(mat_raw), length(selective))

# Line-to-line spread of each gene after lineage, library and Y status (raw scale).
# Host genes for the simulations: selective, measured in every line, and outside the
# focus family, whose members carry real Y effects.
pan_design <- build_design(frame, tested = "exposed", covariates = base_covariates)
stopifnot(length(pan_design$rows) == nrow(frame))
pooled_df <- nrow(frame) - qr(pan_design$design)$rank
hosts_pool <- setdiff(selective[rowSums(is.na(mat_raw[selective, ])) == 0], focus)
gene_sd <- sqrt(colSums(qr.resid(qr(pan_design$design), t(mat_raw[hosts_pool, ]))^2) / pooled_df)
sd_quantiles <- quantile(gene_sd, c(0.25, 0.5, 0.75, 0.9))
log_msg("residual SD of %d host genes: median %.3f (quartiles %.3f, %.3f; 90th percentile %.3f)", length(hosts_pool), sd_quantiles[2], sd_quantiles[1], sd_quantiles[3], sd_quantiles[4])

context_value <- function(id) context[[id]][match(frame$ModelID, context$ModelID)]
class_critical <- screen[type == "binary", .(interaction_crit = bh_next_critical(interaction_p, alpha), inside_crit = bh_next_critical(inside_p, alpha), class_tests = .N), by = class]

# --- 1. census ---------------------------------------------------------------------
level_columns <- c(lineage = "OncotreeLineage", disease = "OncotreePrimaryDisease", subtype = "OncotreeSubtype")
tested_ids <- unique(screen[class %in% cancer_classes, .(context_id, class, label)])
types <- rbindlist(lapply(names(level_columns), function(level) {
  column <- level_columns[[level]]
  frame[!is.na(get(column)) & get(column) != "", .(n_negative = sum(exposed == 1), n_positive = sum(exposed == 0)), by = .(label = get(column))][, level := level]
}))
types <- merge(types, tested_ids[, .(level = class, label, context_id)], by = c("level", "label"), all.x = TRUE)
# A subtype with the same lines as a disease was fitted once, under the disease.
same_lines <- definitions[class %in% cancer_classes & !is.na(duplicate_of) & duplicate_of != "", .(level = class, label, same_lines_as = duplicate_of)]
types <- merge(types, same_lines, by = c("level", "label"), all.x = TRUE)
gate <- params$context$min_cell
types[, smaller_group := pmin(n_negative, n_positive)]
types[, status := fifelse(!is.na(context_id), "tested",
  fifelse(smaller_group >= gate & same_lines_as %in% tested_ids$context_id, "tested under another name",
    fifelse(smaller_group == 0, "one Y group missing", "below the gate")))]
stopifnot(types[smaller_group >= gate, all(status %in% c("tested", "tested under another name"))], types[status == "tested", all(smaller_group >= gate)])
types[, floor_p := exact_floor(n_negative, n_positive)]
census <- types[, .(
  with_lines = .N, both_groups = sum(smaller_group >= 1), at_least_3 = sum(smaller_group >= 3),
  at_least_5 = sum(smaller_group >= gate), at_least_10 = sum(smaller_group >= 10),
  tested = sum(status == "tested"), tested_under_another_name = sum(status == "tested under another name"),
  lines_in_testable = sum((n_negative + n_positive)[smaller_group >= gate]), lines_total = sum(n_negative + n_positive)
), by = level][match(names(level_columns), level)]
print(census)
write_result(census, "power_census", dir = out_dir)

# A trial run may restrict the simulations to a few contexts.
keep_context <- function(ids) if (is.null(power_params$only_contexts)) rep(TRUE, length(ids)) else ids %in% power_params$only_contexts
tested_ids <- tested_ids[keep_context(context_id)]

# --- thresholds and the screen's own fit, per context ---------------------------------
# Reproducing the stored p-values proves that the simulation runs the screen's code.
# limma warns about partially missing rows in every fit, as it did in the screen.
prepare_context <- function(id) {
  value <- context_value(id)
  ctx_frame <- with_context(frame, value, "binary")
  baseline <- suppressWarnings(fit_context_tests(ctx_frame, mat_rank, base_covariates))
  stored <- screen[context_id == id]
  at <- match(stored$feature, baseline$feature)
  gap <- max(abs(log10(baseline$interaction_p[at]) - log10(stored$interaction_p)), abs(log10(baseline$inside_p[at]) - log10(stored$inside_p)), na.rm = TRUE)
  if (!is.finite(gap) || gap > power_params$reproduce_log10_max) stop(sprintf("%s: stored p-values not reproduced (gap %.2g)", id, gap))
  class_row <- class_critical[class == stored$class[1]]
  inside_design <- build_design(ctx_frame[ctx == 1], tested = "exposed", covariates = base_covariates)
  list(
    id = id, class = stored$class[1], label = stored$label[1], ctx_frame = ctx_frame, baseline = baseline,
    parts = context_fit_parts(ctx_frame, mat_rank, base_covariates),
    in_context = value %in% 1, target = frame$exposed == 1 & value %in% 1, reproduce_gap = gap,
    crit = c(candidate = alpha, one_type = bh_next_critical(stored$inside_p, alpha), class_inside = class_row$inside_crit, class_interaction = class_row$interaction_crit),
    inside_se_unit = sqrt(solve(crossprod(inside_design$design))["exposed", "exposed"]),
    inside_df = median(baseline$inside_df[at]),
    # The inside effect is estimated within lineage, so only lineages holding both Y groups inform it.
    lineages = uniqueN(ctx_frame[ctx == 1, OncotreeLineage]),
    lineages_with_both = ctx_frame[ctx == 1, .(negative = sum(exposed == 1), positive = sum(exposed == 0)), by = OncotreeLineage][negative > 0 & positive > 0, .N]
  )
}
bars <- c("candidate", "one_type", "class_inside", "as_run")
# Each simulated gene is tested as the one gene added to the screen. Simulated rows
# live on every line of the screen; the fits take the context's own lines.
simulate_context <- function(ctx, simulated_raw) {
  results <- one_more_tests(ctx$parts, simulated_raw)
  results[, `:=`(
    candidate = inside_p < ctx$crit[["candidate"]], one_type = inside_p < ctx$crit[["one_type"]],
    class_inside = inside_p < ctx$crit[["class_inside"]],
    as_run = interaction_p < ctx$crit[["class_interaction"]] & inside_p < ctx$crit[["class_inside"]]
  )]
  results
}
detection_rates <- function(results, by) results[, lapply(.SD, mean), by = by, .SDcols = bars]

# --- 3. benchmarks ---------------------------------------------------------------------
eligible <- cohort[ModelID %in% rownames(effect)]
feature_of <- function(symbol) {
  feature <- colnames(effect)[gene_symbol(colnames(effect)) == symbol]
  stopifnot(length(feature) == 1)
  feature
}
biomarker_of <- function(source, column, lines = eligible) {
  switch(source,
    context = context[[column]][match(lines$ModelID, context$ModelID)],
    cohort = as.numeric(lines[[column]]),
    y_status = as.numeric(lines[[column]] == "Y_negative"),
    stop("unknown biomarker source: ", source)
  )
}
benchmarks <- as.data.table(power_params$benchmarks)
benchmarks[, feature := vapply(gene, feature_of, character(1))]
benchmark_fits <- lapply(seq_len(nrow(benchmarks)), function(i) {
  benchmark_fit(eligible, effect, biomarker_of(benchmarks$source[i], benchmarks$biomarker[i]), benchmarks$feature[i], base_covariates)
})
benchmark_table <- cbind(benchmarks, rbindlist(lapply(benchmark_fits, function(fit) {
  data.table(
    effect = fit$coefficient, se = fit$se, p = fit$p, n_altered = fit$n_altered, n_reference = fit$n_reference,
    residual_sd = fit$residual_sd, sd_altered = fit$sd_altered, sd_reference = fit$sd_reference,
    effect_in_residual_sd = fit$coefficient / fit$residual_sd, effect_in_altered_sd = fit$coefficient / fit$sd_altered
  )
})))
benchmark_table[, small_pool := n_altered < power_params$benchmark_small_pool]
print(benchmark_table[, .(class, label, effect = round(effect, 2), p = signif(p, 2), n_altered, n_reference, sd_altered = round(sd_altered, 2), sd_reference = round(sd_reference, 2))])
write_result(benchmark_table, "power_benchmarks", dir = out_dir)

# The same estimator on the stage 08 exposures must give the stage 08 estimates.
stage08_check <- as.data.table(power_params$stage08_checks)
stage08_check[, estimate_here := vapply(seq_len(.N), function(i) {
  benchmark_fit(eligible, effect, biomarker_of(source[i], biomarker[i]), feature_of(gene[i]), base_covariates)$coefficient
}, numeric(1))]
stage08_check[, estimate_stage08 := controls[dataset == "crispr_effect"][match(stage08_check$control, control), estimate]]
print(stage08_check)
stopifnot(!anyNA(stage08_check$estimate_stage08))

# How dependent each benchmark gene already is in each tested cancer type: its fitted
# lineage and library level there, against the median line. A copy of the gene starts
# from that level, and lines already at the bottom of the ranking have little room to move.
benchmark_levels <- rbindlist(lapply(seq_len(nrow(benchmarks)), function(i) {
  level <- benchmark_nuisance(benchmark_fits[[i]], frame)
  data.table(benchmark = benchmarks$label[i], context_id = tested_ids$context_id,
    gene_level = vapply(tested_ids$context_id, function(id) mean(level[context_value(id) %in% 1]), numeric(1)) - median(level))
}))
write_result(benchmark_levels, "power_benchmark_levels", dir = out_dir)

# --- 2 and 4. simulations per cancer type ------------------------------------------------
set.seed(power_params$seed)
shift_hosts <- sample(hosts_pool, power_params$shift_genes)
shift_grid <- rbind(
  data.table(unit = "sd", shift = power_params$shift_sd_grid),
  data.table(unit = "chronos", shift = power_params$shift_chronos_grid)
)
# Clean shifts towards more dependent, one block of host genes per grid point.
shifted_rows <- function(target) {
  do.call(rbind, lapply(seq_len(nrow(shift_grid)), function(g) {
    amount <- if (shift_grid$unit[g] == "sd") shift_grid$shift[g] * gene_sd[shift_hosts] else shift_grid$shift[g]
    spike_cell(mat_raw, shift_hosts, target, -amount)
  }))
}
shift_rates <- function(results) {
  results[, grid := rep(seq_len(nrow(shift_grid)), each = length(shift_hosts))]
  rates <- detection_rates(results, "grid")
  cbind(shift_grid[rates$grid], rates[, ..bars])
}

# Two ways of moving a benchmark into a cancer type's Y-negative lines.
#   typical gene  the benchmark's mean effect, and the spread its altered group has
#                 beyond its reference group, added to real selective genes: what
#                 would happen to an ordinary gene hit that hard.
#   own gene      a resampled copy of the benchmark gene with its own lineage levels
#                 and noise: what would happen to a gene just like it.
# For the Y paralogs the effect and residuals come from the lines outside the cancer
# type, so that the prediction for a type does not use that type.
is_y_benchmark <- benchmarks$source == "y_status"
benchmark_for_context <- function(i, ctx) {
  fit <- benchmark_fits[[i]]
  if (!is_y_benchmark[i]) return(fit)
  outside <- !(eligible$ModelID %in% frame$ModelID[ctx$in_context])
  held_out <- benchmark_fit(eligible[outside], effect, biomarker_of("y_status", "y_call", eligible[outside]), benchmarks$feature[i], base_covariates)
  kept <- c("coefficient", "altered", "reference", "sd_altered", "sd_reference")
  fit[kept] <- held_out[kept]
  fit
}
on_typical_genes <- function(fit, target, reps, with_spread) {
  hosts <- hosts_pool[sample.int(length(hosts_pool), reps, replace = TRUE)]
  out <- spike_cell(mat_raw, hosts, target, fit$coefficient)
  # The host brings its own noise, so only the altered group's excess spread is added.
  excess <- sqrt(max(0, 1 - fit$sd_reference^2 / fit$sd_altered^2))
  if (with_spread && excess > 0) {
    out[, target] <- out[, target] + excess * fit$altered[sample.int(length(fit$altered), reps * sum(target), replace = TRUE)]
  }
  out
}
# "own gene, no effect" is the null of the own-gene design: the same resampled gene
# with nothing added in the cell. It measures how often each benchmark gene passes by chance.
variants <- rbind(
  CJ(benchmark = seq_len(nrow(benchmarks)), design = "typical gene", with_spread = c(TRUE, FALSE), sorted = FALSE),
  CJ(benchmark = seq_len(nrow(benchmarks)), design = "own gene", with_spread = TRUE, sorted = FALSE),
  CJ(benchmark = seq_len(nrow(benchmarks)), design = "own gene, no effect", with_spread = FALSE, sorted = FALSE)
)

run_cancer_type <- function(index) {
  id <- tested_ids$context_id[index]
  set.seed(power_params$seed + index)
  ctx <- prepare_context(id)
  shift_fit <- simulate_context(ctx, shifted_rows(ctx$target))
  reps <- power_params$benchmark_replicates
  fits <- lapply(seq_len(nrow(benchmarks)), benchmark_for_context, ctx = ctx)
  moved <- do.call(rbind, lapply(seq_len(nrow(variants)), function(v) {
    fit <- fits[[variants$benchmark[v]]]
    if (variants$design[v] == "own gene, no effect") fit$coefficient <- 0
    if (variants$design[v] == "typical gene") {
      on_typical_genes(fit, ctx$target, reps, variants$with_spread[v])
    } else {
      benchmark_replicates(fit, frame, ctx$target, reps, variants$with_spread[v])
    }
  }))
  colnames(moved) <- colnames(mat_rank)
  bench_fit <- simulate_context(ctx, moved)
  bench_fit[, variant := rep(seq_len(nrow(variants)), each = reps)]
  rates <- detection_rates(bench_fit, "variant")
  bench_rates <- cbind(variants[rates$variant], rates[, ..bars])
  bench_rates[, benchmark_label := benchmarks$label[benchmark]]
  bench_rates[, effect_used := fifelse(design == "own gene, no effect", 0, vapply(fits, `[[`, numeric(1), "coefficient")[benchmark])]
  log_msg("%s: %d vs %d lines, %d simulated genes", id, sum(ctx$target), sum(ctx$in_context & frame$exposed == 0), nrow(shift_fit) + nrow(bench_fit))
  list(
    info = data.table(
      context_id = id, level = ctx$class, label = ctx$label, n_negative = sum(ctx$target), n_positive = sum(ctx$in_context & frame$exposed == 0),
      crit_one_type = ctx$crit[["one_type"]], crit_class_inside = ctx$crit[["class_inside"]], crit_class_interaction = ctx$crit[["class_interaction"]],
      inside_se_unit = ctx$inside_se_unit, inside_df = ctx$inside_df, reproduce_gap = ctx$reproduce_gap
    ),
    shift = cbind(context_id = id, shift_rates(shift_fit)),
    bench = cbind(context_id = id, bench_rates[, !"benchmark"])
  )
}
collect <- function(runs, ids, what) {
  failed <- vapply(runs, function(x) !is.list(x) || inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(what, " failed: ", paste(ids[failed], collapse = ", "), "\n", paste(unlist(runs[failed]), collapse = "\n"))
  runs
}
info_columns <- c("context_id", "level", "label", "n_negative", "n_positive", "crit_one_type", "crit_class_inside", "crit_class_interaction", "inside_se_unit", "inside_df", "reproduce_gap")
if (resume) {
  shift_detection <- read_result("power_shift_detection", dir = out_dir)
  benchmark_detection <- read_result("power_benchmark_detection", dir = out_dir)
  context_info <- read_result("power_cancer_types", dir = out_dir)[status == "tested" & context_id %in% tested_ids$context_id, ..info_columns]
  stopifnot(setequal(context_info$context_id, tested_ids$context_id), setequal(shift_detection$context_id, tested_ids$context_id), setequal(benchmark_detection$context_id, tested_ids$context_id))
  type_columns <- context_info[, .(context_id, level, label, n_negative, n_positive)]
  log_msg("resumed: simulations of %d cancer types read back from %s", nrow(context_info), out_dir)
} else {
  cancer_runs <- collect(parallel::mclapply(seq_len(nrow(tested_ids)), run_cancer_type, mc.cores = n_cores, mc.preschedule = FALSE), tested_ids$context_id, "cancer-type simulations")
  context_info <- rbindlist(lapply(cancer_runs, `[[`, "info"))[, ..info_columns]
  type_columns <- context_info[, .(context_id, level, label, n_negative, n_positive)]
  shift_detection <- merge(type_columns, rbindlist(lapply(cancer_runs, `[[`, "shift")), by = "context_id")
  benchmark_detection <- merge(type_columns, rbindlist(lapply(cancer_runs, `[[`, "bench")), by = "context_id")
  setnames(benchmark_detection, "benchmark_label", "benchmark")
  write_result(shift_detection, "power_shift_detection", dir = out_dir)
  write_result(benchmark_detection, "power_benchmark_detection", dir = out_dir)
}

# Limits: the shift giving 50% and 80% detection, the ceiling, and the normal-theory limit.
bar_critical <- function(info, bar_name) switch(bar_name, candidate = alpha, one_type = info$crit_one_type, class_inside = info$crit_class_inside)
limits <- rbindlist(lapply(seq_len(nrow(context_info)), function(i) {
  info <- context_info[i]
  rbindlist(lapply(bars, function(bar_name) {
    rbindlist(lapply(c("sd", "chronos"), function(grid_unit) {
      curve <- shift_detection[context_id == info$context_id & unit == grid_unit]
      out <- detection_curve(curve$shift, curve[[bar_name]], targets)
      # Normal theory covers the inside comparison only; the two-part rule has no closed form.
      best <- if (grid_unit == "sd" && bar_name != "as_run") {
        vapply(targets, function(target) analytic_mde(info$inside_se_unit, info$inside_df, bar_critical(info, bar_name), target), numeric(1))
      } else {
        rep(NA_real_, length(targets))
      }
      out[, `:=`(context_id = info$context_id, bar = bar_name, unit = grid_unit, best_case = best)]
      out
    }))
  }))
}))
# The ceiling is detection at the largest shift of the SD grid.
limits[, ceiling := ceiling[unit == "sd"][1], by = .(context_id, bar)]
limits <- merge(type_columns, limits, by = "context_id")
write_result(limits, "power_limits", dir = out_dir)

# One row per cancer type, tested or not: counts, exact-test floor, normal-theory limits.
types <- merge(types, context_info[, !c("level", "label", "n_negative", "n_positive")], by = "context_id", all.x = TRUE)
untested <- types$status != "tested" & types$smaller_group >= power_params$floor_min_group
types[untested, inside_se_unit := sqrt(1 / n_negative + 1 / n_positive)]
types[untested, inside_df := power_params$pooled_df_for_untested]
single_discovery <- alpha / (length(selective) + 1)
types[, bar_one_type := fifelse(is.na(crit_one_type), single_discovery, crit_one_type)]
types[, best_case_candidate_sd := mapply(function(se, df) if (is.na(se)) NA_real_ else analytic_mde(se, df, alpha, 0.8), inside_se_unit, inside_df)]
types[, best_case_one_type_sd := mapply(function(se, df, crit) if (is.na(se)) NA_real_ else analytic_mde(se, df, crit, 0.8), inside_se_unit, inside_df, bar_one_type)]
types[, floor_reaches_one_type := floor_p < bar_one_type]
setorder(types, level, -smaller_group, -n_negative)
write_result(types, "power_cancer_types", dir = out_dir)

# --- 5. positive controls ----------------------------------------------------------------
positive <- screen[class %in% cancer_classes][, own_fdr := bh(inside_p), by = context_id][feature %in% pan_hits$feature & context_id %in% context_info$context_id,
  .(context_id, level = class, label, symbol, n_negative = n_neg_in, n_positive = n_pos_in, inside, inside_p, own_fdr, inside_fdr_class, interaction_p, interaction_fdr_class)]
positive[, result := fifelse(inside_fdr_class < alpha, "significant in the screen's family",
  fifelse(own_fdr < alpha, "genome-wide in this cancer type alone", fifelse(inside_p < alpha, "nominal only", "not nominal")))]
predicted <- benchmark_detection[design == "own gene" & benchmark %in% benchmarks$label[is_y_benchmark],
  .(context_id, symbol = sub(" .*$", "", benchmark), effect_held_out = effect_used, predicted_candidate = candidate, predicted_one_type = one_type, predicted_class_inside = class_inside)]
positive <- merge(positive, predicted, by = c("context_id", "symbol"), all.x = TRUE)
stopifnot(!anyNA(positive$predicted_one_type))
# The real effects sit in every Y-negative line, not in one cancer type, and that changes
# the ranks. For a prediction like for like, each gene is resampled with its effect in all
# Y-negative lines (estimated without the cancer type's own lines) and tested inside the type.
predict_everywhere <- function(index) {
  id <- tested_ids$context_id[index]
  set.seed(power_params$seed + 2000 + index)
  ctx <- prepare_context(id)
  reps <- power_params$benchmark_replicates
  y_rows <- which(is_y_benchmark)
  moved <- do.call(rbind, lapply(y_rows, function(i) benchmark_replicates(benchmark_for_context(i, ctx), frame, frame$exposed == 1, reps, with_spread = TRUE)))
  fit <- simulate_context(ctx, moved)
  fit[, symbol := rep(benchmarks$gene[y_rows], each = reps)]
  cbind(context_id = id, fit[, .(predicted_everywhere_candidate = mean(candidate), predicted_everywhere_one_type = mean(one_type), predicted_everywhere_class_inside = mean(class_inside)), by = symbol])
}
everywhere <- rbindlist(collect(parallel::mclapply(seq_len(nrow(tested_ids)), predict_everywhere, mc.cores = n_cores, mc.preschedule = FALSE), tested_ids$context_id, "positive-control predictions"))
positive <- merge(positive, everywhere, by = c("context_id", "symbol"), all.x = TRUE)
stopifnot(!anyNA(positive$predicted_everywhere_one_type))
write_result(positive[order(level, -pmin(n_negative, n_positive), symbol)], "power_positive_controls", dir = out_dir)
positive_summary <- positive[, .(
  chances = .N, nominal = sum(inside_p < alpha), one_type = sum(own_fdr < alpha), class_inside = sum(inside_fdr_class < alpha),
  types_recovering_none = uniqueN(context_id) - uniqueN(context_id[own_fdr < alpha]),
  predicted_one_type = sum(predicted_one_type), predicted_class_inside = sum(predicted_class_inside),
  predicted_everywhere_one_type = sum(predicted_everywhere_one_type), predicted_everywhere_class_inside = sum(predicted_everywhere_class_inside)
), by = level]
print(positive_summary)

# Every gene that passes with a cancer type analysed alone: the bar the benchmarks are held to,
# applied to the real data.
alone_hits <- screen[class %in% cancer_classes][, own_fdr := bh(inside_p), by = context_id][own_fdr < alpha,
  .(context_id, level = class, label, feature, symbol, n_negative = n_neg_in, n_positive = n_pos_in, inside, inside_p, own_fdr, inside_fdr_class, interaction_p, interaction_fdr_class, established = feature %in% pan_hits$feature)]
write_result(alone_hits[order(level, label, inside_p)], "power_alone_hits", dir = out_dir)

# --- the fast path against refitting the whole screen ----------------------------------------
# one_more_tests() computes what limma would give with one gene added. For a sample of
# simulated genes in a few contexts, the whole screen is refitted with that one gene
# appended (embed_and_fit with one gene per fit) and the two p-values are compared.
reference_check <- NULL
reference_ids <- intersect(power_params$reference_check_contexts, c(context_info$context_id, screen[type == "binary", unique(context_id)]))
reference_ids <- reference_ids[keep_context(reference_ids)]
if (length(reference_ids) > 0) {
  gaps <- parallel::mclapply(seq_along(reference_ids), function(k) {
    set.seed(power_params$seed + 900 + k)
    ctx <- prepare_context(reference_ids[k])
    half <- power_params$reference_check_rows %/% 2
    sample_rows <- rbind(
      spike_cell(mat_raw, sample(hosts_pool, half), ctx$target, -sample(power_params$shift_chronos_grid, half, replace = TRUE)),
      benchmark_replicates(benchmark_fits[[sample.int(nrow(benchmarks), 1)]], frame, ctx$target, half, with_spread = TRUE)
    )
    fast <- one_more_tests(ctx$parts, sample_rows)
    slow <- suppressWarnings(embed_and_fit(mat_rank, sample_rows, ctx$ctx_frame, base_covariates, batch = 1, baseline = ctx$baseline))$results
    max(abs(log10(fast$inside_p) - log10(slow$inside_p)), abs(log10(fast$interaction_p) - log10(slow$interaction_p)))
  }, mc.cores = n_cores)
  failed <- !vapply(gaps, is.numeric, logical(1))
  if (any(failed)) stop("reference check failed: ", paste(reference_ids[failed], collapse = ", "), "\n", paste(unlist(gaps[failed]), collapse = "\n"))
  gaps <- unlist(gaps)
  reference_check <- list(contexts = reference_ids, genes_per_context = 2 * (power_params$reference_check_rows %/% 2), largest_gap = max(gaps))
  log_msg("fast path against refitting the screen: largest difference in log10 p %.2g over %d contexts", reference_check$largest_gap, length(reference_ids))
}

# --- supplementary: clean shifts in the genotype-defined contexts ----------------------------
genotype_ids <- unique(screen[type == "binary" & !class %in% cancer_classes, .(context_id, class, label, n_neg_in, n_pos_in)])
genotype_ids <- genotype_ids[keep_context(context_id)]
genotype_limits <- NULL
if (nrow(genotype_ids) > 0) {
  set.seed(power_params$seed + 500)
  genotype_hosts <- sample(hosts_pool, power_params$genotype_genes)
  run_genotype <- function(index) {
    id <- genotype_ids$context_id[index]
    # Where the screen could not fit the Y effect inside the context, it could make no call there.
    if (!screen[context_id == id, any(is.finite(inside_p))]) {
      return(data.table(target = rep(targets, 2), limit = NA_real_, reachable = FALSE, ceiling = 0, context_id = id, bar = rep(c("one_type", "as_run"), each = length(targets)), reproduce_gap = NA_real_, inside_fitted = FALSE,
        inside_se_unit = NA_real_, lineages = NA_integer_, lineages_with_both = NA_integer_))
    }
    set.seed(power_params$seed + 1000 + index)
    ctx <- prepare_context(id)
    grid <- power_params$shift_chronos_grid
    fit <- simulate_context(ctx, do.call(rbind, lapply(grid, function(amount) spike_cell(mat_raw, genotype_hosts, ctx$target, -amount))))
    fit[, grid := rep(seq_along(grid), each = length(genotype_hosts))]
    rates <- detection_rates(fit, "grid")
    rbindlist(lapply(c("one_type", "as_run"), function(bar_name) {
      out <- detection_curve(grid[rates$grid], rates[[bar_name]], targets)
      out[, `:=`(context_id = id, bar = bar_name, reproduce_gap = ctx$reproduce_gap, inside_fitted = TRUE, inside_se_unit = ctx$inside_se_unit, lineages = ctx$lineages, lineages_with_both = ctx$lineages_with_both)]
      out
    }))
  }
  genotype_runs <- collect(parallel::mclapply(seq_len(nrow(genotype_ids)), run_genotype, mc.cores = n_cores, mc.preschedule = FALSE), genotype_ids$context_id, "genotype-context simulations")
  genotype_limits <- merge(genotype_ids, rbindlist(genotype_runs), by = "context_id")
  setnames(genotype_limits, c("limit", "ceiling"), c("limit_chronos", "detection_at_largest_shift"))
  write_result(genotype_limits[order(class, bar, target, limit_chronos)], "power_genotype_context_limits", dir = out_dir)
}

# --- 6. drugs ----------------------------------------------------------------------------
drug_ann <- read_result("drug_feature_annotation")
drug_screen <- readRDS(file.path(derived_dir, "context_drug_screen.rds"))
setDT(drug_screen)
drug_family <- drug_screen[class == "lineage", .(compounds = uniqueN(feature)), by = dataset]
braf_controls <- controls[grepl("^BRAF altered -> BRAF inhibitor", control)]
drug_limits <- rbindlist(lapply(drug_family$dataset, function(name) {
  mat <- load_matrix(name)
  counts <- analysis_frame(cohort, mat)[, .(n_negative = sum(exposed == 1), n_positive = sum(exposed == 0)), by = .(lineage = OncotreeLineage)]
  family_size <- drug_family[dataset == name, compounds]
  counts[, `:=`(dataset = name, compounds = family_size, testable = pmin(n_negative, n_positive) >= gate)]
  # The drug screen used a classical test on each context's own lines, so the degrees of freedom are its own.
  counts[testable == TRUE, best_case_candidate_sd := mapply(function(a, b) analytic_mde(sqrt(1 / a + 1 / b), a + b - 2, alpha, 0.8), n_negative, n_positive)]
  counts[testable == TRUE, best_case_dataset_sd := mapply(function(a, b) analytic_mde(sqrt(1 / a + 1 / b), a + b - 2, alpha / family_size, 0.8), n_negative, n_positive)]
  counts[, floor_p := exact_floor(n_negative, n_positive)]
  counts[, floor_reaches_dataset := floor_p < alpha / family_size]
  # The known drug benchmark in the same units: BRAF inhibitor response in BRAF-altered lines.
  control <- braf_controls[dataset == name]
  if (nrow(control) == 1 && control$feature %in% colnames(mat)) {
    lines <- cohort[ModelID %in% rownames(mat)]
    fit <- benchmark_fit(lines, mat, context$driver__BRAF[match(lines$ModelID, context$ModelID)], control$feature)
    counts[, `:=`(braf_inhibitor = drug_ann[dataset == name & feature == control$feature, compound_name][1], braf_effect_sd = fit$coefficient / fit$residual_sd, braf_n_altered = fit$n_altered)]
  }
  counts
}), fill = TRUE)
setcolorder(drug_limits, c("dataset", "lineage"))
write_result(drug_limits[order(dataset, -pmin(n_negative, n_positive))], "power_drug_limits", dir = out_dir)

# --- summary -------------------------------------------------------------------------------
finite_or_na <- function(x, f) if (any(is.finite(x))) f(x[is.finite(x)]) else NA_real_
limit_summary <- limits[target == 0.8, .(
  types = .N, reachable = sum(reachable), median_limit = finite_or_na(limit, median), smallest_limit = finite_or_na(limit, min), largest_limit = finite_or_na(limit, max),
  never_reach = sum(!reachable), lowest_ceiling = min(ceiling)
), by = .(level, bar, unit)]
benchmark_summary <- benchmark_detection[, .(
  types = .N, at_least_80 = sum(one_type >= 0.8), at_least_50 = sum(one_type >= 0.5), as_run_at_least_80 = sum(as_run >= 0.8), best = max(one_type), best_as_run = max(as_run)
), by = .(benchmark, design, with_spread, level)]
benchmark_summary <- merge(benchmarks[, .(benchmark = label, class)], benchmark_summary, by = "benchmark", sort = FALSE)
write_result(benchmark_summary, "power_benchmark_summary", dir = out_dir)
null_shift <- shift_detection[unit == "sd" & shift == 0, .(candidate = mean(candidate), one_type = max(one_type), as_run = max(as_run))]
# Each benchmark gene with no effect, pooled over cancer types.
null_benchmark <- benchmark_detection[design == "own gene, no effect", .(candidate = mean(candidate), one_type = mean(one_type), class_inside = mean(class_inside), as_run = mean(as_run)), by = benchmark]
write_result(null_benchmark, "power_benchmark_null", dir = out_dir)
rising <- shift_detection[order(context_id, unit, shift), .(ok = all(diff(one_type) > -0.05)), by = .(context_id, unit)]
all_gaps <- c(context_info$reproduce_gap, genotype_limits$reproduce_gap)
checks <- list(
  stored_p_values_reproduced = max(all_gaps, na.rm = TRUE) <= power_params$reproduce_log10_max,
  fast_path_equals_refitting_the_screen = !is.null(reference_check) && reference_check$largest_gap <= power_params$reference_log10_max,
  null_shift_near_alpha = abs(null_shift$candidate - alpha) < 0.02,
  null_shift_never_genome_wide = null_shift$one_type <= 0.01 && null_shift$as_run <= 0.01,
  null_benchmark_near_alpha = abs(mean(null_benchmark$candidate) - alpha) < 0.02,
  null_benchmark_never_genome_wide = max(null_benchmark$one_type) <= 0.01 && max(null_benchmark$as_run) <= 0.01,
  curves_rise_with_shift = all(rising$ok),
  stage08_estimates_match = all(abs(stage08_check$estimate_here - stage08_check$estimate_stage08) < power_params$stage08_tolerance)
)
print(unlist(checks))
jsonlite::write_json(list(
  lines = nrow(frame), rows_fitted = nrow(mat_raw), selective_genes = length(selective), host_genes = length(hosts_pool), pooled_df = pooled_df,
  gene_sd = as.list(setNames(as.numeric(sd_quantiles), c("q25", "median", "q75", "q90"))),
  class_critical = class_critical, census = census, positive_controls = positive_summary, limit_summary = limit_summary,
  cancer_types_simulated = nrow(context_info), genotype_contexts_simulated = nrow(genotype_ids),
  simulation = list(shift_genes = power_params$shift_genes, benchmark_replicates = power_params$benchmark_replicates, genotype_genes = power_params$genotype_genes),
  max_reproduce_gap = max(all_gaps, na.rm = TRUE), resumed = resume, reference_check = reference_check,
  null_shift = as.list(null_shift),
  null_benchmark = list(candidate_mean = mean(null_benchmark$candidate), candidate_range = range(null_benchmark$candidate), one_type_max = max(null_benchmark$one_type), as_run_max = max(null_benchmark$as_run)),
  stage08 = stage08_check, checks = checks,
  parameters_sha256 = file_sha256(parameter_file)
), file.path(out_dir, "power_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 8)
if (!all(unlist(checks))) stop("power stage checks failed: ", paste(names(checks)[!unlist(checks)], collapse = ", "))
log_msg("power stage done: %d cancer types and %d genotype contexts simulated", nrow(context_info), nrow(genotype_ids))
