# Stage 36: driver status inside cancer types, what each split could detect.
# The tests of stage 35 found what they found; this stage asks how large a Y-loss effect
# confined to the altered lines of one cancer type would have had to be, and how many
# times larger the collection would have to be for effects of known size.
#
# 1. Same code: the stored p-values of every split are reproduced first.
# 2. Clean shifts: real selective genes with a fixed shift in the Y-negative altered
#    lines of the split, each tested as the one gene added to the split's fit.
# 3. Benchmarks: the known vulnerabilities of stage 24 (and the four paired genes) moved
#    into those lines, as a typical gene hit that hard and as a copy of the gene itself.
# 4. Enlarged copies: the same shifts tested on 2, 4, 8 and 16 copies of the design,
#    each copy carrying the spread of a different real gene; from these, the multiple of
#    today's lines at which an effect of a reference size is found.
# 5. Check on data: the largest splits cut to a half or a quarter of their lines and then
#    copied back to full size; the limit of the copies is set against that of the real split.
# 6. Organoid and combined splits: limits from the fitted designs only (best cases).
#
# Bars of significance, each the p-value at which one more test would pass BH in the
# stored family, so they describe the layer as it stands:
#   candidate    p < 0.05 for one gene named in advance
#   split_alone  the split analysed alone, over its own genes
#   family       the Y effect among altered lines, over all splits of the layer
#   specificity  the layer's rule: family, and the difference from the unaltered lines
#                within the family of all such differences (where it can be tested)
#   named_difference  the difference from the unaltered lines for one gene named in
#                advance, p < 0.05 (where it can be tested)
#   family_raw   the family bar met by plain least squares on Chronos scores: a test
#                without the ceiling that ranks impose, and not calibrated in cells this
#                small (see 12_context_genes.R); given for comparison only
# Simulation grids, the number of host genes and the benchmarks are those of stage 24,
# read from power_parameters.json. LOY_EXTENSION_TRIAL=1 runs a few splits and genes on
# the trial output of stage 35 (use with LOY_EXTENSION_OUTPUT).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("driver_type")

power_params <- jsonlite::fromJSON(file.path(script_dir, "power_parameters.json"))
dp <- ext_params$driver_type
sp <- ext_params$shared
pw <- dp$power
alpha <- sp$fdr
targets <- power_params$power_targets
gate <- sp$min_cell
covariates <- dp$covariates
trial <- Sys.getenv("LOY_EXTENSION_TRIAL") == "1"
observed_shuffle <- if (trial) 900000 else 0
tests_summary <- jsonlite::fromJSON(file.path(ext_out_dir, "driver_type_tests_summary.json"))
stopifnot(identical(isTRUE(tests_summary$trial), trial), identical(tests_summary$parameters_sha256, ext_section_hash(c("shared", "driver_type"))))
rank_rows <- function(mat_t) {
  out <- t(apply(mat_t, 1, rank_normalise))
  dimnames(out) <- dimnames(mat_t)
  out
}
collect <- function(runs, ids, what) {
  failed <- vapply(runs, function(x) !is.list(x) || inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(what, " failed: ", paste(ids[failed], collapse = ", "), "\n", paste(unlist(runs[failed]), collapse = "\n"))
  runs
}

splits <- read_ext("driver_type_splits")[fitted == TRUE]
split_lines <- read_ext("driver_type_split_lines")
families <- read_ext("driver_type_genes")
screen <- readRDS(file.path(ext_out_dir, "driver_type_screen.rds"))
setDT(screen)
census <- read_ext("driver_type_census")

# --- the rows of stage 35 ------------------------------------------------------------------------
cohort <- read_result("cohort")
context <- read_result("context_matrix")
effect <- load_matrix("crispr_effect")
frame <- analysis_frame(cohort, effect)
cell_genes <- families[source == "Cell lines"]
selective <- cell_genes[selective == TRUE, feature]
focus <- cell_genes[focused == TRUE, feature]
universe <- cell_genes$feature
if (trial) {
  universe <- union(head(selective, sp$trial_genes), focus)
  selective <- intersect(selective, universe)
}
mat_raw <- t(effect[frame$ModelID, universe])
mat_rank <- rank_rows(mat_raw)
split_rows <- function(id, shuffle = 0) {
  lines <- split_lines[split_id == id]
  rows <- frame[match(lines$ModelID, ModelID)]
  stopifnot(nrow(lines) > 0, !anyNA(rows$ModelID), identical(as.numeric(rows$exposed), as.numeric(lines$exposed)))
  if (shuffle > 0) {
    set.seed(dp$seed + 1000 * shuffle + match(id, splits$split_id))
    rows[, exposed := permute_within(exposed, lines$ctx)]
  }
  split_frame(rows, rep(TRUE, nrow(rows)), lines$ctx)
}
crit_family <- c(inside = bh_next_critical(screen$inside_p, alpha), interaction = bh_next_critical(screen$interaction_p, alpha))
stopifnot(abs(crit_family[["inside"]] / tests_summary$bars$family_inside - 1) < 1e-10, abs(crit_family[["interaction"]] / tests_summary$bars$family_interaction - 1) < 1e-10)
log_msg("power of the driver-by-type layer: %d splits, %d rows fitted, family bars %.3g (among altered) and %.3g (difference); %d cores", nrow(splits), nrow(mat_raw), crit_family[["inside"]], crit_family[["interaction"]], n_cores)

# Line-to-line spread of each gene after lineage, library and Y status, and the host genes
# of the simulations, as in stage 24.
pan_design <- build_design(frame, tested = "exposed", covariates = covariates)
stopifnot(length(pan_design$rows) == nrow(frame))
pooled_df <- nrow(frame) - qr(pan_design$design)$rank
hosts_pool <- setdiff(selective[rowSums(is.na(mat_raw[selective, , drop = FALSE])) == 0], focus)
gene_sd <- sqrt(colSums(qr.resid(qr(pan_design$design), t(mat_raw[hosts_pool, ]))^2) / pooled_df)

multiples <- if (trial) pw$multiples[pw$multiples <= 4] else pw$multiples
n_shift_hosts <- if (trial) pw$trial_genes else power_params$shift_genes
n_enlarged <- if (trial) 20 else pw$enlarged_genes
n_cut_hosts <- if (trial) 20 else pw$cut_genes
cut_repeats <- if (trial) 2 else pw$cut_repeats
benchmark_reps <- if (trial) 20 else power_params$benchmark_replicates
stopifnot(multiples[1] == 1, n_enlarged * max(multiples) <= length(hosts_pool))
set.seed(dp$seed)
shift_hosts <- sample(hosts_pool, n_shift_hosts)
enlarged_hosts <- matrix(sample(hosts_pool, n_enlarged * max(multiples)), n_enlarged)
shift_grid <- rbind(data.table(unit = "sd", shift = power_params$shift_sd_grid), data.table(unit = "chronos", shift = power_params$shift_chronos_grid))
# Clean shifts towards more dependent, one block of host genes per grid point.
shifted_rows <- function(target, hosts, grid) {
  do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
    amount <- if (grid$unit[g] == "sd") grid$shift[g] * gene_sd[hosts] else grid$shift[g]
    spike_cell(mat_raw, hosts, target, -amount)
  }))
}

# --- benchmarks, estimated as in stage 24 -----------------------------------------------------------
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
  benchmark_fit(eligible, effect, biomarker_of(benchmarks$source[i], benchmarks$biomarker[i]), benchmarks$feature[i], covariates)
})
stored_benchmarks <- read_result("power_benchmarks")
benchmarks[, effect := vapply(benchmark_fits, `[[`, numeric(1), "coefficient")]
benchmarks[, effect_sd := effect / vapply(benchmark_fits, `[[`, numeric(1), "residual_sd")]
stopifnot(identical(benchmarks$label, stored_benchmarks$label), all(abs(benchmarks$effect - stored_benchmarks$effect) < 1e-8))
# The gene whose status defines each benchmark, to flag a benchmark moved into a split of its own driver.
benchmarks[, biomarker_gene := fifelse(source == "y_status", NA_character_, fifelse(source == "cohort", toupper(sub("_altered$", "", biomarker)), sub("_.*$", "", sub("^[a-z_]+__", "", biomarker))))]
is_y_benchmark <- benchmarks$source == "y_status"
variants <- rbind(
  CJ(benchmark = seq_len(nrow(benchmarks)), design = "typical gene", sorted = FALSE),
  CJ(benchmark = seq_len(nrow(benchmarks)), design = "own gene", sorted = FALSE)
)
if (trial) variants <- variants[benchmark %in% c(2, 11, 22)]
on_typical_genes <- function(fit, target, reps) {
  hosts <- hosts_pool[sample.int(length(hosts_pool), reps, replace = TRUE)]
  out <- spike_cell(mat_raw, hosts, target, fit$coefficient)
  # The host brings its own noise, so only the altered group's excess spread is added.
  excess <- sqrt(max(0, 1 - fit$sd_reference^2 / fit$sd_altered^2))
  if (excess > 0) out[, target] <- out[, target] + excess * fit$altered[sample.int(length(fit$altered), reps * sum(target), replace = TRUE)]
  out
}

# --- one split -------------------------------------------------------------------------------------------
bars <- c("candidate", "split_alone", "family", "specificity", "named_difference", "family_raw")
path_probes <- 30
path_tolerance <- 0.02
prepare_split <- function(id) {
  sframe <- split_rows(id, observed_shuffle)
  parts <- suppressWarnings(split_fit_parts(sframe, mat_rank, covariates, gate))
  baseline <- suppressWarnings(fit_split(sframe, mat_rank, covariates, gate))
  stored <- screen[split_id == id]
  at <- match(stored$feature, baseline$feature)
  gap <- max(abs(log10(baseline$inside_p[at]) - log10(stored$inside_p)), abs(log10(baseline$interaction_p[at]) - log10(stored$interaction_p)), na.rm = TRUE)
  if (!is.finite(gap) || gap > power_params$reproduce_log10_max) stop(sprintf("%s: stored p-values not reproduced (gap %.2g)", id, gap))
  # The simulation's own path: real genes sent through it as if each were the one gene added must
  # come back with their stored p-values, up to what one more copy of a gene does to the variance prior.
  set.seed(dp$seed + 7000 + match(id, splits$split_id))
  probe <- sample(intersect(stored$feature, hosts_pool), min(path_probes, length(intersect(stored$feature, hosts_pool))))
  through <- one_more_split_tests(parts, mat_raw[probe, , drop = FALSE])
  kept <- stored[match(probe, feature)]
  path_gap <- max(abs(log10(through$inside_p) - log10(kept$inside_p)), abs(log10(through$interaction_p) - log10(kept$interaction_p)), na.rm = TRUE)
  if (!is.finite(path_gap) || path_gap > path_tolerance) stop(sprintf("%s: the simulation path does not return the stored p-values (gap %.2g)", id, path_gap))
  list(
    id = id, sframe = sframe, parts = parts, reproduce_gap = gap, path_gap = path_gap,
    target = frame$ModelID %in% sframe[exposed == 1 & ctx == 1, ModelID], in_type = frame$ModelID %in% sframe$ModelID,
    crit = c(candidate = alpha, split_alone = bh_next_critical(stored$inside_p, alpha), family = crit_family[["inside"]], family_interaction = crit_family[["interaction"]]),
    has_specificity = !is.null(parts$interaction), inside_unit = parts$inside$unscaled, inside_df = median(baseline$inside_df[at]), df_new = parts$inside$df_new
  )
}
# `raw_p`: least-squares p-values of the same simulated genes on the Chronos scale, where they were computed.
detected <- function(results, prep, raw_p = NULL) {
  out <- results[, .(
    candidate = inside_p < prep$crit[["candidate"]], split_alone = inside_p < prep$crit[["split_alone"]], family = inside_p < prep$crit[["family"]],
    specificity = if (prep$has_specificity) interaction_p < prep$crit[["family_interaction"]] & inside_p < prep$crit[["family"]] else NA,
    named_difference = if (prep$has_specificity) interaction_p < prep$crit[["candidate"]] else NA
  )]
  out[, family_raw := if (is.null(raw_p)) NA else raw_p < prep$crit[["family"]]]
  out[]
}
# Share detected at each bar within each block of simulated genes (a grid point or a variant).
rates_by <- function(flags, block) cbind(flags, block = block)[, lapply(.SD, mean), by = block, .SDcols = bars]
benchmark_for_split <- function(i, prep) {
  fit <- benchmark_fits[[i]]
  if (!is_y_benchmark[i]) return(fit)
  # For the Y paralogs the effect and residuals come from the lines outside the cancer type.
  outside <- !(eligible$ModelID %in% frame$ModelID[prep$in_type])
  held_out <- benchmark_fit(eligible[outside], effect, biomarker_of("y_status", "y_call", eligible[outside]), benchmarks$feature[i], covariates)
  kept <- c("coefficient", "altered", "reference", "sd_altered", "sd_reference")
  fit[kept] <- held_out[kept]
  fit
}
run_split <- function(index) {
  id <- splits$split_id[index]
  prep <- prepare_split(id)
  set.seed(dp$seed + index)
  # 2. clean shifts
  shift_rows <- shifted_rows(prep$target, shift_hosts, shift_grid)
  shift_fit <- one_more_split_tests(prep$parts, shift_rows)
  shift_rates <- rates_by(detected(shift_fit, prep, ols_split_tests(prep$parts, shift_rows)$inside_p), rep(seq_len(nrow(shift_grid)), each = length(shift_hosts)))
  shift_out <- cbind(split_id = id, shift_grid[shift_rates$block], shift_rates[, ..bars])
  # 3. benchmarks
  fits <- lapply(seq_len(nrow(benchmarks)), benchmark_for_split, prep = prep)
  moved <- do.call(rbind, lapply(seq_len(nrow(variants)), function(v) {
    fit <- fits[[variants$benchmark[v]]]
    if (variants$design[v] == "typical gene") on_typical_genes(fit, prep$target, benchmark_reps) else benchmark_replicates(fit, frame, prep$target, benchmark_reps, with_spread = TRUE)
  }))
  colnames(moved) <- colnames(mat_rank)
  bench_rates <- rates_by(detected(one_more_split_tests(prep$parts, moved), prep, ols_split_tests(prep$parts, moved)$inside_p), rep(seq_len(nrow(variants)), each = benchmark_reps))
  bench_out <- cbind(split_id = id, variants[bench_rates$block], bench_rates[, ..bars])
  bench_out[, effect_used := vapply(fits, `[[`, numeric(1), "coefficient")[benchmark]]
  # 4. enlarged copies, SD grid only
  sd_grid <- data.table(unit = "sd", shift = pw$enlarged_sd_grid)
  enlarged_out <- rbindlist(lapply(multiples, function(k) {
    hosts <- as.vector(t(enlarged_hosts[, seq_len(k), drop = FALSE]))
    rows <- shifted_rows(prep$target, hosts, sd_grid)
    group <- rep(seq_len(nrow(sd_grid) * n_enlarged), each = k)
    rates <- rates_by(detected(enlarged_split_tests(prep$parts, rows, group), prep), rep(seq_len(nrow(sd_grid)), each = n_enlarged))
    cbind(split_id = id, multiple = k, shift = sd_grid$shift[rates$block], rates[, ..bars])
  }))
  log_msg("%s: %d against %d altered lines, %d simulated tests", id, sum(prep$target), prep$parts$cells[["pos_in"]], nrow(shift_fit) + nrow(moved) + nrow(enlarged_out) * n_enlarged)
  list(
    info = data.table(
      split_id = id, crit_split_alone = prep$crit[["split_alone"]], crit_family = prep$crit[["family"]], crit_family_interaction = prep$crit[["family_interaction"]],
      has_specificity = prep$has_specificity, inside_unit = prep$inside_unit, inside_df = prep$inside_df, df_new = prep$df_new, reproduce_gap = prep$reproduce_gap, path_gap = prep$path_gap
    ),
    shift = shift_out, bench = bench_out, enlarged = enlarged_out
  )
}

chosen <- seq_len(nrow(splits))
if (trial) {
  # A trial takes the largest split, two with a specificity test and one without.
  by_size <- splits[, order(-pmin(neg_alt, pos_alt))]
  chosen <- unique(c(by_size[1], head(intersect(by_size, which(splits$label == split_labels[1])), 2), head(intersect(by_size, which(splits$label == split_labels[2])), 2)))[seq_len(pw$trial_splits)]
}
runs <- collect(parallel::mclapply(chosen, run_split, mc.cores = n_cores, mc.preschedule = FALSE), splits$split_id[chosen], "split simulations")
info <- merge(splits[, .(split_id, level, type, gene, label, neg_alt, pos_alt, neg_wt, pos_wt, kept_inside, shared_group, altered_set, robust_possible, exact_floor)], rbindlist(lapply(runs, `[[`, "info")), by = "split_id", sort = FALSE)
shift_detection <- rbindlist(lapply(runs, `[[`, "shift"))
benchmark_detection <- rbindlist(lapply(runs, `[[`, "bench"))
enlarged_detection <- rbindlist(lapply(runs, `[[`, "enlarged"))
benchmark_detection <- cbind(benchmark_detection, benchmarks[benchmark_detection$benchmark, .(label, class, gene, biomarker_gene)])
benchmark_detection[, own_driver := !is.na(biomarker_gene) & biomarker_gene == info$gene[match(split_id, info$split_id)]]
benchmark_detection[, benchmark := NULL]
write_ext(shift_detection, "driver_type_shift_detection")
write_ext(benchmark_detection, "driver_type_benchmark_detection")
write_ext(enlarged_detection, "driver_type_enlarged_detection")
info[, followup_floor := permutation_floor(neg_alt, pos_alt, sp$followup_permutations)]
info[, robust_reachable := followup_floor < sp$robust_permutation_p]
log_msg("stored p-values reproduced in %d splits; largest gap in log10 p %.2g by refitting, %.2g through the simulation path", nrow(info), max(info$reproduce_gap), max(info$path_gap))

# --- limits: the shift found 50% and 80% of the time, and the best case from the design ----------------
bar_critical <- function(row, bar_name) switch(bar_name, candidate = alpha, split_alone = row$crit_split_alone, family = row$crit_family, family_raw = row$crit_family, NA_real_)
needs_four_cells <- c("specificity", "named_difference")
limits <- rbindlist(lapply(seq_len(nrow(info)), function(i) {
  row <- info[i]
  rbindlist(lapply(bars, function(bar_name) {
    if (bar_name %in% needs_four_cells && !row$has_specificity) return(NULL)
    rbindlist(lapply(c("sd", "chronos"), function(grid_unit) {
      curve <- shift_detection[split_id == row$split_id & unit == grid_unit]
      out <- detection_curve(curve$shift, curve[[bar_name]], targets)
      # Normal theory covers the Y effect among altered lines only; the rules that use the difference have
      # no closed form here. Least squares on the Chronos scale has the residual degrees of freedom alone.
      best <- if (grid_unit == "sd" && !bar_name %in% needs_four_cells) {
        vapply(targets, function(target) analytic_mde(row$inside_unit, if (bar_name == "family_raw") row$df_new else row$inside_df, bar_critical(row, bar_name), target), numeric(1))
      } else {
        rep(NA_real_, length(targets))
      }
      out[, `:=`(split_id = row$split_id, bar = bar_name, unit = grid_unit, best_case = best)]
      out
    }))
  }))
}))
limits[, ceiling := ceiling[unit == "sd"][1], by = .(split_id, bar)]
limits <- merge(info[, .(split_id, level, type, gene, label, neg_alt, pos_alt, neg_wt, pos_wt)], limits, by = "split_id", sort = FALSE)
write_ext(limits, "driver_type_limits")
ceilings <- dcast(limits[unit == "sd" & target == 0.8], split_id ~ bar, value.var = "ceiling")
setnames(ceilings, setdiff(names(ceilings), "split_id"), paste0("ceiling_", setdiff(names(ceilings), "split_id")))
info <- merge(info, ceilings, by = "split_id", sort = FALSE)
write_ext(info, "driver_type_power_splits")
print(limits[unit == "sd" & target == 0.8, .(splits = .N, reachable = sum(reachable), median_limit = median(limit, na.rm = TRUE), median_ceiling = median(ceiling), median_best_case = median(best_case, na.rm = TRUE)), by = bar])

# --- lines needed: the multiple of today's lines at which a reference effect is found ---------------------
# What "as large as a known effect" means depends on the gene it lands on. Own spread: the known effect
# over the spread of its own gene, a signal-to-noise ratio (the paired genes are far more variable than a
# typical gene). Typical gene: the same Chronos shift on a gene of typical spread, which is a larger
# number of SDs. Plain sizes are listed as well.
host_sd <- median(gene_sd)
plain_sizes <- c(0.75, 1, 1.5, 2, 3, 4, 6)
reference <- rbind(
  benchmarks[, .(reference = label, class, convention = "own spread", effect_chronos = effect, effect_sd = abs(effect_sd))],
  benchmarks[, .(reference = label, class, convention = "typical gene", effect_chronos = effect, effect_sd = abs(effect) / host_sd)],
  data.table(reference = sprintf("%g SD", plain_sizes), class = "size", convention = "size", effect_chronos = NA_real_, effect_sd = plain_sizes)
)
# The copies are tested with the rank-based test only.
enlarged_bars <- setdiff(bars, "family_raw")
enlarged_limits <- rbindlist(lapply(seq_len(nrow(info)), function(i) {
  row <- info[i]
  rbindlist(lapply(enlarged_bars, function(bar_name) {
    if (bar_name %in% needs_four_cells && !row$has_specificity) return(NULL)
    rbindlist(lapply(multiples, function(k) {
      curve <- enlarged_detection[split_id == row$split_id & multiple == k]
      cbind(split_id = row$split_id, bar = bar_name, multiple = k, detection_curve(curve$shift, curve[[bar_name]], 0.8)[, .(limit, ceiling)])
    }))
  }))
}))
write_ext(enlarged_limits, "driver_type_enlarged_limits")
# The same question by normal theory on the design scaled by the multiple, for comparison.
analytic_multiple <- function(row, bar_name, size) {
  crit <- bar_critical(row, bar_name)
  if (is.na(crit)) return(NA_real_)
  prior_df <- max(0, row$inside_df - row$df_new)
  reaches <- function(k) analytic_power(size, row$inside_unit / sqrt(k), min(k * row$df_new + (k - 1) + prior_df, 1e6), crit) - 0.8
  if (reaches(1) >= 0) return(1)
  if (reaches(1e4) < 0) return(Inf)
  uniroot(reaches, c(1, 1e4), tol = 1e-6)$root
}
lines_needed_table <- rbindlist(lapply(seq_len(nrow(info)), function(i) {
  row <- info[i]
  rbindlist(lapply(unique(enlarged_limits[split_id == row$split_id, bar]), function(bar_name) {
    curve <- enlarged_limits[split_id == row$split_id & bar == bar_name][order(multiple)]
    out <- cbind(split_id = row$split_id, bar = bar_name, reference)
    out[, multiple := lines_needed(curve$multiple, curve$limit, effect_sd)]
    # Where the limit at the next smaller size does not exist, the multiple is an upper bound.
    out[, at_most := lines_needed_bound(curve$multiple, curve$limit, effect_sd)]
    out[, multiple_normal_theory := vapply(effect_sd, function(size) analytic_multiple(row, bar_name, size), numeric(1))]
    out
  }))
}))
# A reference effect outside the simulated shifts cannot be placed on the curves.
lines_needed_table[, below_grid := effect_sd < min(pw$enlarged_sd_grid)]
lines_needed_table[, above_grid := effect_sd > max(pw$enlarged_sd_grid)]
lines_needed_table[below_grid == TRUE | above_grid == TRUE, c("multiple", "at_most") := .(NA_real_, NA)]
lines_needed_table[, beyond_largest := is.infinite(multiple)]
lines_needed_table <- merge(lines_needed_table, info[, .(split_id, level, type, gene, label, neg_alt, pos_alt, neg_wt, pos_wt)], by = "split_id", sort = FALSE)
lines_needed_table[, `:=`(neg_alt_needed = ceiling(neg_alt * multiple), pos_alt_needed = ceiling(pos_alt * multiple))]
write_ext(lines_needed_table, "driver_type_lines_needed")
print(lines_needed_table[class == "Y paralog" & bar %in% c("family", "specificity"), .(splits = .N, size_sd = round(effect_sd[1], 2), found_today = sum(multiple == 1), median_multiple = median(multiple), beyond = sum(beyond_largest), at_most = sum(at_most)), by = .(bar, convention, reference)])
stopifnot(lines_needed_table[class == "Y paralog", !any(below_grid | above_grid)])

# --- check on data: cut a real split, copy it back to full size, compare with the real split --------------
# If copies stand in for lines, a split cut to 1/k and copied k times should find what the real split finds.
sd_grid <- data.table(unit = "sd", shift = pw$enlarged_sd_grid)
cut_sets <- info[!duplicated(altered_set)]
cut_plan <- rbindlist(lapply(pw$cut_fractions, function(fraction) cut_sets[pmin(neg_alt, pos_alt) * fraction >= gate, .(split_id, fraction)]))
cut_block <- rep(seq_len(nrow(sd_grid)), each = n_cut_hosts)
cut_one <- function(index) {
  id <- cut_plan$split_id[index]
  fraction <- cut_plan$fraction[index]
  k <- round(1 / fraction)
  prep <- prepare_split(id)
  set.seed(dp$seed + 5000 + index)
  one_host <- enlarged_hosts[seq_len(n_cut_hosts), 1]
  k_hosts <- as.vector(t(enlarged_hosts[seq_len(n_cut_hosts), seq_len(k), drop = FALSE]))
  real <- rates_by(detected(one_more_split_tests(prep$parts, shifted_rows(prep$target, one_host, sd_grid)), prep), cut_block)
  runs <- rbindlist(lapply(seq_len(cut_repeats), function(r) {
    kept <- prep$sframe[, .SD[sample.int(.N, ceiling(.N * fraction))], by = .(exposed, ctx)]
    kept <- split_frame(kept, rep(TRUE, nrow(kept)), kept$ctx)
    parts <- tryCatch(suppressWarnings(split_fit_parts(kept, mat_rank, covariates, gate)), error = function(e) NULL)
    if (is.null(parts)) return(NULL)
    target <- frame$ModelID %in% kept[exposed == 1 & ctx == 1, ModelID]
    cut_prep <- prep
    cut_prep$has_specificity <- prep$has_specificity && !is.null(parts$interaction)
    alone <- rates_by(detected(one_more_split_tests(parts, shifted_rows(target, one_host, sd_grid)), cut_prep), cut_block)
    copied <- rates_by(detected(enlarged_split_tests(parts, shifted_rows(target, k_hosts, sd_grid), rep(seq_len(nrow(sd_grid) * n_cut_hosts), each = k)), cut_prep), cut_block)
    rbind(cbind(alone, what = "cut"), cbind(copied, what = "copied"))[, repeat_id := r]
  }))
  if (nrow(runs) == 0) return(list(table = NULL))
  means <- runs[, lapply(.SD, mean, na.rm = TRUE), by = .(what, block), .SDcols = bars]
  limit_of <- function(rates, bar_name) if (all(is.na(rates[[bar_name]]))) NA_real_ else detection_curve(sd_grid$shift[rates$block], rates[[bar_name]], 0.8)$limit
  list(table = rbindlist(lapply(setdiff(bars, "family_raw"), function(bar_name) {
    data.table(split_id = id, fraction = fraction, copies = k, bar = bar_name, repeats = uniqueN(runs$repeat_id),
      limit_real = limit_of(real, bar_name), limit_cut = limit_of(means[what == "cut"], bar_name), limit_copied = limit_of(means[what == "copied"], bar_name))
  })))
}
cut_check <- if (nrow(cut_plan) > 0) rbindlist(lapply(collect(parallel::mclapply(seq_len(nrow(cut_plan)), cut_one, mc.cores = n_cores, mc.preschedule = FALSE), cut_plan$split_id, "cut check"), `[[`, "table")) else data.table()
if (nrow(cut_check) > 0) {
  # Below 0: the copies find smaller shifts than the real split does, so they flatter a larger collection.
  cut_check[, copied_against_real_log2 := log2(limit_copied / limit_real)]
  cut_check[, cut_against_real_log2 := log2(limit_cut / limit_real)]
  # A limit at the smallest simulated shift is cut off by the grid and cannot be compared.
  cut_check[, at_grid_floor := (!is.na(limit_copied) & limit_copied <= min(sd_grid$shift)) | (!is.na(limit_real) & limit_real <= min(sd_grid$shift))]
  cut_check[, within_tolerance := abs(copied_against_real_log2) <= pw$cut_tolerance_log2]
  cut_check <- merge(cut_check, info[, .(split_id, level, type, gene, neg_alt, pos_alt)], by = "split_id", sort = FALSE)
}
write_ext(cut_check, "driver_type_cut_check")
print(cut_check)

# --- organoid and combined splits: best cases from the fitted designs ---------------------------------------
organoid_splits <- read_ext("driver_type_organoid_splits")
organoid_screen <- read_ext("driver_type_organoid_screen")
combined <- read_ext("driver_type_combined")
organoid_family <- bh_next_critical(organoid_screen$inside_p, alpha)
organoid_limits <- rbindlist(lapply(seq_len(nrow(organoid_splits)), function(i) {
  row <- organoid_splits[i]
  crit <- c(candidate = alpha, split_alone = bh_next_critical(organoid_screen[split_id == row$split_id, inside_p], alpha), family = organoid_family)
  data.table(
    source = "Sanger organoids", split_id = row$split_id, type = row$type, gene = row$gene, label = row$label, neg_alt = row$neg_alt, pos_alt = row$pos_alt,
    bar = names(crit), critical_p = crit, unit = row$inside_unit, df = row$residual_df,
    best_case_sd = vapply(crit, function(p_crit) analytic_mde(row$inside_unit, row$residual_df, p_crit, 0.8), numeric(1))
  )
}))
combined_limits <- combined[, .(cell_unit = cell_unit[1], organoid_unit = organoid_unit[1], crit_pair = bh_next_critical(combined_p, alpha)), by = .(pair, gene)]
combined_limits[, unit := 1 / sqrt(1 / cell_unit^2 + 1 / organoid_unit^2)]
combined_family <- bh_next_critical(combined$combined_p, alpha)
combined_limits <- rbindlist(lapply(seq_len(nrow(combined_limits)), function(i) {
  row <- combined_limits[i]
  crit <- c(candidate = alpha, split_alone = row$crit_pair, family = combined_family)
  data.table(
    source = "Cell lines and organoids", split_id = paste(row$pair, row$gene), type = row$pair, gene = row$gene, label = NA_character_, neg_alt = NA_real_, pos_alt = NA_real_,
    bar = names(crit), critical_p = crit, unit = row$unit, df = Inf, best_case_sd = vapply(crit, function(p_crit) (qnorm(1 - p_crit / 2) + qnorm(0.8)) * row$unit, numeric(1))
  )
}))
write_ext(rbind(organoid_limits, combined_limits), "driver_type_organoid_limits")

# --- how far the untested pairs are from the gate ----------------------------------------------------------
gate_distance <- census[, .(
  pairs = .N, both_y_groups = sum(pmin(neg_alt, pos_alt) >= 1), at_least_2 = sum(pmin(neg_alt, pos_alt) >= 2), at_least_3 = sum(pmin(neg_alt, pos_alt) >= 3),
  at_least_5 = sum(pmin(neg_alt, pos_alt) >= gate), at_least_10 = sum(pmin(neg_alt, pos_alt) >= 10), at_least_20 = sum(pmin(neg_alt, pos_alt) >= 20),
  all_four_cells = sum(pmin(neg_alt, pos_alt, neg_wt, pos_wt) >= gate)
), by = .(source, level)]
write_ext(gate_distance, "driver_type_gate_distance")

# --- summary ------------------------------------------------------------------------------------------------
at_80 <- limits[unit == "sd" & target == 0.8]
by_bar <- at_80[, .(splits = .N, reachable = sum(reachable), median_limit_sd = median(limit, na.rm = TRUE), smallest_limit_sd = min(limit, na.rm = TRUE), median_ceiling = median(ceiling), median_best_case_sd = median(best_case, na.rm = TRUE)), by = bar]
paired_needed <- lines_needed_table[class == "Y paralog", .(splits = .N, size_sd = effect_sd[1], found_today = sum(multiple == 1), within_16 = sum(is.finite(multiple)), median_multiple = median(multiple), smallest_multiple = min(multiple)), by = .(bar, convention, reference)]
size_needed <- lines_needed_table[class == "size", .(splits = .N, found_today = sum(multiple == 1), within_16 = sum(is.finite(multiple)), median_multiple = median(multiple), at_most = sum(at_most)), by = .(bar, reference)]
compared <- if (nrow(cut_check) > 0) cut_check[is.finite(copied_against_real_log2) & at_grid_floor == FALSE] else data.table()
summary_list <- list(
  splits = nrow(info), specificity_splits = sum(info$has_specificity), reproduce_gap = max(info$reproduce_gap),
  bars = list(family_inside = crit_family[["inside"]], family_interaction = crit_family[["interaction"]], split_alone_median = median(info$crit_split_alone)),
  simulation = list(shift_genes = n_shift_hosts, benchmark_replicates = benchmark_reps, enlarged_genes = n_enlarged, multiples = multiples, cut_genes = n_cut_hosts, cut_repeats = cut_repeats, host_pool = length(hosts_pool), host_sd_median = median(gene_sd)),
  limits_by_bar = by_bar, paired_lines_needed = paired_needed, size_lines_needed = size_needed, host_sd = host_sd, path_gap = max(info$path_gap), path_probes = path_probes,
  robust_not_reachable = sum(!info$robust_reachable),
  cut_check = if (nrow(cut_check) > 0) list(rows = nrow(cut_check), sets_halved = cut_check[fraction == 0.5, uniqueN(split_id)], sets_quartered = cut_check[fraction == 0.25, uniqueN(split_id)],
    compared = nrow(compared), within_tolerance = sum(compared$within_tolerance), by_bar = compared[, .(comparisons = .N, mean_log2 = mean(copied_against_real_log2), min_log2 = min(copied_against_real_log2), max_log2 = max(copied_against_real_log2)), by = bar]) else list(rows = 0),
  organoid_best_case_sd = organoid_limits[bar == "family", range(best_case_sd, na.rm = TRUE)], combined_best_case_sd = combined_limits[bar == "family", range(best_case_sd)],
  trial = trial, labels_shuffled = trial, power_parameters_sha256 = file_sha256(file.path(script_dir, "power_parameters.json")),
  parameters_sha256 = ext_section_hash(c("shared", "driver_type"))
)
jsonlite::write_json(summary_list, file.path(ext_out_dir, "driver_type_power_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
log_msg("driver-by-type power stage done: %d splits simulated", nrow(info))
