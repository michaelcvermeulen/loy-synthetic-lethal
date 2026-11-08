# Stage 27: context-specific Y-loss dependencies in RNAi (DEMETER2).
# The procedure of 12_context_genes.R on the RNAi matrix: for every binary
# context of the frozen panel, the Y-by-context interaction and the Y effect
# inside the context, on rank-normalised scores, BH within context class.
# RNAi genes are measured in different sets of lines, so genes are fitted in
# groups that share their lines and a gene-by-context test is kept only where
# every Y-by-context cell holds the minimum number of that gene's own lines.
# 260 of the 307 lines also have CRISPR data: agreement with CRISPR is
# replication across assays, not across independent lines.
# The four continuous contexts are not screened here (see the progress ledger).
# LOY_EXTENSION_TRIAL=1 runs a few hundred genes (use with LOY_EXTENSION_OUTPUT).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
dir.create(ext_out_dir, showWarnings = FALSE, recursive = TRUE)

rp <- ext_params$rnai
sp <- ext_params$shared
fdr_cut <- sp$fdr
trial <- Sys.getenv("LOY_EXTENSION_TRIAL") == "1"
cohort <- read_result("cohort")
context <- read_result("context_matrix")
definitions <- read_result("context_definitions")[kept == TRUE & type == "binary"]
cell_line_genes <- read_result("gene_annotation")
pan_hits <- read_result("gene_pan_hits")
rnai <- load_matrix("rnai")

frame <- analysis_frame(cohort, rnai)
scores <- rnai[frame$ModelID, ]
y_genes <- cell_line_genes[is_y_gene == TRUE, symbol]
dependent <- colSums(scores <= rp$dependent_cut, na.rm = TRUE)
not_dependent <- colSums(scores > rp$dependent_cut, na.rm = TRUE)
selective <- setdiff(colnames(scores)[dependent >= rp$selective_min_lines & not_dependent >= rp$selective_min_lines], y_genes)
focus <- intersect(union(pan_hits$symbol, params$xy_paralogs$x_partner), colnames(scores))
if (trial) selective <- head(selective, sp$trial_genes)
universe <- union(selective, focus)
mat_raw <- t(scores[, universe])
mat_rank <- t(apply(mat_raw, 1, rank_normalise))
dimnames(mat_rank) <- dimnames(mat_raw)
groups <- pattern_groups(mat_rank)
log_msg("RNAi context screen: %d lines (%d Y-negative, %d Y-positive), %d genes (%d selective, %d focus: %s) in %d line sets, %d binary contexts in the panel",
  nrow(frame), sum(frame$exposed == 1), sum(frame$exposed == 0), length(universe), length(selective), length(focus), paste(focus, collapse = ", "), length(groups), nrow(definitions))

context_value <- function(id) context[[id]][match(frame$ModelID, context$ModelID)]
supported <- vapply(definitions$context_id, function(id) context_supported(context_cells(frame$exposed, context_value(id), "binary"), "binary"), logical(1))
definitions <- definitions[supported]
log_msg("%d contexts pass the gate among RNAi lines: %s", nrow(definitions), paste(sprintf("%s %d", names(table(definitions$class)), table(definitions$class)), collapse = ", "))

# shuffle = 0 fits the real labels with the per-gene gate. shuffle > 0 shuffles Y
# status within lineage and context level; its tests are later restricted to the
# gene-by-context pairs that were tested with the real labels.
fit_context <- function(i, shuffle = 0) {
  id <- definitions$context_id[i]
  ctx_frame <- with_context(frame, context_value(id), "binary")
  if (shuffle > 0) {
    set.seed(ext_params$seed + 1000 * shuffle + i)
    ctx_frame[, exposed := permute_within(exposed, paste(OncotreeLineage, ctx))]
    ctx_frame[, interaction := exposed * ctx]
  }
  gate <- if (shuffle == 0) sp$min_cell else 1
  ranked <- suppressWarnings(fit_context_patterns(ctx_frame, mat_rank, character(), min_cell = gate, min_residual_df = sp$min_residual_df))
  out <- ranked[, .(feature, context_id = id, interaction_rank = interaction_estimate, interaction_p, inside_rank = inside_estimate, inside_p, inside_unit, n, n_neg_in, n_pos_in, n_neg_out, n_pos_out)]
  if (shuffle == 0) {
    raw <- suppressWarnings(fit_context_patterns(ctx_frame, mat_raw, character(), min_cell = gate, min_residual_df = sp$min_residual_df))
    out[, c("interaction", "interaction_se", "inside", "inside_se") := raw[, .(interaction_estimate, interaction_se, inside_estimate, inside_se)]]
  }
  out
}
run_scan <- function(shuffle) rbindlist(parallel::mclapply(seq_len(nrow(definitions)), fit_context, shuffle = shuffle, mc.cores = n_cores))
call_hits <- function(dt) {
  dt <- merge(dt, definitions[, .(context_id, class, label, type)], by = "context_id", sort = FALSE)
  dt <- add_family_fdr(dt, "interaction_p", "interaction")
  dt <- add_family_fdr(dt, "inside_p", "inside")
  dt[, context_specific := interaction_fdr_class < fdr_cut & inside_fdr_class < fdr_cut]
  dt[, detected_inside_only := !(context_specific %in% TRUE) & inside_fdr_class < fdr_cut]
  dt[]
}
cell_bin <- function(dt) cut(pmin(dt$n_neg_in, dt$n_pos_in, dt$n_neg_out, dt$n_pos_out), c(4, 6, 9, 14, 24, 49, Inf), labels = c("5-6", "7-9", "10-14", "15-24", "25-49", "50+"))
tail_table <- function(dt, run) {
  dt[, .(
    run = run, tests = .N, below_0.05 = mean(interaction_p < 0.05, na.rm = TRUE), below_0.001 = mean(interaction_p < 0.001, na.rm = TRUE),
    inside_below_0.001 = mean(inside_p < 0.001, na.rm = TRUE),
    interaction_fdr_class_hits = sum(interaction_fdr_class < fdr_cut, na.rm = TRUE), inside_fdr_class_hits = sum(inside_fdr_class < fdr_cut, na.rm = TRUE),
    context_specific = sum(context_specific, na.rm = TRUE)
  ), keyby = .(smallest_cell = cell_bin(dt))]
}

results <- run_scan(0)
results <- results[!is.na(interaction_p)]
results[, symbol := feature]
results[, selective := feature %in% selective]
results[, focused := feature %in% focus]
tested_pairs <- results[, paste(feature, context_id)]
log_msg("gene-by-context tests kept by the per-gene gate: %d (%d in the screen, %d in the focus family)", nrow(results), sum(results$selective), sum(results$focused))

# Families: the genome-wide screen holds the selective genes; the focus family the cell-line hits and X partners.
screen <- call_hits(results[selective == TRUE])
focused <- merge(results[focused == TRUE], definitions[, .(context_id, class, label, type)], by = "context_id", sort = FALSE)
focused[, interaction_fdr_focused := bh(interaction_p)]
focused[, inside_fdr_focused := bh(inside_p)]
focused[, context_specific := interaction_fdr_focused < fdr_cut & inside_fdr_focused < fdr_cut]

# --- the same procedure under shuffled Y status ---------------------------------------------------
calibration <- rbind(
  tail_table(screen, "observed"),
  rbindlist(lapply(seq_len(sp$null_runs), function(b) {
    null <- run_scan(b)
    null <- null[paste(feature, context_id) %in% tested_pairs & feature %in% selective]
    tail_table(call_hits(null), paste0("shuffled_", b))
  }))
)
write_ext(calibration, "rnai_context_calibration")
null_summary <- calibration[, .(tests = sum(tests), interaction_fdr_class_hits = sum(interaction_fdr_class_hits), inside_fdr_class_hits = sum(inside_fdr_class_hits), context_specific = sum(context_specific), below_0.001 = weighted.mean(below_0.001, tests)), by = run]
print(null_summary)

saveRDS(screen, file.path(ext_out_dir, "rnai_context_screen.rds"))
fwrite(screen, file.path(ext_out_dir, "rnai_context_screen.csv.gz"))
write_ext(focused[order(inside_p)], "rnai_context_focused")
hits <- rbind(
  screen[context_specific == TRUE][, family := "screen"],
  focused[context_specific == TRUE][, family := "focused"],
  fill = TRUE
)
write_ext(hits[order(interaction_p)], "rnai_context_hits")
write_ext(screen[detected_inside_only == TRUE][order(inside_p)], "rnai_context_inside_only")
summary_tab <- screen[, .(
  tests = .N, median_lines = as.numeric(median(n)), min_interaction_p = min(interaction_p, na.rm = TRUE), min_inside_p = min(inside_p, na.rm = TRUE),
  interaction_fdr_class_hits = sum(interaction_fdr_class < fdr_cut, na.rm = TRUE), inside_fdr_class_hits = sum(inside_fdr_class < fdr_cut, na.rm = TRUE),
  context_specific = sum(context_specific, na.rm = TRUE)
), by = .(class, context_id, label)]
write_ext(summary_tab[order(class, min_interaction_p)], "rnai_context_summary")
print(summary_tab[, .(contexts = .N, tests = sum(tests), interaction_hits = sum(interaction_fdr_class_hits), inside_hits = sum(inside_fdr_class_hits), context_specific = sum(context_specific)), by = class])
print(focused[order(inside_p)][1:min(.N, 15), .(symbol, label, class, interaction_p, inside, inside_p, inside_fdr_focused, n_neg_in, n_pos_in)])

# --- follow-up of every hit: the Y effect inside the context --------------------------------------------
# Checks may only downgrade a hit. Line sets differ by source library, so the gene's
# coverage group enters as a covariate in one of the checks.
exposures <- c(primary = "y_call", rna_only = "y_call_rna", dna_only = "y_call_dna", concordant_only = "y_call_concordant", lenient_rna_rule = "y_call_lenient")
# Which of the largest gene line sets a line belongs to: a stand-in for its source libraries.
largest <- groups[order(-vapply(groups, function(g) length(g$features), integer(1)))][seq_len(min(6, length(groups)))]
line_group <- setNames(apply(vapply(largest, function(g) frame$ModelID %in% g$lines, logical(nrow(frame))), 1, function(v) paste(as.integer(v), collapse = "")), frame$ModelID)
follow <- function(feature, id) {
  values_all <- setNames(rank_normalise(scores[, feature]), rownames(scores))
  inside_frame <- function(exposure) {
    out <- analysis_frame(cohort, rnai, exposure)
    out[, ctx := context[[id]][match(ModelID, context$ModelID)]]
    out[ctx %in% 1 & is.finite(values_all[ModelID])]
  }
  base <- inside_frame("y_call")
  design <- build_design(base, tested = "exposed")
  base <- base[design$rows]
  checks <- followup_checks(values_all[base$ModelID], base$exposed, design$design[, colnames(design$design) != "exposed", drop = FALSE], base$OncotreeLineage,
    n_perm = sp$followup_permutations, seed = ext_params$seed + 27)
  by_definition <- rbindlist(lapply(names(exposures), function(name) {
    alt <- inside_frame(exposures[[name]])
    alt_design <- build_design(alt, tested = "exposed")
    if (!alt_design$estimable) return(data.table(check = "Y-status definition", variant = name, estimate = NA_real_, p = NA_real_, n = nrow(alt)))
    alt <- alt[alt_design$rows]
    fit <- fit_ols(values_all[alt$ModelID], alt$exposed, alt_design$design[, colnames(alt_design$design) != "exposed", drop = FALSE])
    data.table(check = "Y-status definition", variant = name, estimate = fit$estimate, p = fit$p, n = fit$n)
  }))
  summary_row <- data.table(
    check = c("permutation within lineage", "leave one line out", "leave one line out"), variant = c(sprintf("%d shuffles", sp$followup_permutations), "smallest |effect|", "largest p"),
    estimate = c(checks$estimate, sign(checks$estimate) * checks$loo_min_abs, NA_real_), p = c(checks$permutation_p, NA_real_, checks$loo_max_p), n = checks$n
  )
  nuisance <- design$design[, colnames(design$design) != "exposed", drop = FALSE]
  with_coverage <- full_rank(cbind(nuisance, lineage_design(line_group[base$ModelID])[, -1, drop = FALSE]))
  covered <- fit_ols(values_all[base$ModelID], base$exposed, with_coverage)
  out <- rbind(by_definition, summary_row, data.table(check = "covariates", variant = "plus line-set group", estimate = covered$estimate, p = covered$p, n = covered$n))
  sign_holds <- function(name) isTRUE(sign(by_definition[variant == name, estimate]) == sign(checks$estimate))
  out[, robust := checks$permutation_p < sp$robust_permutation_p & checks$loo_sign_stable & checks$loo_max_p < sp$robust_leave_one_out_p & sign_holds("rna_only") & sign_holds("dna_only")]
  cbind(data.table(symbol = feature, context_id = id), out)
}
targets <- unique(rbind(hits[, .(feature, context_id)], screen[detected_inside_only == TRUE, .(feature, context_id)]))
followup <- if (nrow(targets) > 0) {
  rbindlist(lapply(seq_len(nrow(targets)), function(i) follow(targets$feature[i], targets$context_id[i])))
} else {
  data.table(symbol = character(), context_id = character(), check = character(), variant = character(), estimate = numeric(), p = numeric(), n = integer(), robust = logical())
}
write_ext(followup, "rnai_hit_followup")
if (nrow(followup) > 0) print(followup[, .(symbol, context_id, check, variant, estimate = round(estimate, 3), p = signif(p, 3), robust)][1:min(.N, 40)])

null_hits <- null_summary[run != "observed"]
jsonlite::write_json(list(
  lines = nrow(frame), y_negative = sum(frame$exposed == 1), y_positive = sum(frame$exposed == 0),
  genes_screen = length(selective), genes_focused = length(focus), line_sets = length(groups),
  contexts_fitted = nrow(definitions), contexts_by_class = as.list(table(definitions$class)),
  screen_tests = nrow(screen), focused_tests = nrow(focused),
  screen_context_specific = sum(screen$context_specific, na.rm = TRUE), focused_context_specific = sum(focused$context_specific, na.rm = TRUE),
  screen_interaction_fdr_class_hits = sum(screen$interaction_fdr_class < fdr_cut, na.rm = TRUE), inside_only = sum(screen$detected_inside_only, na.rm = TRUE),
  shuffled_interaction_fdr_class_hits = null_hits$interaction_fdr_class_hits, shuffled_inside_fdr_class_hits = null_hits$inside_fdr_class_hits, shuffled_context_specific = null_hits$context_specific,
  trial = trial, parameters_sha256 = ext_section_hash(c("shared", "rnai"))
), file.path(ext_out_dir, "rnai_context_run_summary.json"), auto_unbox = TRUE, pretty = TRUE)
log_msg("RNAi context screen done: %d contexts, %d tests in the screen, %d context-specific (screen), %d (focused), %d inside only; shuffled runs gave %s context-specific and %s inside detections",
  nrow(definitions), nrow(screen), sum(screen$context_specific, na.rm = TRUE), sum(focused$context_specific, na.rm = TRUE), sum(screen$detected_inside_only, na.rm = TRUE),
  paste(null_hits$context_specific, collapse = ", "), paste(null_hits$inside_fdr_class_hits, collapse = ", "))
