# Stage 12 (plan Stage 5, genes): context-specific Y-loss dependencies.
# For every supported context: the Y-by-context interaction, the Y effect inside
# the context and the Y effect outside it, on selective dependencies plus the
# pan-cancer hits and X-Y partners (the focused family).
#
# Inference uses rank-normalised Chronos scores. Selective dependencies are
# heavy-tailed, and with 5-9 lines in a cell a single strongly dependent line
# drives a normal-theory test: on the raw scale the smallest cells gave 2.6
# times the expected number of p < 0.001. Effect sizes are still reported on
# the raw Chronos scale. A within-lineage-and-context shuffle of Y status is run
# alongside to show how many hits the procedure yields under the null.
# The screen is fitted without donor blocking (11 multi-model donors of 557).
# Hits are re-examined in 19_hit_robustness.R (permutation, leave-one-out,
# alternative Y-status definitions, RNAi), also without donor blocking.
# A continuous context has no "inside" group, so its calls rest on the
# interaction alone; a binary context also needs the effect inside it.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
assert_freeze()

fdr_cut <- params$gene_scan$fdr
cohort <- read_result("cohort")
context <- read_result("context_matrix")
definitions <- read_result("context_definitions")[kept == TRUE]
gene_ann <- read_result("gene_annotation")
effect <- load_matrix("crispr_effect")
pan_hits <- read_result("gene_pan_hits")
base_covariates <- c("LibraryKY", "LibraryHumagneCD")
n_null <- 5

focus <- union(pan_hits$feature, gene_ann[in_scan == TRUE & !is.na(y_partner), feature])
universe <- union(gene_ann[selective == TRUE, feature], focus)
frame <- analysis_frame(cohort, effect)
mat_raw <- t(effect[frame$ModelID, universe])
mat_rank <- t(apply(mat_raw, 1, rank_normalise))
dimnames(mat_rank) <- dimnames(mat_raw)
log_msg("context scan: %d lines, %d genes (%d selective, %d focus), %d candidate contexts", nrow(frame), length(universe), sum(gene_ann$selective), length(focus), nrow(definitions))

scan_subset <- function(sub_frame, tested, covariates, mat, trend) {
  design <- build_design(sub_frame, tested = tested, covariates = covariates)
  # Needs an estimable coefficient and a few residual degrees of freedom.
  if (!design$estimable || nrow(design$design) - ncol(design$design) < 3) return(NULL)
  sub_frame <- sub_frame[design$rows]
  scan <- tryCatch(limma_scan(mat[, sub_frame$ModelID, drop = FALSE], design$design, tested, trend = trend), error = function(e) {
    log_msg("fit skipped (%s): %s", tested, conditionMessage(e))
    NULL
  })
  if (is.null(scan)) return(NULL)
  scan[, se := estimate / t]
  scan
}

# shuffle = 0 fits the real labels; shuffle > 0 shuffles Y status within lineage
# and context level first, which keeps every cell count and removes any Y effect.
fit_context <- function(i, shuffle = 0) {
  id <- definitions$context_id[i]
  type <- definitions$type[i]
  value <- context[[id]][match(frame$ModelID, context$ModelID)]
  cells <- context_cells(frame$exposed, value, type)
  if (!context_supported(cells, type)) return(NULL)
  ctx_frame <- with_context(frame, value, type)
  if (shuffle > 0) {
    set.seed(params$seed + 1000 * shuffle + i)
    strata <- if (type == "binary") paste(ctx_frame$OncotreeLineage, ctx_frame$ctx) else ctx_frame$OncotreeLineage
    ctx_frame[, exposed := permute_within(exposed, strata)]
    ctx_frame[, interaction := exposed * ctx]
  }
  interaction_cov <- c("exposed", "ctx", base_covariates)
  ranked <- scan_subset(ctx_frame, "interaction", interaction_cov, mat_rank, trend = FALSE)
  if (is.null(ranked)) return(NULL)
  out <- ranked[, .(feature, interaction_rank = estimate, interaction_p = p, n = n)]
  if (type == "binary") {
    inside_rank <- scan_subset(ctx_frame[ctx == 1], "exposed", base_covariates, mat_rank, trend = FALSE)
    if (!is.null(inside_rank)) out[, c("inside_rank", "inside_p") := inside_rank[, .(estimate, p)]]
  }
  if (shuffle == 0) {
    raw <- scan_subset(ctx_frame, "interaction", interaction_cov, mat_raw, trend = TRUE)
    if (!is.null(raw)) out[, c("interaction", "interaction_se", "interaction_p_raw") := raw[, .(estimate, se, p)]]
    if (type == "binary") {
      inside <- scan_subset(ctx_frame[ctx == 1], "exposed", base_covariates, mat_raw, trend = TRUE)
      outside <- scan_subset(ctx_frame[ctx == 0], "exposed", base_covariates, mat_raw, trend = TRUE)
      if (!is.null(inside)) out[, c("inside", "inside_se", "inside_p_raw") := inside[, .(estimate, se, p)]]
      if (!is.null(outside)) out[, c("outside", "outside_se", "outside_p_raw") := outside[, .(estimate, se, p)]]
    }
  }
  out[, context_id := id]
  out[, c("n_neg_in", "n_pos_in", "n_neg_out", "n_pos_out") := as.list(cells)]
  out
}
run_scan <- function(shuffle) {
  out <- rbindlist(parallel::mclapply(seq_len(nrow(definitions)), fit_context, shuffle = shuffle, mc.cores = n_cores), fill = TRUE)
  for (col in c("inside_rank", "inside_p")) if (!col %in% names(out)) out[, (col) := NA_real_]
  out <- merge(out, definitions[, .(context_id, class, label, type)], by = "context_id", sort = FALSE)
  merge(out, gene_ann[, .(feature, symbol, chromosome, y_partner, selective)], by = "feature", sort = FALSE)
}
call_hits <- function(dt) {
  dt <- add_family_fdr(dt, "interaction_p", "interaction")
  dt <- add_family_fdr(dt, "inside_p", "inside")
  # A context-specific claim needs the interaction and, for a binary context, the effect inside it.
  dt[, context_specific := interaction_fdr_class < fdr_cut & (type == "continuous" | inside_fdr_class < fdr_cut)]
  dt[]
}
cell_bin <- function(dt) cut(pmin(dt$n_neg_in, dt$n_pos_in, dt$n_neg_out, dt$n_pos_out, na.rm = TRUE), c(4, 6, 9, 14, 24, 49, Inf), labels = c("5-6", "7-9", "10-14", "15-24", "25-49", "50+"))

results <- run_scan(0)
for (col in c("interaction", "interaction_se", "interaction_p_raw", "inside", "inside_se", "inside_p_raw", "outside", "outside_se", "outside_p_raw")) {
  if (!col %in% names(results)) results[, (col) := NA_real_]
}
results[, focused := feature %in% focus]
# Smallest raw-scale interaction detectable with 80% power at a nominal two-sided 0.05.
results[, interaction_mde := power_multiplier * abs(interaction_se)]

# Families: the genome-wide screen holds selective genes only; the focused family holds hits and X-Y partners.
screen <- call_hits(results[selective == TRUE])
screen[, detected_inside_only := !(context_specific %in% TRUE) & type == "binary" & inside_fdr_class < fdr_cut]
focused <- results[focused == TRUE]
focused[, interaction_fdr_focused := bh(interaction_p)]
focused[, inside_fdr_focused := bh(inside_p)]
focused[, context_specific := interaction_fdr_focused < fdr_cut & (type == "continuous" | inside_fdr_focused < fdr_cut)]

# --- the same procedure under shuffled Y status ------------------------------
tail_table <- function(dt, run) {
  dt[, .(
    run = run, tests = .N, below_0.05 = mean(interaction_p < 0.05, na.rm = TRUE), below_0.001 = mean(interaction_p < 0.001, na.rm = TRUE),
    interaction_fdr_class_hits = sum(interaction_fdr_class < fdr_cut, na.rm = TRUE), context_specific = sum(context_specific, na.rm = TRUE)
  ), keyby = .(smallest_cell = cell_bin(dt))]
}
calibration <- rbind(
  tail_table(screen, "observed"),
  rbindlist(lapply(seq_len(n_null), function(b) tail_table(call_hits(run_scan(b)[selective == TRUE]), paste0("shuffled_", b))))
)
write_result(calibration, "context_gene_calibration")
null_summary <- calibration[, .(interaction_fdr_class_hits = sum(interaction_fdr_class_hits), context_specific = sum(context_specific), below_0.001 = weighted.mean(below_0.001, tests)), by = run]
print(null_summary)

saveRDS(screen, file.path(derived_dir, "context_gene_screen.rds"))
fwrite(screen, file.path(derived_dir, "context_gene_screen.csv.gz"))
write_result(focused[order(interaction_p)], "context_gene_focused")
hits <- rbind(
  screen[context_specific == TRUE][, family := "screen"],
  focused[context_specific == TRUE][, family := "focused"],
  fill = TRUE
)
write_result(hits[order(interaction_p)], "context_gene_hits")
write_result(screen[detected_inside_only == TRUE][order(inside_p)], "context_gene_inside_only")

summary_tab <- screen[, .(
  tests = .N, min_interaction_p = min(interaction_p, na.rm = TRUE),
  interaction_fdr_class_hits = sum(interaction_fdr_class < fdr_cut, na.rm = TRUE),
  interaction_fdr_global_hits = sum(interaction_fdr_global < fdr_cut, na.rm = TRUE),
  inside_fdr_class_hits = sum(inside_fdr_class < fdr_cut, na.rm = TRUE),
  context_specific = sum(context_specific, na.rm = TRUE),
  median_interaction_mde = median(interaction_mde, na.rm = TRUE)
), by = .(class, context_id, label, type, n_neg_in, n_pos_in, n_neg_out, n_pos_out)]
write_result(summary_tab[order(class, -context_specific, min_interaction_p)], "context_gene_summary")
print(summary_tab[, .(contexts = .N, tests = sum(tests), interaction_hits = sum(interaction_fdr_class_hits), inside_hits = sum(inside_fdr_class_hits), context_specific = sum(context_specific), median_mde = median(median_interaction_mde)), by = class])

# Heterogeneity of the Y effect across lineages (Cochran's Q on within-lineage raw-scale estimates).
lineage_rows <- results[class == "lineage" & is.finite(inside) & is.finite(inside_se)]
heterogeneity <- lineage_rows[, {
  q <- cochran_q(inside, abs(inside_se))
  .(lineages = .N, pooled = q$pooled, q = q$q, q_df = q$df, q_p = q$p, i2 = q$i2)
}, by = .(feature, symbol, focused, selective)]
heterogeneity[selective == TRUE, q_fdr := bh(q_p)]
heterogeneity[focused == TRUE, q_fdr_focused := bh(q_p)]
write_result(heterogeneity[order(q_p)], "context_gene_lineage_heterogeneity")
per_lineage <- lineage_rows[focused == TRUE, .(feature, symbol, lineage = label, estimate = inside, se = abs(inside_se), p = inside_p, n_neg_in, n_pos_in)]
write_result(per_lineage, "context_gene_focused_by_lineage")

null_hits <- null_summary[run != "observed"]
log_msg("context scan done: %d contexts fitted, %d interaction tests in the screen, %d context-specific (screen), %d (focused); shuffled runs gave %s",
  uniqueN(screen$context_id), nrow(screen), sum(screen$context_specific, na.rm = TRUE), sum(focused$context_specific, na.rm = TRUE), paste(null_hits$context_specific, collapse = ", "))
print(hits[order(interaction_p)][1:min(.N, 25), .(family, symbol, label, class, interaction, interaction_p, inside, inside_p, outside, n_neg_in, n_pos_in)])
jsonlite::write_json(list(
  lines = nrow(frame), genes_screen = sum(gene_ann$selective), genes_focused = length(focus),
  contexts_fitted = uniqueN(screen$context_id), screen_tests = nrow(screen), focused_tests = nrow(focused),
  screen_context_specific = sum(screen$context_specific, na.rm = TRUE), focused_context_specific = sum(focused$context_specific, na.rm = TRUE),
  screen_interaction_fdr_class_hits = sum(screen$interaction_fdr_class < fdr_cut, na.rm = TRUE),
  screen_interaction_fdr_global_hits = sum(screen$interaction_fdr_global < fdr_cut, na.rm = TRUE),
  inside_only = sum(screen$detected_inside_only, na.rm = TRUE),
  heterogeneity_fdr_hits = sum(heterogeneity$q_fdr < fdr_cut, na.rm = TRUE),
  shuffled_interaction_fdr_class_hits = null_hits$interaction_fdr_class_hits, shuffled_context_specific = null_hits$context_specific
), file.path(derived_dir, "context_gene_summary.json"), auto_unbox = TRUE, pretty = TRUE)
