# Stage 28: combination of the CRISPR sources and the table across assays.
# 1. Pan-cancer: cell lines (stored scan) and organoids, combined per gene from each
#    source's own p-value and sign, weighted by the design (fixed-effect pooling of
#    standardised effects). Scores are on different scales, so effect sizes are
#    never pooled.
# 2. Within cancer type: the rank-based test of the Y effect inside the matching
#    cell-line context and organoid type. Colorectal carries MSI as a covariate on
#    both sides, so the cell-line side is refitted here after the stored p-values
#    have been reproduced without it.
# 3. CRISPR against RNAi, context by context: side by side, no combined p-value
#    (most lines are shared).
# 4. The X partners across sources, and what each test could detect at best.
# A combined hit is "detected across sources"; "specific" stays with the
# interaction rule inside one source. Sources enter as columns, so a third
# CRISPR source can be added.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("organoid")

cp <- ext_params$combine
sp <- ext_params$shared
fdr_cut <- sp$fdr
base_covariates <- c("LibraryKY", "LibraryHumagneCD")
cohort <- read_result("cohort")
context <- read_result("context_matrix")
gene_ann <- read_result("gene_annotation")
effect <- load_matrix("crispr_effect")
rnai <- load_matrix("rnai")
pan_hits <- read_result("gene_pan_hits")
cell_pan <- read_result("gene_pan_results")
rnai_pan <- read_result("gene_pan_rnai")
screen <- readRDS(file.path(derived_dir, "context_gene_screen.rds"))
setDT(screen)
organoid_cohort <- read_ext("organoid_cohort")[in_primary == TRUE]
organoid_genes <- read_ext("organoid_gene_annotation")
organoid_pan <- read_ext("organoid_pan_results")
organoid_types <- read_ext("organoid_type_screen")
unit_of <- function(design, coef) sqrt(solve(crossprod(design))[coef, coef])
replicated <- function(z_other) !is.na(z_other) & 2 * pnorm(-abs(z_other)) < sp$replication_p

# --- 1. pan-cancer -----------------------------------------------------------------------------
frame <- analysis_frame(cohort, effect)
pan_design <- build_design(frame, tested = "exposed", covariates = base_covariates)
stopifnot(pan_design$estimable, length(pan_design$rows) == nrow(frame))
cell_unit <- unit_of(pan_design$design, "exposed")
organoid_unit <- unique(round(organoid_pan$unit, 12))
stopifnot(length(organoid_unit) == 1)
shared <- merge(
  cell_pan[, .(feature, symbol, chromosome, y_partner, cell_estimate = estimate, cell_p = p, cell_fdr = fdr, cell_hit = associated)],
  organoid_pan[!is.na(feature), .(feature, organoid_symbol = symbol, organoid_estimate = estimate, organoid_p = p, organoid_fdr = fdr, organoid_hit = associated, seen_before)],
  by = "feature"
)
assert_unique(shared$feature, "genes shared by the two scans")
z <- cbind(cell_lines = signed_z(shared$cell_p, shared$cell_estimate), organoids = signed_z(shared$organoid_p, shared$organoid_estimate))
units <- c(cell_lines = cell_unit, organoids = organoid_unit)
pooled <- stouffer(z, 1 / units)
shared[, c("cell_z", "organoid_z") := .(z[, "cell_lines"], z[, "organoids"])]
shared[, c("cell_smd", "organoid_smd") := .(cell_z * cell_unit, organoid_z * organoid_unit)]
shared[, c("combined_z", "combined_p") := .(pooled$z, pooled$p)]
shared[, combined_fdr := bh(combined_p)]
shared[, combined_hit := combined_fdr < fdr_cut]
difference <- source_difference(cbind(shared$cell_smd, shared$organoid_smd), units)
shared[, c("smd_difference", "smd_difference_low", "smd_difference_high", "smd_difference_p") := difference[, .(difference, ci_low, ci_high, p)]]
shared[, evidence := NA_character_]
if (any(shared$combined_hit)) shared[combined_hit == TRUE, evidence := evidence_class(z[shared$combined_hit, , drop = FALSE], 1 / units, alpha = fdr_cut)]
shared[, cell_hit_replicated := cell_hit & sign(organoid_z) == sign(cell_z) & replicated(organoid_z)]
shared[, organoid_hit_replicated := organoid_hit & sign(organoid_z) == sign(cell_z) & replicated(cell_z)]
setorder(shared, combined_p)
write_ext(shared, "combined_pan")
log_msg("pan-cancer: %d genes in both scans; weights %.1f (cell lines, unit %.4f) and %.1f (organoids, unit %.4f), organoids carry %.0f%% of the weight; %d combined at FDR < %.2f",
  nrow(shared), 1 / cell_unit, cell_unit, 1 / organoid_unit, organoid_unit, 100 * (1 / organoid_unit^2) / sum(1 / units^2), sum(shared$combined_hit), fdr_cut)
print(shared[combined_hit == TRUE | cell_hit == TRUE | organoid_hit == TRUE, .(symbol, cell_estimate, cell_p, organoid_estimate, organoid_p, combined_p, combined_fdr, evidence, smd_difference, smd_difference_p, cell_hit, organoid_hit)])
# Genes hit in one scan and absent from the other.
only_cell <- cell_pan[associated == TRUE & !feature %in% shared$feature, symbol]
log_msg("cell-line hits not in the organoid gene set: %s", if (length(only_cell)) paste(only_cell, collapse = ", ") else "none")

# The three genes seen in the organoid comparison on disk, in every source.
seen <- rbindlist(lapply(ext_params$organoid$seen_before, function(gene) {
  data.table(
    symbol = gene,
    organoid_estimate = organoid_pan[symbol == gene, estimate], organoid_p = organoid_pan[symbol == gene, p], organoid_fdr = organoid_pan[symbol == gene, fdr],
    cell_estimate = cell_pan[symbol == gene, estimate][1], cell_p = cell_pan[symbol == gene, p][1], cell_fdr = cell_pan[symbol == gene, fdr][1],
    rnai_estimate = rnai_pan[symbol == gene, estimate][1], rnai_p = rnai_pan[symbol == gene, p][1], rnai_fdr = rnai_pan[symbol == gene, fdr][1]
  )
}))
write_ext(seen, "combined_seen_before")
print(seen)

# --- 2. within cancer type ------------------------------------------------------------------------
# The screen's rows, as in stage 12: selective genes plus the focus family, ranked over all lines.
focus <- union(pan_hits$feature, gene_ann[in_scan == TRUE & !is.na(y_partner), feature])
universe <- union(gene_ann[selective == TRUE, feature], focus)
mat_rank <- t(apply(t(effect[frame$ModelID, universe]), 1, rank_normalise))
dimnames(mat_rank) <- list(universe, frame$ModelID)
frame[, MSI := as.numeric(cohort$MSIScore[match(ModelID, cohort$ModelID)] >= params$context$msi_score_min)]
pairs <- as.data.table(cp$pairs)
reproduce_gap <- c()
cell_inside_units <- c()
within_type <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  id <- pairs$cell_line_context[i]
  ctx_frame <- with_context(frame, context[[id]][match(frame$ModelID, context$ModelID)], "binary")
  stored <- screen[context_id == id]
  baseline <- suppressWarnings(fit_context_tests(ctx_frame, mat_rank, base_covariates))
  at <- match(stored$feature, baseline$feature)
  gap <- max(abs(log10(baseline$interaction_p[at]) - log10(stored$interaction_p)), abs(log10(baseline$inside_p[at]) - log10(stored$inside_p)), na.rm = TRUE)
  if (!is.finite(gap) || gap > cp$reproduce_log10_max) stop(sprintf("%s: stored p-values not reproduced (gap %.2g)", id, gap))
  reproduce_gap[[id]] <<- gap
  covariates <- if (pairs$cell_line_msi_covariate[i]) c(base_covariates, "MSI") else base_covariates
  inside_lines <- ctx_frame[ctx == 1]
  if (pairs$cell_line_msi_covariate[i]) stopifnot(!anyNA(inside_lines$MSI))
  fitted <- if (pairs$cell_line_msi_covariate[i]) suppressWarnings(fit_context_tests(ctx_frame[!is.na(MSI)], mat_rank, covariates)) else baseline
  inside_design <- build_design(inside_lines, tested = "exposed", covariates = covariates)
  unit_cell <- unit_of(inside_design$design, "exposed")
  cell_inside_units[[id]] <<- unit_cell
  cell <- merge(fitted[, .(feature, cell_inside_rank = inside_estimate, cell_inside_p = inside_p, cell_interaction_p = interaction_p)],
    baseline[, .(feature, cell_inside_p_stored_model = inside_p)], by = "feature", sort = FALSE)
  cell <- cell[feature %in% gene_ann[selective == TRUE, feature]]
  organoid <- organoid_types[type == pairs$organoid_type[i] & !is.na(feature) & !is.na(inside_p),
    .(feature, organoid_symbol = symbol, organoid_inside_rank = inside_rank, organoid_inside_p = inside_p, organoid_unit = inside_unit, organoid_interaction_p = interaction_p, organoid_selective = selective, seen_before)]
  both <- merge(cell, organoid, by = "feature")
  both[, cell_unit := unit_cell]
  z <- cbind(signed_z(both$cell_inside_p, both$cell_inside_rank), signed_z(both$organoid_inside_p, both$organoid_inside_rank))
  weights <- cbind(1 / both$cell_unit, 1 / both$organoid_unit)
  pooled <- stouffer(z, weights)
  both[, c("cell_z", "organoid_z", "combined_z", "combined_p") := .(z[, 1], z[, 2], pooled$z, pooled$p)]
  both[, combined_fdr := bh(combined_p)]
  both[, combined_hit := combined_fdr < fdr_cut]
  difference <- source_difference(cbind(both$cell_z * both$cell_unit, both$organoid_z * both$organoid_unit), cbind(both$cell_unit, both$organoid_unit))
  both[, c("smd_difference", "smd_difference_low", "smd_difference_high", "smd_difference_p") := difference[, .(difference, ci_low, ci_high, p)]]
  both[, evidence := NA_character_]
  if (any(both$combined_hit)) both[combined_hit == TRUE, evidence := evidence_class(z[both$combined_hit, , drop = FALSE], weights[both$combined_hit, , drop = FALSE], alpha = fdr_cut)]
  both[, pair := pairs$label[i]]
  both[, c("cell_n_negative", "cell_n_positive") := .(sum(inside_lines$exposed == 1), sum(inside_lines$exposed == 0))]
  type_cohort <- organoid_cohort[OncotreeLineage == pairs$organoid_type[i]]
  both[, c("organoid_n_negative", "organoid_n_positive") := .(sum(type_cohort$y_negative == 1), sum(type_cohort$y_negative == 0))]
  both[, msi_adjusted_cell_lines := pairs$cell_line_msi_covariate[i]]
  both[, symbol := gene_symbol(feature)]
  both
}))
setorder(within_type, combined_p)
write_ext(within_type, "combined_within_type")
print(within_type[, .(genes = .N, combined_hits = sum(combined_hit), min_combined_p = min(combined_p), cell_weight_share = round(mean((1 / cell_unit^2) / (1 / cell_unit^2 + 1 / organoid_unit^2)), 2),
  cell = sprintf("%d vs %d", cell_n_negative[1], cell_n_positive[1]), organoid = sprintf("%d vs %d", organoid_n_negative[1], organoid_n_positive[1])), by = pair])
print(within_type[, .SD[1:8], by = pair][, .(pair, symbol, cell_inside_rank, cell_inside_p, organoid_inside_rank, organoid_inside_p, combined_p, combined_fdr, evidence)])
log_msg("stored cell-line p-values reproduced before the refit: largest gap %.2g log10 units", max(unlist(reproduce_gap)))

# --- 3. CRISPR against RNAi, context by context -------------------------------------------------------
rnai_screen <- rbind(
  readRDS(file.path(ext_out_dir, "rnai_context_screen.rds"))[, .(symbol = feature, context_id, rnai_interaction_rank = interaction_rank, rnai_interaction_p = interaction_p, rnai_inside_rank = inside_rank, rnai_inside_p = inside_p, rnai_inside = inside, rnai_lines = n)],
  read_ext("rnai_context_focused")[, .(symbol = feature, context_id, rnai_interaction_rank = interaction_rank, rnai_interaction_p = interaction_p, rnai_inside_rank = inside_rank, rnai_inside_p = inside_p, rnai_inside = inside, rnai_lines = n)]
)
rnai_screen <- unique(rnai_screen, by = c("symbol", "context_id"))
crispr_screen <- rbind(
  screen[, .(symbol, context_id, class, label, crispr_interaction_rank = interaction_rank, crispr_interaction_p = interaction_p, crispr_inside_rank = inside_rank, crispr_inside_p = inside_p, crispr_inside = inside)],
  read_result("context_gene_focused")[, .(symbol, context_id, class, label, crispr_interaction_rank = interaction_rank, crispr_interaction_p = interaction_p, crispr_inside_rank = inside_rank, crispr_inside_p = inside_p, crispr_inside = inside)]
)
crispr_screen <- unique(crispr_screen, by = c("symbol", "context_id"))
cross <- merge(crispr_screen, rnai_screen, by = c("symbol", "context_id"))
cross[, inside_both_nominal := crispr_inside_p < 0.05 & rnai_inside_p < 0.05 & sign(crispr_inside_rank) == sign(rnai_inside_rank)]
cross[, interaction_both_nominal := crispr_interaction_p < 0.05 & rnai_interaction_p < 0.05 & sign(crispr_interaction_rank) == sign(rnai_interaction_rank)]
cross[, cancer_type := class %in% c("lineage", "disease", "subtype")]
setorder(cross, crispr_inside_p)
saveRDS(cross, file.path(ext_out_dir, "combined_cross_assay.rds"))
fwrite(cross, file.path(ext_out_dir, "combined_cross_assay.csv.gz"))
# Agreement expected if the assays were unrelated: each nominal at 0.05 with matching sign, 0.05 * 0.05 / 2.
cross_summary <- cross[, .(
  pairs = .N, inside_tested = sum(!is.na(crispr_inside_p) & !is.na(rnai_inside_p)), inside_both_nominal = sum(inside_both_nominal, na.rm = TRUE),
  inside_expected_by_chance = sum(!is.na(crispr_inside_p) & !is.na(rnai_inside_p)) * 0.05 * 0.05 / 2,
  interaction_both_nominal = sum(interaction_both_nominal, na.rm = TRUE), interaction_expected_by_chance = .N * 0.05 * 0.05 / 2,
  inside_rank_correlation = cor(crispr_inside_rank, rnai_inside_rank, use = "complete.obs", method = "spearman")
), by = .(cancer_type)]
write_ext(cross_summary, "combined_cross_assay_summary")
print(cross_summary)
# The cell-line leads of the context screen, looked up in RNAi.
leads <- rbind(read_result("context_gene_hits")[, .(symbol, context_id)], read_result("context_gene_inside_only")[, .(symbol, context_id)], read_result("power_alone_hits")[, .(symbol, context_id)], fill = TRUE)
leads <- unique(leads[!is.na(symbol)])
lead_lookup <- merge(leads, cross, by = c("symbol", "context_id"), all.x = TRUE)
lead_lookup[, tested_in_rnai := !is.na(rnai_lines)]
write_ext(lead_lookup, "combined_cell_line_leads_in_rnai")
print(lead_lookup[, .(symbol, context_id, crispr_inside, crispr_inside_p, rnai_inside, rnai_inside_p, tested_in_rnai)])

# --- 4. the X partners across sources ------------------------------------------------------------------
rnai_frame <- analysis_frame(cohort, rnai)
partner_rows <- rbindlist(lapply(params$xy_paralogs$x_partner, function(gene) {
  rows <- list()
  cell <- cell_pan[symbol == gene]
  if (nrow(cell) == 1) rows$cell <- data.table(source = "Cell lines, CRISPR", z = signed_z(cell$p, cell$estimate), unit = cell_unit, p = cell$p, n_negative = cell$n_negative, n_positive = cell$n_positive)
  if (gene %in% colnames(rnai)) {
    lines <- rnai_frame[is.finite(rnai[ModelID, gene])]
    design <- build_design(lines, tested = "exposed")
    row <- rnai_pan[symbol == gene]
    if (nrow(row) == 1) rows$rnai <- data.table(source = "Cell lines, RNAi", z = signed_z(row$p, row$estimate), unit = unit_of(design$design, "exposed"), p = row$p, n_negative = sum(lines$exposed == 1), n_positive = sum(lines$exposed == 0))
  }
  organoid <- organoid_pan[symbol == gene]
  if (nrow(organoid) == 1) rows$organoid <- data.table(source = "Organoids, CRISPR", z = signed_z(organoid$p, organoid$estimate), unit = organoid$unit, p = organoid$p, n_negative = organoid$n_negative, n_positive = organoid$n_positive)
  both <- shared[symbol == gene]
  if (nrow(both) == 1) rows$combined <- data.table(source = "CRISPR combined", z = both$combined_z, unit = 1 / sqrt(sum(1 / units^2)), p = both$combined_p, n_negative = NA_integer_, n_positive = NA_integer_)
  if (length(rows) == 0) return(NULL)
  cbind(data.table(symbol = gene, y_partner = params$xy_paralogs$y_gene[match(gene, params$xy_paralogs$x_partner)]), rbindlist(rows))
}))
partner_rows[, smd := z * unit]
partner_rows[, c("smd_low", "smd_high") := .(smd - qnorm(0.975) * unit, smd + qnorm(0.975) * unit)]
write_ext(partner_rows, "combined_xy_partners")
print(dcast(partner_rows, symbol ~ source, value.var = "smd"))

# --- 5. what each test could detect at best ----------------------------------------------------------------
# Standardised effect found with the stated power, for a named gene (p < 0.05) and for
# one more discovery in the family as it stands; normal theory on each design.
power_row <- function(analysis, source, unit, df, family_p, n_negative, n_positive) {
  critical <- bh_next_critical(family_p, fdr_cut)
  data.table(
    analysis = analysis, source = source, n_negative = n_negative, n_positive = n_positive, unit = unit, df = df, family = sum(is.finite(family_p)),
    named = analytic_mde(unit, df, 0.05, sp$power), family_critical_p = critical, genome_wide = analytic_mde(unit, df, critical, sp$power)
  )
}
pooled_unit <- function(u) 1 / sqrt(sum(1 / u^2))
# Residual degrees of freedom of each organoid design (the moderated tests add prior degrees of freedom on top).
organoid_frame <- copy(organoid_cohort)[, exposed := y_negative]
organoid_covariates <- c("LibraryV1", "MSI")
organoid_df <- function(rows, lineage) nrow(rows) - ncol(build_design(rows, tested = "exposed", covariates = organoid_covariates, lineage = lineage)$design)
organoid_pan_df <- organoid_df(organoid_frame, TRUE)
pan_rows <- rbind(
  power_row("Pan-cancer", "Cell lines", cell_unit, nrow(frame) - ncol(pan_design$design), cell_pan$p, sum(frame$exposed == 1), sum(frame$exposed == 0)),
  power_row("Pan-cancer", "Organoids", organoid_unit, organoid_pan_df, organoid_pan$p, sum(organoid_frame$y_negative == 1), sum(organoid_frame$y_negative == 0)),
  power_row("Pan-cancer", "Combined", pooled_unit(units), nrow(frame) - ncol(pan_design$design) + organoid_pan_df, shared$combined_p, sum(frame$exposed == 1) + sum(organoid_frame$y_negative == 1), sum(frame$exposed == 0) + sum(organoid_frame$y_negative == 0))
)
type_rows <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  rows <- within_type[pair == pairs$label[i]]
  id <- pairs$cell_line_context[i]
  unit_cell <- rows$cell_unit[1]
  unit_organoid <- median(rows$organoid_unit)
  df_cell <- rows$cell_n_negative[1] + rows$cell_n_positive[1] - 2 - length(base_covariates) - as.numeric(pairs$cell_line_msi_covariate[i])
  df_organoid <- organoid_df(organoid_frame[OncotreeLineage == pairs$organoid_type[i]], FALSE)
  rbind(
    power_row(pairs$label[i], "Cell lines", unit_cell, df_cell, screen[context_id == id, inside_p], rows$cell_n_negative[1], rows$cell_n_positive[1]),
    power_row(pairs$label[i], "Organoids", unit_organoid, df_organoid, organoid_types[selective == TRUE, inside_p], rows$organoid_n_negative[1], rows$organoid_n_positive[1]),
    power_row(pairs$label[i], "Combined", pooled_unit(c(unit_cell, unit_organoid)), df_cell + df_organoid, rows$combined_p, rows$cell_n_negative[1] + rows$organoid_n_negative[1], rows$cell_n_positive[1] + rows$organoid_n_positive[1])
  )
}))
rnai_design <- build_design(rnai_frame, tested = "exposed")
rnai_row <- power_row("Pan-cancer", "RNAi (lines with every library)", unit_of(rnai_design$design, "exposed"), nrow(rnai_frame) - ncol(rnai_design$design), rnai_pan$p, sum(rnai_frame$exposed == 1), sum(rnai_frame$exposed == 0))
power <- rbind(pan_rows, rnai_row, type_rows)
write_ext(power, "extension_power")
print(power[, .(analysis, source, n_negative, n_positive, unit = round(unit, 3), family, named = round(named, 2), genome_wide = round(genome_wide, 2))])

jsonlite::write_json(list(
  pan_genes = nrow(shared), pan_combined_hits = sum(shared$combined_hit), pan_hit_genes = shared[combined_hit == TRUE, symbol], pan_evidence = as.list(table(shared[combined_hit == TRUE, evidence])),
  organoid_weight_share = (1 / organoid_unit^2) / sum(1 / units^2), cell_unit = cell_unit, organoid_unit = organoid_unit,
  cell_hits_not_in_organoids = only_cell,
  within_type = lapply(split(within_type, by = "pair"), function(d) list(genes = nrow(d), combined_hits = sum(d$combined_hit), hit_genes = d[combined_hit == TRUE, symbol], min_combined_p = min(d$combined_p))),
  reproduce_gap_log10 = as.list(unlist(reproduce_gap)), cross_assay_pairs = nrow(cross),
  parameters_sha256 = ext_section_hash(c("shared", "combine"))
), file.path(ext_out_dir, "combined_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("combination done: %d pan-cancer hits; within type %s", sum(shared$combined_hit), paste(sprintf("%s %d", pairs$label, vapply(pairs$label, function(l) sum(within_type[pair == l, combined_hit]), numeric(1))), collapse = "; "))
