# Stage 33: the three CRISPR sources together (cell lines, Sanger organoids, NextGen).
# 1. All cancer types: each source's stored p-value and sign, combined per gene with
#    weights from each design. A gene missing from one source is combined over the
#    others. The two-source result of stage 28 stays as it is and is carried beside.
# 2. Cell lines against NextGen on one scale: both are Chronos scores, so the two
#    effect sizes are compared directly.
# 3. Inside a cancer type, source by source, where NextGen meets the gate.
# 4. The pooled test for cancer types that no source can test alone (listed by
#    stage 31): ranks within source, Y groups compared within strata (library in
#    the cell lines, screen format in NextGen), p-values from every arrangement
#    of Y status within strata. The four paired genes are a named family; for the
#    rest the table is a ranking and no gene is called.
# 5. The same test for the paired genes in three more lineages (secondary, no calls).
# 6. Donors screened at both institutes, the X partners across sources, the
#    follow-up of every combined hit, and what each test could detect at best.
# LOY_EXTENSION_TRIAL=1 reads a trial run of stage 32 and shuffles Y status
# within strata in the pooled tests (use with LOY_EXTENSION_OUTPUT).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("organoid", dir = ext_dir)
assert_extension_freeze("nextgen")

np <- ext_params$nextgen
sp <- ext_params$shared
fdr_cut <- sp$fdr
trial <- Sys.getenv("LOY_EXTENSION_TRIAL") == "1"
scan_summary <- jsonlite::fromJSON(file.path(ext_out_dir, "nextgen_scan_summary.json"))
if (!identical(isTRUE(scan_summary$trial), trial)) stop("the NextGen scan in the output folder and this run differ in trial mode")
unit_of <- function(design, coef) sqrt(solve(crossprod(design))[coef, coef])
pooled_unit <- function(u) 1 / sqrt(sum(1 / u^2))
replicated <- function(z_other) !is.na(z_other) & 2 * pnorm(-abs(z_other)) < sp$replication_p
interval <- qnorm(0.975)

cell_cohort <- read_result("cohort")
cell_effect <- load_matrix("crispr_effect")
cell_genes <- read_result("gene_annotation")
cell_pan <- read_result("gene_pan_results")
screen <- readRDS(file.path(derived_dir, "context_gene_screen.rds"))
setDT(screen)
organoid_cohort <- read_ext("organoid_cohort", dir = ext_dir)
organoid_primary <- organoid_cohort[in_primary == TRUE]
organoid_genes <- read_ext("organoid_gene_annotation", dir = ext_dir)
organoid_pan <- read_ext("organoid_pan_results", dir = ext_dir)
organoid_types <- read_ext("organoid_type_screen", dir = ext_dir)
organoid_effect <- readRDS(file.path(ext_dir, "organoid_effect.rds"))
two_pan <- read_ext("combined_pan", dir = ext_dir)
two_type <- read_ext("combined_within_type", dir = ext_dir)
two_partners <- read_ext("combined_xy_partners", dir = ext_dir)
two_power <- read_ext("extension_power", dir = ext_dir)
ng_cohort <- read_ext("nextgen_cohort")
ng_primary <- ng_cohort[in_primary == TRUE][order(ModelID)]
ng_genes <- read_ext("nextgen_gene_annotation")
ng_pan <- read_ext("nextgen_pan_results")
ng_types <- read_ext("nextgen_type_screen")
ng_effect <- readRDS(file.path(ext_out_dir, "nextgen_effect.rds"))
type_counts <- read_ext("nextgen_pooled_types")
test_counts <- read_ext("nextgen_test_counts")
named_features <- ng_genes[in_scan == TRUE][match(np$named_genes, symbol), feature]

# --- 1. all cancer types -------------------------------------------------------------------------
cell_frame <- analysis_frame(cell_cohort, cell_effect)
cell_covariates <- c("LibraryKY", "LibraryHumagneCD")
pan_design <- build_design(cell_frame, tested = "exposed", covariates = cell_covariates)
stopifnot(pan_design$estimable, length(pan_design$rows) == nrow(cell_frame))
cell_unit <- unit_of(pan_design$design, "exposed")
# The sources share no model, and no donor as far as the release's model table knows the NextGen models.
model_table <- fread(input_path("model"), select = c("ModelID", "PatientID"))
ng_patients <- model_table[ModelID %in% ng_primary$ModelID, PatientID]
stopifnot(!any(ng_cohort$ModelID %in% cell_frame$ModelID), !any(ng_patients %in% cell_frame$PatientID), !any(ng_primary$hcmi_donor %in% hcmi_donor(organoid_cohort$ModelID)[!is.na(hcmi_donor(organoid_cohort$ModelID))]))
log_msg("independence: no NextGen model among the %d cell lines; %d NextGen cohort models are in the 26Q1 model table and none shares a donor with a cell line; none shares a donor with a screened Sanger organoid", nrow(cell_frame), length(ng_patients))
organoid_unit <- unique(round(organoid_pan$unit, 12))
ng_unit <- unique(round(ng_pan$unit, 12))
stopifnot(length(organoid_unit) == 1, length(ng_unit) == 1)
units <- c(cell_lines = cell_unit, organoids = organoid_unit, nextgen = ng_unit)
three <- merge(
  cell_pan[, .(feature, symbol, chromosome, y_partner, cell_estimate = estimate, cell_se = (ci_high - ci_low) / (2 * qt(0.975, df)), cell_p = p, cell_fdr = fdr, cell_hit = associated)],
  organoid_pan[!is.na(feature), .(feature, organoid_estimate = estimate, organoid_p = p, organoid_fdr = fdr, organoid_hit = associated)],
  by = "feature", all.x = TRUE
)
three <- merge(three, ng_pan[, .(feature, nextgen_estimate = estimate, nextgen_se = se, nextgen_p = p, nextgen_fdr = fdr, nextgen_hit = associated, nextgen_discovery = discovery)], by = "feature", all.x = TRUE)
# The family: genes in the cell-line scan and in at least one other scan.
three <- three[!is.na(organoid_p) | !is.na(nextgen_p)]
assert_unique(three$feature, "genes of the three-source combination")
if (!trial) stopifnot(nrow(three) == test_counts$combination_pan_genes[1])
z <- cbind(cell_lines = signed_z(three$cell_p, three$cell_estimate), organoids = signed_z(three$organoid_p, three$organoid_estimate), nextgen = signed_z(three$nextgen_p, three$nextgen_estimate))
# Same code path as stage 28: the cell-line and organoid statistics equal the stored ones.
at <- match(two_pan$feature, three$feature)
stopifnot(
  !anyNA(at), max(abs(z[at, "cell_lines"] - two_pan$cell_z)) < 1e-12, max(abs(z[at, "organoids"] - two_pan$organoid_z)) < 1e-12,
  abs(cell_unit - two_power[analysis == "Pan-cancer" & source == "Cell lines", unit]) < 1e-12, abs(organoid_unit - two_power[analysis == "Pan-cancer" & source == "Organoids", unit]) < 1e-10
)
pooled <- stouffer(z, 1 / units)
three[, c("cell_z", "organoid_z", "nextgen_z") := .(z[, "cell_lines"], z[, "organoids"], z[, "nextgen"])]
three[, c("combined_z", "combined_p", "n_sources") := .(pooled$z, pooled$p, pooled$n_sources)]
three[, combined_fdr := bh(combined_p)]
three[, combined_hit := combined_fdr < fdr_cut]
# The p-value a gene at this rank needs to pass BH in this family.
three[, bh_cut := frank(combined_p, ties.method = "max") * fdr_cut / .N]
three[, c("n_supporting", "evidence") := .(NA_integer_, NA_character_)]
if (any(three$combined_hit)) {
  support <- evidence_summary(z[three$combined_hit, , drop = FALSE], 1 / units, alpha = fdr_cut)
  three[combined_hit == TRUE, c("n_supporting", "evidence") := .(support$n_supporting, support$label)]
}
at <- match(three$feature, two_pan$feature)
three[, c("two_source_p", "two_source_fdr", "two_source_hit", "two_source_evidence") := two_pan[at, .(combined_p, combined_fdr, combined_hit, evidence)]]
three[, cell_hit_in_nextgen := cell_hit & sign(nextgen_z) == sign(cell_z) & replicated(nextgen_z)]
three[, nextgen_hit_in_cell_lines := nextgen_hit & sign(nextgen_z) == sign(cell_z) & replicated(cell_z)]
# --- 2. cell lines against NextGen on the Chronos scale --------------------------------------------
three[, chronos_difference := nextgen_estimate - cell_estimate]
three[, chronos_difference_se := sqrt(nextgen_se^2 + cell_se^2)]
three[, c("chronos_difference_low", "chronos_difference_high") := .(chronos_difference - interval * chronos_difference_se, chronos_difference + interval * chronos_difference_se)]
three[, chronos_difference_p := 2 * pnorm(-abs(chronos_difference / chronos_difference_se))]
setorder(three, combined_p)
write_ext(three, "nextgen_combined_pan")
share <- (1 / units^2) / sum(1 / units^2)
log_msg("all cancer types: %d genes in the cell-line scan and at least one other (%d in all three, %d without NextGen, %d without the organoids); weight shares %.0f%% cell lines, %.0f%% organoids, %.0f%% NextGen; %d combined at FDR < %.2f (two sources: %d)",
  nrow(three), sum(three$n_sources == 3), sum(is.na(three$nextgen_p)), sum(is.na(three$organoid_p)), 100 * share[1], 100 * share[2], 100 * share[3], sum(three$combined_hit), fdr_cut, sum(three$two_source_hit, na.rm = TRUE))
print(three[combined_hit == TRUE | cell_hit == TRUE | organoid_hit %in% TRUE | nextgen_hit %in% TRUE, .(symbol, cell_estimate, cell_p, organoid_p, nextgen_estimate, nextgen_p, combined_p, combined_fdr, evidence, two_source_hit, chronos_difference, chronos_difference_p)])

# --- 3. inside a cancer type, source by source --------------------------------------------------------
pairs <- as.data.table(np$within_type)
two_pairs <- as.data.table(ext_params$combine$pairs)
within_type <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  two <- two_pairs[label == pairs$label[i]]
  # The cell-line side is the stored screen; a pair refitted in stage 28 (MSI) would need that refit here.
  stopifnot(nrow(two) == 1, !two$cell_line_msi_covariate)
  stored <- two_type[pair == pairs$label[i]]
  cell <- screen[context_id == two$cell_line_context & selective == TRUE & !is.na(inside_p), .(feature, symbol, cell_inside_rank = inside_rank, cell_inside_p = inside_p, cell_n_negative = n_neg_in, cell_n_positive = n_pos_in)]
  organoid <- organoid_types[type == two$organoid_type & !is.na(feature) & !is.na(inside_p), .(feature, organoid_inside_rank = inside_rank, organoid_inside_p = inside_p, organoid_unit = inside_unit)]
  nextgen <- ng_types[type == pairs$nextgen_lineage[i] & !is.na(inside_p), .(feature, nextgen_inside_rank = inside_rank, nextgen_inside_p = inside_p, nextgen_unit = inside_unit, nextgen_inside = inside, nextgen_selective = selective, nextgen_n_negative = n_neg_in, nextgen_n_positive = n_pos_in)]
  all <- merge(merge(cell, organoid, by = "feature", all.x = TRUE), nextgen, by = "feature", all.x = TRUE)
  all <- all[!is.na(organoid_inside_p) | !is.na(nextgen_inside_p)]
  all[, cell_unit := unique(stored$cell_unit)]
  z <- cbind(signed_z(all$cell_inside_p, all$cell_inside_rank), signed_z(all$organoid_inside_p, all$organoid_inside_rank), signed_z(all$nextgen_inside_p, all$nextgen_inside_rank))
  weights <- cbind(1 / all$cell_unit, 1 / all$organoid_unit, 1 / all$nextgen_unit)
  weights[is.na(z)] <- 0
  at <- match(stored$feature, all$feature)
  stopifnot(!anyNA(at), max(abs(z[at, 1] - stored$cell_z)) < 1e-12, max(abs(z[at, 2] - stored$organoid_z)) < 1e-12)
  pooled <- stouffer(z, weights)
  all[, c("cell_z", "organoid_z", "nextgen_z", "combined_z", "combined_p", "n_sources") := .(z[, 1], z[, 2], z[, 3], pooled$z, pooled$p, pooled$n_sources)]
  all[, combined_fdr := bh(combined_p)]
  all[, combined_hit := combined_fdr < fdr_cut]
  all[, c("n_supporting", "evidence") := .(NA_integer_, NA_character_)]
  if (any(all$combined_hit)) {
    support <- evidence_summary(z[all$combined_hit, , drop = FALSE], weights[all$combined_hit, , drop = FALSE], alpha = fdr_cut)
    all[combined_hit == TRUE, c("n_supporting", "evidence") := .(support$n_supporting, support$label)]
  }
  at <- match(all$feature, stored$feature)
  all[, c("two_source_p", "two_source_fdr", "two_source_hit") := stored[at, .(combined_p, combined_fdr, combined_hit)]]
  all[, c("organoid_n_negative", "organoid_n_positive") := .(stored$organoid_n_negative[1], stored$organoid_n_positive[1])]
  all[, pair := pairs$label[i]]
  all
}))
setorder(within_type, combined_p)
write_ext(within_type, "nextgen_combined_within_type")
print(within_type[, .(genes = .N, in_all_three = sum(n_sources == 3), combined_hits = sum(combined_hit), min_combined_p = min(combined_p), min_combined_fdr = min(combined_fdr), two_source_hits = sum(two_source_hit, na.rm = TRUE)), by = pair])
print(within_type[, .SD[1:8], by = pair][, .(pair, symbol, cell_inside_p, organoid_inside_p, nextgen_inside_p, combined_p, combined_fdr, evidence, two_source_fdr)])

# --- 4. the pooled test ----------------------------------------------------------------------------------
# Strata as stage 31 counted them: the library in the cell lines and at Sanger, the screen format in NextGen.
cell_frame[, pooled_stratum := paste("Cell lines", fifelse(LibraryKY == 1, "KY", "no KY"), fifelse(LibraryHumagneCD == 1, "Humagne-CD", "no Humagne-CD"), sep = " | ")]
organoid_primary[, lineage := np$sanger_lineage$lineage[match(OncotreeLineage, np$sanger_lineage$type)]]
organoid_symbol <- organoid_genes[in_scan == TRUE & !is.na(feature), setNames(symbol, feature)]
source_features <- list(`Cell lines` = cell_genes[in_scan == TRUE, feature], `Sanger organoids` = names(organoid_symbol), NextGen = if (trial) ng_pan$feature else ng_genes[in_scan == TRUE, feature])
samples_of <- function(type) rbind(
  cell_frame[OncotreeLineage == type, .(source = "Cell lines", ModelID, stratum = pooled_stratum, exposed)],
  organoid_primary[lineage == type, .(source = "Sanger organoids", ModelID, stratum = paste("Sanger organoids", library, sep = " | "), exposed = y_negative)],
  ng_primary[OncotreeLineage == type, .(source = "NextGen", ModelID, stratum = paste("NextGen", ScreenType, sep = " | "), exposed = y_negative)]
)
# Sources that hold both Y groups in at least one stratum; the others add nothing to the test.
informative_sources <- function(samples) samples[, .(both = any(exposed == 1) && any(exposed == 0)), by = .(source, stratum)][both == TRUE, unique(source)]
scores_of <- function(samples, features) {
  parts <- lapply(split(samples, by = "source"), function(rows) switch(rows$source[1],
    `Cell lines` = t(cell_effect[rows$ModelID, features, drop = FALSE]),
    `Sanger organoids` = t(organoid_effect[rows$ModelID, organoid_symbol[features], drop = FALSE]),
    NextGen = t(ng_effect[rows$ModelID, features, drop = FALSE])
  ))
  out <- do.call(cbind, parts[unique(samples$source)])
  dimnames(out) <- list(features, samples$ModelID)
  out
}
# Features every listed source scans, with a score in each of the samples.
complete_features <- function(samples) {
  features <- Reduce(intersect, source_features[unique(samples$source)])
  features[rowSums(is.na(scores_of(samples, features))) == 0]
}
stratified_estimate <- function(mat, exposed, strata) {
  centred <- qr.resid(qr(full_rank(lineage_design(strata))), exposed)
  drop(mat %*% centred) / sum(centred^2)
}
pooled_test <- function(samples, features, seed) {
  raw <- scores_of(samples, features)
  ranked <- rank_within_source(raw, samples$source)
  test <- stratified_permutation_test(ranked, samples$exposed, samples$stratum, max_arrangements = np$pooled$max_arrangements, n_random = np$pooled$random_arrangements, seed = seed)
  test$table[, estimate_chronos := stratified_estimate(raw, samples$exposed, samples$stratum)]
  test
}
describe <- function(test) data.table(
  exact = test$exact, arrangements = test$arrangements, floor = test$floor, strata = nrow(test$strata),
  n_negative = sum(test$strata$n_exposed), n_positive = sum(test$strata$n - test$strata$n_exposed)
)
# One named gene over the samples that have a score for it; nothing when no stratum holds both Y groups.
named_test <- function(samples, feature, seed) {
  has <- samples[vapply(source, function(one) feature %in% source_features[[one]], logical(1))]
  if (nrow(has) > 0) has <- has[!is.na(scores_of(has, feature)[1, ])]
  has <- has[source %in% informative_sources(has)]
  if (nrow(has) == 0) return(NULL)
  test <- pooled_test(has, feature, seed)
  # The Sanger scores are not Chronos scores, so a pool that holds them has no estimate in Chronos units.
  if ("Sanger organoids" %in% has$source) test$table[, estimate_chronos := NA_real_]
  cbind(data.table(symbol = gene_symbol(feature), sources = paste(unique(has$source), collapse = ", ")), test$table[, .(feature, estimate_rank = estimate, estimate_chronos, statistic, p)], describe(test))
}
shuffled_in_trial <- function(samples) {
  if (trial) {
    set.seed(ext_params$seed)
    samples[, exposed := permute_within(exposed, stratum)]
  }
  samples
}
pooled_names <- type_counts[pooled == TRUE, unique(type)]
pooled_runs <- lapply(setNames(pooled_names, pooled_names), function(type) {
  samples <- samples_of(type)
  samples <- shuffled_in_trial(samples[source %in% informative_sources(samples)])
  list(samples = samples, everything = pooled_test(samples, complete_features(samples), ext_params$seed + 33))
})
pooled_table <- rbindlist(lapply(pooled_names, function(type) {
  run <- pooled_runs[[type]]
  out <- copy(run$everything$table)
  setnames(out, c("estimate", "p_family"), c("estimate_rank", "p_family_wise"))
  out[, symbol := gene_symbol(feature)]
  out[, named := feature %in% named_features]
  out[, fdr := bh(p)]
  cbind(type = type, out, describe(run$everything))
}))
pooled_table <- pooled_table[order(type, p, -abs(statistic))]
write_ext(pooled_table, "nextgen_pooled_type")
# The four paired genes: pooled, and each source on its own. BH among the four within each row set.
pooled_named <- rbindlist(lapply(pooled_names, function(type) {
  samples <- pooled_runs[[type]]$samples
  rows <- rbindlist(c(
    lapply(seq_along(named_features), function(k) cbind(source = "Pooled", named_test(samples, named_features[k], ext_params$seed + 33 + k))),
    unlist(lapply(unique(samples$source), function(one) lapply(seq_along(named_features), function(k) {
      out <- named_test(samples[source == one], named_features[k], ext_params$seed + 33 + k)
      if (is.null(out)) NULL else cbind(source = one, out)
    })), recursive = FALSE)
  ))
  rows[, named_fdr := bh(p), by = source]
  cbind(type = type, rows)
}))
write_ext(pooled_named, "nextgen_pooled_named")
for (type in pooled_names) {
  test <- pooled_runs[[type]]$everything
  log_msg("pooled %s: %d genes, %d against %d models in %d strata, %s arrangements, smallest possible p %.3g; smallest p %.3g, smallest family-wise p %.3g; BH critical p for one discovery %.3g",
    type, nrow(test$table), sum(test$strata$n_exposed), sum(test$strata$n - test$strata$n_exposed), nrow(test$strata),
    format(test$arrangements, big.mark = ","), test$floor, min(test$table$p), min(test$table$p_family), bh_next_critical(test$table$p, fdr_cut))
}
print(pooled_named[, .(type, source, symbol, estimate_rank = round(estimate_rank, 3), estimate_chronos = round(estimate_chronos, 3), p = signif(p, 3), named_fdr = signif(named_fdr, 3), floor = signif(floor, 3), n_negative, n_positive)])

# --- 5. the paired genes in lineages that the sources test alone (secondary, no calls) -----------------------
sensitivity <- rbindlist(lapply(np$pooled$sensitivity_lineages, function(type) {
  samples <- shuffled_in_trial(samples_of(type))
  cbind(type = type, rbindlist(lapply(seq_along(named_features), function(k) named_test(samples, named_features[k], ext_params$seed + 33 + k))))
}))
write_ext(sensitivity, "nextgen_pooled_sensitivity")
print(sensitivity[, .(type, symbol, sources, estimate_rank = round(estimate_rank, 3), estimate_chronos = round(estimate_chronos, 3), p = signif(p, 3), exact, n_negative, n_positive)])

# --- 6. donors screened at both institutes ---------------------------------------------------------------------
shared <- ng_cohort[shared_with_sanger == TRUE]
at <- match(shared$sanger_model, organoid_cohort$ModelID)
stopifnot(!anyNA(at))
donors <- shared[, .(ModelID, name, sanger_model, same_name = name == sanger_model, sex = Sex, lineage = OncotreeLineage, screen_format = ScreenType, nextgen_rna = rna_call, nextgen_dna = dna_call, nextgen_tier = y_tier, nextgen_call = y_call,
  nextgen_rna_y_mean = rna_y_mean, nextgen_dna_y_ratio = dna_y_ratio)]
donors[, c("sanger_rna", "sanger_dna", "sanger_tier", "sanger_call", "sanger_rna_y_mean", "sanger_dna_y_ratio", "sanger_in_cohort") := organoid_cohort[at, .(rna_call, dna_call, y_tier, y_call, rna_y_mean, dna_y_ratio, in_primary)]]
sanger_donor <- hcmi_donor(organoid_cohort$ModelID)
donors[, sanger_models_of_donor := vapply(shared$hcmi_donor, function(one) paste(organoid_cohort$ModelID[which(sanger_donor == one)], collapse = "; "), character(1))]
donors[, calls_agree := fifelse(sex == "Male" & !is.na(nextgen_call) & !is.na(sanger_call), nextgen_call == sanger_call, NA)]
# Dependency profiles: each gene centred within its own source, then the rank correlation over genes
# between a NextGen model and a Sanger organoid, against pairs of the same lineage from different donors.
common <- intersect(ng_genes[in_scan == TRUE, feature], names(organoid_symbol))
ng_models <- ng_cohort[malignant == TRUE]
organoid_models <- organoid_cohort[malignant == TRUE]
organoid_models[, lineage := np$sanger_lineage$lineage[match(OncotreeLineage, np$sanger_lineage$type)]]
centred_ranks <- function(mat) apply(scale(mat, center = TRUE, scale = FALSE), 1, rank)
profile <- cor(centred_ranks(ng_effect[ng_models$ModelID, common]), centred_ranks(organoid_effect[organoid_models$ModelID, organoid_symbol[common]]))
same_lineage <- outer(ng_models$OncotreeLineage, organoid_models$lineage, function(a, b) !is.na(a) & !is.na(b) & a == b)
same_donor <- outer(ng_models$hcmi_donor, hcmi_donor(organoid_models$ModelID), function(a, b) !is.na(a) & !is.na(b) & a == b)
dimnames(same_lineage) <- dimnames(same_donor) <- dimnames(profile)
unrelated <- profile[same_lineage & !same_donor]
donors[, profile_correlation := profile[cbind(match(ModelID, rownames(profile)), match(sanger_model, colnames(profile)))]]
donors[, unrelated_above := vapply(profile_correlation, function(value) mean(unrelated >= value), numeric(1))]
donors[, c("unrelated_pairs", "unrelated_median", "unrelated_q95", "unrelated_max") := .(length(unrelated), median(unrelated), quantile(unrelated, 0.95), max(unrelated))]
donors[, above_unrelated_q95 := profile_correlation > unrelated_q95]
write_ext(donors, "nextgen_shared_donors")
print(donors[, .(name, sanger_model, same_name, sex, nextgen_tier, nextgen_call, sanger_tier, sanger_call, calls_agree, profile_correlation = round(profile_correlation, 3), unrelated_above = signif(unrelated_above, 2), sanger_models_of_donor)])
log_msg("shared donors: %d (%d male); Y calls agree in %d of %d males called at both; profile correlation median %.3f against %.3f for %d unrelated pairs of the same lineage (largest %.3f)",
  nrow(donors), sum(donors$sex == "Male"), sum(donors$calls_agree, na.rm = TRUE), sum(!is.na(donors$calls_agree)), median(donors$profile_correlation), median(unrelated), length(unrelated), max(unrelated))

# --- the X partners across sources ---------------------------------------------------------------------------------
partner_rows <- rbind(
  two_partners[, .(symbol, y_partner, source = fifelse(source == "CRISPR combined", "CRISPR combined, two sources", fifelse(source == "Organoids, CRISPR", "Sanger organoids, CRISPR", source)), z, unit, p, n_negative, n_positive)],
  ng_pan[!is.na(y_partner), .(symbol, y_partner, source = "NextGen, CRISPR", z = signed_z(p, estimate), unit = ng_unit, p, n_negative, n_positive)],
  three[!is.na(y_partner), .(symbol, y_partner, source = "CRISPR combined, three sources", z = combined_z, unit = vapply(seq_len(.N), function(i) pooled_unit(units[!is.na(c(cell_z[i], organoid_z[i], nextgen_z[i]))]), numeric(1)), p = combined_p, n_negative = NA_real_, n_positive = NA_real_)]
)
partner_rows[, smd := z * unit]
partner_rows[, c("smd_low", "smd_high") := .(smd - interval * unit, smd + interval * unit)]
setorder(partner_rows, symbol, source)
write_ext(partner_rows, "nextgen_xy_partners")
print(dcast(partner_rows[symbol %in% np$named_genes], symbol ~ source, value.var = "smd"))

# --- follow-up of every combined hit ---------------------------------------------------------------------------------
# Designed after the combined hits of part 1 were seen: two of them are carried by no single source,
# and the plan had not said how such a hit is followed up. The checks of the single-source stages are
# applied to the combined statistic: Y status shuffled within stratum in every source at once, each
# model left out in turn, the RNA-only and the DNA-confirmed calls, and each source's adjusted model.
# Least squares in each source (weighted where its scan was), pooled with the weights of part 1.
# A hit is one of thousands of genes, so each check is also read against the p-value the hit needs
# to stay one at its rank (`bh_cut`), which takes far more shuffles than the 20,000 of a named test.
# Checks may only downgrade a hit.
source_label <- c(cell_lines = "Cell lines", organoids = "Sanger organoids", nextgen = "NextGen")
organoid_weights <- read_ext("organoid_weights", dir = ext_dir)[, setNames(weight, ModelID)]
ng_weights <- read_ext("nextgen_weights")[, setNames(weight, ModelID)]
source_frames <- function(definition) {
  organoid <- organoid_cohort[malignant == TRUE & male == TRUE & !is.na(get(definition))][order(-qc_order)][!duplicated(PatientID)][order(ModelID)]
  organoid[, exposed := as.numeric(get(definition) == "Y_negative")]
  nextgen <- ng_cohort[eligible == TRUE & !is.na(get(definition))]
  nextgen <- nextgen[best_of_group(ModelID, PatientID, ScreenNNMD)][order(ModelID)]
  nextgen[, exposed := as.numeric(get(definition) == "Y_negative")]
  nextgen[, OncotreeLineage := stratum]
  if (trial) {
    # A trial never shows a real NextGen result: its Y status is shuffled here as well.
    set.seed(ext_params$seed)
    nextgen[, exposed := permute_within(exposed, OncotreeLineage)]
  }
  list(cell_lines = analysis_frame(cell_cohort, cell_effect, exposure = definition), organoids = organoid, nextgen = nextgen)
}
primary_frames <- source_frames("y_call")
stopifnot(
  identical(primary_frames$cell_lines$ModelID, cell_frame$ModelID), setequal(primary_frames$organoids$ModelID, organoid_primary$ModelID),
  identical(primary_frames$nextgen$ModelID, ng_primary$ModelID)
)
# What one source's least-squares fit of one gene needs; nothing when the source did not enter the combination.
source_inputs <- function(name, frame, target) {
  values <- switch(name, cell_lines = cell_effect[frame$ModelID, target], organoids = organoid_effect[frame$ModelID, organoid_symbol[[target]]], nextgen = ng_effect[frame$ModelID, target])
  design <- build_design(frame, tested = "exposed", covariates = switch(name, cell_lines = cell_covariates, organoids = c("LibraryV1", "MSI"), nextgen = character()))
  if (!design$estimable) return(NULL)
  rows <- frame[design$rows]
  list(
    y = unname(values[design$rows]), x = rows$exposed, z = design$design[, colnames(design$design) != "exposed", drop = FALSE], strata = rows$OncotreeLineage,
    weights = switch(name, cell_lines = NULL, organoids = unname(organoid_weights[rows$ModelID]), nextgen = unname(ng_weights[rows$ModelID]))
  )
}
plain_z <- function(inputs) {
  root <- if (is.null(inputs$weights)) rep(1, length(inputs$y)) else sqrt(inputs$weights)
  fit <- fit_ols(inputs$y * root, inputs$x * root, inputs$z * root)
  c(z = signed_z(fit$p, fit$estimate), n = fit$n)
}
pool_z <- function(z, names_in) sum(z / units[names_in]) / sqrt(sum(1 / units[names_in]^2))
# `target` is a gene header; the name avoids the tables' own `feature` column inside their brackets.
follow_combined <- function(target) {
  row <- three[which(three$feature == target)]
  stopifnot(nrow(row) == 1)
  available <- names(units)[!is.na(c(row$cell_z, row$organoid_z, row$nextgen_z))]
  parts <- lapply(setNames(available, available), function(name) {
    inputs <- source_inputs(name, primary_frames[[name]], target)
    followup_statistics(inputs$y, inputs$x, inputs$z, inputs$strata, weights = inputs$weights, n_perm = if (trial) 2000 else np$combined_permutations, seed = ext_params$seed + 33 + match(name, names(units)))
  })
  combined <- combine_followup(parts, as.list(1 / units))
  by_source <- data.table(check = "least squares in one source", variant = source_label[available], z = unname(combined$source_z[available]), p = 2 * pnorm(-abs(combined$source_z[available])), n = vapply(parts, function(part) part$n, numeric(1)))
  by_definition <- rbindlist(lapply(c(rna_only = "y_call_rna", dna_confirmed = "y_call_concordant"), function(definition) {
    frames <- source_frames(definition)
    each <- vapply(available, function(name) plain_z(source_inputs(name, frames[[name]], target)), numeric(2))
    data.table(check = "Y-status definition", z = pool_z(each["z", ], available), n = sum(each["n", ]))
  }), idcol = "variant")
  by_definition[, p := 2 * pnorm(-abs(z))]
  stored_z <- function(table) {
    at <- which(table$feature == target)
    if (length(at) == 1) signed_z(table$p_adjusted[at], table$estimate_adjusted[at]) else NA_real_
  }
  stored <- c(cell_lines = stored_z(cell_pan), organoids = stored_z(organoid_pan), nextgen = stored_z(ng_pan))[available]
  adjusted <- data.table(check = "covariates", variant = "each source's adjusted model", z = pool_z(stored, available), n = NA_real_)
  adjusted[, p := 2 * pnorm(-abs(z))]
  summary_rows <- data.table(
    check = c("all sources, least squares", "permutation within stratum", "leave one model out", "leave one model out"),
    variant = c("pooled", sprintf("%d shuffles in every source", length(parts[[1]]$perm_t)), "smallest |z|", "largest p"),
    z = c(combined$z, combined$z, sign(combined$z) * combined$loo_min_abs, NA_real_), p = c(combined$p, combined$permutation_p, NA_real_, combined$loo_max_p), n = sum(by_source$n)
  )
  out <- rbind(by_source, summary_rows, by_definition[, .(check, variant, z, p, n)], adjusted[, .(check, variant, z, p, n)])
  sign_holds <- function(name) isTRUE(sign(by_definition[variant == name, z]) == sign(combined$z))
  out[, robust := combined$permutation_p < sp$robust_permutation_p & combined$loo_sign_stable & combined$loo_max_p < sp$robust_leave_one_out_p & sign_holds("rna_only") & sign_holds("dna_confirmed")]
  out[, bh_cut := row$bh_cut]
  out[, below_cut := p < bh_cut]
  cbind(data.table(feature = target, symbol = row$symbol, named = target %in% named_features, sources = paste(source_label[available], collapse = ", ")), out)
}
hit_features <- three[combined_hit == TRUE, feature]
followed <- parallel::mclapply(hit_features, follow_combined, mc.cores = max(1, min(n_cores, length(hit_features))))
failed <- !vapply(followed, is.data.table, logical(1))
if (any(failed)) stop("follow-up of combined hits failed: ", paste(utils::head(unlist(lapply(followed[failed], as.character)), 3), collapse = " | "))
combined_followup <- rbindlist(followed)
write_ext(combined_followup, "nextgen_combined_followup")
print(combined_followup[check != "least squares in one source", .(symbol, check, variant, z = round(z, 2), p = signif(p, 3), n, bh_cut = signif(bh_cut, 3), below_cut, robust)])
print(dcast(combined_followup[check == "least squares in one source"], symbol ~ variant, value.var = "z"))

# One row per combined hit: what it needs at its rank, what the shuffles and the leave-one-out give,
# what becomes of it when a source's adjusted model stands in for its primary one, and whether the
# gene is a dependency at all.
combine_as <- function(cell = cell_pan[, .(feature, p, estimate)], organoid = organoid_pan[!is.na(feature), .(feature, p, estimate)], nextgen = ng_pan[, .(feature, p, estimate)]) {
  pick <- function(table) table[match(three$feature, table$feature)]
  z <- cbind(signed_z(pick(cell)$p, pick(cell)$estimate), signed_z(pick(organoid)$p, pick(organoid)$estimate), signed_z(pick(nextgen)$p, pick(nextgen)$estimate))
  pooled <- stouffer(z, 1 / units)
  data.table(feature = three$feature, p = pooled$p, fdr = bh(pooled$p))
}
variants <- list(
  as_reported = combine_as(),
  nextgen_with_screen_quality = combine_as(nextgen = ng_pan[, .(feature, p = p_quality_adjusted, estimate = estimate_quality_adjusted)]),
  nextgen_with_genomic_covariates = combine_as(nextgen = ng_pan[, .(feature, p = p_adjusted, estimate = estimate_adjusted)]),
  every_source_adjusted = combine_as(cell = cell_pan[, .(feature, p = p_adjusted, estimate = estimate_adjusted)], organoid = organoid_pan[!is.na(feature), .(feature, p = p_adjusted, estimate = estimate_adjusted)],
    nextgen = ng_pan[, .(feature, p = p_adjusted, estimate = estimate_adjusted)])
)
stopifnot(max(abs(log10(variants$as_reported$p) - log10(three$combined_p))) < 1e-9)
of_check <- function(pattern, column) vapply(hit_features, function(target) combined_followup[feature == target & variant %like% pattern][[column]][1], numeric(1))
at_hit <- match(hit_features, three$feature)
rnai_pan <- read_result("gene_pan_rnai")
combined_hits <- data.table(
  feature = hit_features, symbol = three$symbol[at_hit], named = hit_features %in% named_features, n_sources = three$n_sources[at_hit], evidence = three$evidence[at_hit],
  rank = frank(three$combined_p, ties.method = "max")[at_hit], combined_p = three$combined_p[at_hit], combined_fdr = three$combined_fdr[at_hit], bh_cut = three$bh_cut[at_hit],
  permutation_p = of_check("shuffles in every source", "p"), permutation_shuffles = if (trial) 2000 else np$combined_permutations, leave_one_out_largest_p = of_check("largest p", "p"),
  two_source_fdr = three$two_source_fdr[at_hit]
)
combined_hits[, permutation_below_cut := permutation_p < bh_cut]
combined_hits[, leave_one_out_below_cut := leave_one_out_largest_p < bh_cut]
for (name in names(variants)[-1]) combined_hits[, (paste0("fdr_", name)) := variants[[name]]$fdr[at_hit]]
in_cell <- cell_pan[match(hit_features, feature)]
in_nextgen <- ng_pan[match(hit_features, feature)]
combined_hits[, c("cell_estimate", "cell_mean_negative", "cell_mean_positive", "cell_dependent_negative", "cell_dependent_positive") := in_cell[, .(estimate, mean_negative, mean_positive, fraction_dependent_negative, fraction_dependent_positive)]]
combined_hits[, c("nextgen_estimate", "nextgen_mean_negative", "nextgen_mean_positive", "nextgen_dependent_negative", "nextgen_dependent_positive") := in_nextgen[, .(estimate, mean_negative, mean_positive, fraction_dependent_negative, fraction_dependent_positive)]]
combined_hits[, rnai_p := rnai_pan$p[match(symbol, rnai_pan$symbol)]]
write_ext(combined_hits, "nextgen_combined_hits")
print(combined_hits[, .(symbol, rank, combined_p = signif(combined_p, 3), bh_cut = signif(bh_cut, 3), permutation_p = signif(permutation_p, 3), permutation_below_cut, loo = signif(leave_one_out_largest_p, 3), leave_one_out_below_cut,
  two = signif(two_source_fdr, 2), quality = signif(fdr_nextgen_with_screen_quality, 2), genomic = signif(fdr_nextgen_with_genomic_covariates, 2), all_adjusted = signif(fdr_every_source_adjusted, 2), cell_dependent_negative, cell_dependent_positive)])

# --- what each test could detect at best ------------------------------------------------------------------------------
# Standardised effect found with the stated power, for a named gene (p < 0.05) and for one more
# discovery in the family as it stands; normal theory on each design. A NextGen discovery must also
# clear the family-wise bar, and a pooled p-value cannot fall below its exact floor.
power_row <- function(analysis, source, unit, df, family_p, n_negative, n_positive, bar = NA_real_, floor = NA_real_) {
  critical <- min(bh_next_critical(family_p, fdr_cut), bar, na.rm = TRUE)
  reachable <- is.na(floor) || floor <= critical
  data.table(
    analysis = analysis, source = source, n_negative = n_negative, n_positive = n_positive, unit = unit, df = df, family = sum(is.finite(family_p)),
    named = analytic_mde(unit, df, 0.05, sp$power), family_critical_p = critical, exact_floor = floor,
    genome_wide = if (reachable) analytic_mde(unit, df, critical, sp$power) else NA_real_
  )
}
residual_df <- function(rows) {
  design <- build_design(copy(rows)[, c("exposed", "OncotreeLineage") := .(y_negative, stratum)], tested = "exposed")
  nrow(design$design) - ncol(design$design)
}
stored_row <- function(name, source_name) two_power[analysis == name & source == source_name]
cell_df <- stored_row("Pan-cancer", "Cell lines")$df
organoid_df <- stored_row("Pan-cancer", "Organoids")$df
ng_df <- residual_df(ng_primary)
pan_rows <- rbind(
  power_row("All cancer types", "NextGen", ng_unit, ng_df, ng_pan$p, sum(ng_primary$y_negative == 1), sum(ng_primary$y_negative == 0), bar = scan_summary$null$pan_bar),
  power_row("All cancer types", "Two sources", stored_row("Pan-cancer", "Combined")$unit, cell_df + organoid_df, two_pan$combined_p, stored_row("Pan-cancer", "Combined")$n_negative, stored_row("Pan-cancer", "Combined")$n_positive),
  power_row("All cancer types", "Three sources", pooled_unit(units), cell_df + organoid_df + ng_df, three$combined_p,
    stored_row("Pan-cancer", "Combined")$n_negative + sum(ng_primary$y_negative == 1), stored_row("Pan-cancer", "Combined")$n_positive + sum(ng_primary$y_negative == 0))
)
type_rows <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  rows <- within_type[pair == pairs$label[i]]
  inside <- ng_primary[OncotreeLineage == pairs$nextgen_lineage[i]]
  inside_unit <- median(rows$nextgen_unit, na.rm = TRUE)
  inside_df <- residual_df(inside)
  two <- stored_row(pairs$label[i], "Combined")
  rbind(
    power_row(pairs$label[i], "NextGen", inside_unit, inside_df, ng_types[type == pairs$nextgen_lineage[i] & selective == TRUE, inside_p], sum(inside$y_negative == 1), sum(inside$y_negative == 0), bar = scan_summary$null$inside_bar),
    power_row(pairs$label[i], "Two sources", two$unit, two$df, two_type[pair == pairs$label[i], combined_p], two$n_negative, two$n_positive),
    power_row(pairs$label[i], "Three sources", pooled_unit(c(two$unit, inside_unit)), two$df + inside_df, rows$combined_p, two$n_negative + sum(inside$y_negative == 1), two$n_positive + sum(inside$y_negative == 0))
  )
}))
pooled_rows <- rbindlist(lapply(pooled_names, function(type) {
  test <- pooled_runs[[type]]$everything
  information <- test$strata[, sum(n_exposed * (n - n_exposed) / n)]
  power_row(sprintf("%s, pooled", type), paste(unique(pooled_runs[[type]]$samples$source), collapse = " and "), 1 / sqrt(information), sum(test$strata$n) - nrow(test$strata) - 1, test$table$p,
    sum(test$strata$n_exposed), sum(test$strata$n - test$strata$n_exposed), floor = test$floor)
}))
power <- rbind(pan_rows, type_rows, pooled_rows)
write_ext(power, "nextgen_power")
print(power[, .(analysis, source, n_negative, n_positive, unit = round(unit, 3), df, family, named = round(named, 2), family_critical_p = signif(family_critical_p, 3), exact_floor = signif(exact_floor, 3), genome_wide = round(genome_wide, 2))])

jsonlite::write_json(list(
  pan_genes = nrow(three), pan_in_all_three = sum(three$n_sources == 3), pan_combined_hits = sum(three$combined_hit), pan_hit_genes = three[combined_hit == TRUE, symbol],
  pan_evidence = as.list(table(three[combined_hit == TRUE, evidence])), two_source_hits = two_pan[combined_hit == TRUE, symbol],
  pan_hit_robust = as.list(unique(combined_followup[, .(symbol, robust)])[, setNames(robust, symbol)]),
  weight_share = as.list(share), units = as.list(units),
  within_type = lapply(split(within_type, by = "pair"), function(d) list(genes = nrow(d), in_all_three = sum(d$n_sources == 3), combined_hits = sum(d$combined_hit), hit_genes = d[combined_hit == TRUE, symbol], min_combined_p = min(d$combined_p), min_combined_fdr = min(d$combined_fdr))),
  pooled = lapply(setNames(pooled_names, pooled_names), function(one_type) {
    test <- pooled_runs[[one_type]]$everything
    list(
      genes = nrow(test$table), exact = test$exact, arrangements = test$arrangements, floor = test$floor, strata = test$strata,
      min_p = min(test$table$p), min_family_wise_p = min(test$table$p_family), bh_critical_p = bh_next_critical(test$table$p, fdr_cut),
      named = pooled_named[type == one_type & source == "Pooled", .(symbol, estimate_rank, p, named_fdr)]
    )
  }),
  shared_donors = list(n = nrow(donors), males = sum(donors$sex == "Male"), called_at_both = sum(!is.na(donors$calls_agree)), calls_agree = sum(donors$calls_agree, na.rm = TRUE), same_name = sum(donors$same_name),
    above_unrelated_q95 = sum(donors$above_unrelated_q95, na.rm = TRUE), compared = sum(!is.na(donors$profile_correlation)),
    profile_median = median(donors$profile_correlation), unrelated_median = median(unrelated), unrelated_max = max(unrelated), unrelated_pairs = length(unrelated)),
  trial = trial, labels_shuffled = trial, parameters_sha256 = ext_section_hash(c("shared", "nextgen"))
), file.path(ext_out_dir, "nextgen_combined_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("three sources done: %d hits over all cancer types (%s); inside a type %s; pooled %s",
  sum(three$combined_hit), paste(three[combined_hit == TRUE, symbol], collapse = ", "),
  paste(sprintf("%s %d", pairs$label, vapply(pairs$label, function(l) sum(within_type[pair == l, combined_hit]), numeric(1))), collapse = "; "), paste(pooled_names, collapse = ", "))
