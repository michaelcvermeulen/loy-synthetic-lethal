# Stage 23: independent checks on the finished outputs. Stops on any failure.
# - headline estimates re-fitted with base R lm, outside the pipeline's own fitting code
# - line counts consistent between tables
# - keys unique, deliverables present and non-empty, figures and report present
# - freeze still intact
# - stages 25 to 30: second freezes, nothing earlier changed, cohorts, headline estimates
#   with plain lm, gates, downloads, calibration, and the statements the report makes in words
# - stages 31 to 33: the NextGen freeze, nothing earlier changed, the cohort, estimates with plain lm,
#   the bar from the shuffles, the pooled test against lm and a brute-force count, and the statements

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()

checks <- list()
record <- function(name, ok, detail = "") {
  checks[[name]] <<- list(pass = isTRUE(ok), detail = detail)
  log_msg("%s %s %s", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)
}

cohort <- read_result("cohort")
genes <- read_result("gene_pan_results")
sizes <- read_result("cohort_sizes")
effect <- readRDS(file.path(matrix_dir, "crispr_effect.rds"))
frame <- cohort[in_primary == TRUE & ModelID %in% rownames(effect)]

# --- headline estimates with plain lm ------------------------------------------------
hits <- genes[associated == TRUE]
refit <- rbindlist(lapply(hits$feature, function(feature) {
  data <- data.frame(y = effect[frame$ModelID, feature], y_negative = frame$y_negative, lineage = factor(frame$OncotreeLineage), ky = frame$LibraryKY, humagne = frame$LibraryHumagneCD)
  fit <- summary(lm(y ~ lineage + ky + humagne + y_negative, data = data))$coefficients["y_negative", ]
  data.table(symbol = gene_symbol(feature), lm = fit[["Estimate"]], lm_p = fit[["Pr(>|t|)"]])
}))
refit[, pipeline := hits$estimate[match(symbol, hits$symbol)]]
write_result_stable(refit, "validation_hit_refit")
print(refit)
record("hit estimates agree with plain lm", all(abs(refit$pipeline - refit$lm) < 0.02) && all(refit$lm_p < 1e-10), sprintf("max difference %.4f", max(abs(refit$pipeline - refit$lm))))

# A raw two-group difference, with no model at all, must point the same way.
raw <- vapply(hits$feature, function(feature) mean(effect[frame[y_negative == 1, ModelID], feature], na.rm = TRUE) - mean(effect[frame[y_negative == 0, ModelID], feature], na.rm = TRUE), numeric(1))
record("hit effects have the same sign as raw group differences", all(sign(raw) == sign(hits$estimate)), paste(sprintf("%s %.2f", hits$symbol, raw), collapse = "; "))

# --- consistency between tables --------------------------------------------------------
crispr <- sizes[dataset == "crispr_effect"]
record("CRISPR group sizes agree between cohort_sizes and the scan", crispr$y_negative == max(genes$n_negative) && crispr$y_positive == max(genes$n_positive),
  sprintf("%d / %d", crispr$y_negative, crispr$y_positive))
calls <- read_result("y_status_calls")
record("one call row per model", !anyDuplicated(calls$ModelID) && nrow(calls) == nrow(fread(input_path("model"), select = "ModelID")))
record("no female or unknown-sex line carries a call", all(is.na(calls[Sex != "Male", y_call])))
record("cohort holds only eligible males", all(cohort$Sex == "Male") && !any(cohort$OncotreePrimaryDisease == "Non-Cancerous") && all(cohort$EngineeredModel == ""))
record("every cohort call matches the call table", all(cohort$y_call == calls$y_call[match(cohort$ModelID, calls$ModelID)], na.rm = TRUE))
drugs <- read_result("drug_pan_results")
record("drug tests respect the minimum group size", all(drugs$n_negative >= params$drug_scan$min_per_group & drugs$n_positive >= params$drug_scan$min_per_group))
record("FDR never below its p-value", all(genes$fdr >= genes$p - 1e-15) && all(drugs$fdr >= drugs$p - 1e-15))
context <- readRDS(file.path(derived_dir, "context_gene_screen.rds"))
binary <- context[type == "binary"]
record("every fitted binary context meets the cell gate", all(pmin(binary$n_neg_in, binary$n_pos_in, binary$n_neg_out, binary$n_pos_out) >= params$context$min_cell))

# --- the TP53-dependent control ---------------------------------------------------------
# MDM2 dependency follows TP53 status; Y-negative lines are more often TP53-altered, so a
# small Y association for MDM2 should weaken once TP53 is in the model.
mdm2 <- genes[symbol == "MDM2", .(estimate, p, estimate_adjusted, p_adjusted)]
write_result_stable(mdm2, "validation_mdm2_control")
record("MDM2 is not a Y hit and moves toward zero after TP53 adjustment", !genes[symbol == "MDM2", associated] && abs(mdm2$estimate_adjusted) <= abs(mdm2$estimate) + 0.01,
  sprintf("base %.3f (p %.2g), adjusted %.3f (p %.2g)", mdm2$estimate, mdm2$p, mdm2$estimate_adjusted, mdm2$p_adjusted))

# --- what the cancer-type tests could detect (stage 24) -----------------------------------
power <- jsonlite::fromJSON(file.path(derived_dir, "power_summary.json"))
record("power stage passed its own checks", all(unlist(power$checks)), paste(names(power$checks)[!unlist(power$checks)], collapse = "; "))
record("power parameters are the ones the stage used", identical(file_sha256(file.path(script_dir, "power_parameters.json")), power$parameters_sha256))
# The census, counted again with base R from the cohort table.
crispr_lines <- cohort[in_primary == TRUE & ModelID %in% rownames(effect)]
testable <- function(column) {
  counts <- table(crispr_lines[[column]][crispr_lines[[column]] != ""], crispr_lines$y_negative[crispr_lines[[column]] != ""])
  c(with_lines = nrow(counts), at_least_5 = sum(counts[, "0"] >= params$context$min_cell & counts[, "1"] >= params$context$min_cell))
}
recount <- vapply(c(lineage = "OncotreeLineage", disease = "OncotreePrimaryDisease", subtype = "OncotreeSubtype"), testable, numeric(2))
census <- read_result("power_census")
record("power census agrees with a recount from the cohort", all(recount["with_lines", census$level] == census$with_lines) && all(recount["at_least_5", census$level] == census$at_least_5),
  paste(sprintf("%s %d of %d", colnames(recount), recount["at_least_5", ], recount["with_lines", ]), collapse = "; "))
# The positive controls, counted again from the stored screen.
screen_lineage <- context[class == "lineage"]
screen_lineage[, own_fdr := p.adjust(inside_p, method = "BH"), by = context_id]
four <- screen_lineage[feature %in% hits$feature]
positive <- read_result("power_positive_controls")[level == "lineage"]
record("positive-control counts agree with the stored screen", nrow(four) == nrow(positive) && sum(four$inside_p < 0.05) == sum(positive$inside_p < 0.05) && sum(four$own_fdr < 0.05) == sum(positive$own_fdr < 0.05) &&
  sum(four$inside_fdr_class < 0.05) == sum(positive$inside_fdr_class < 0.05),
  sprintf("%d chances: %d nominal, %d genome-wide alone, %d in the screen's family", nrow(four), sum(four$inside_p < 0.05), sum(four$own_fdr < 0.05), sum(four$inside_fdr_class < 0.05)))
# One benchmark with plain lm on the eligible males with CRISPR data.
benchmarks <- read_result("power_benchmarks")
braf <- benchmarks[label == "BRAF in BRAF V600E"]
context_matrix <- read_result("context_matrix")
eligible <- cohort[ModelID %in% rownames(effect)]
braf_data <- data.frame(y = effect[eligible$ModelID, braf$feature], altered = context_matrix[[braf$biomarker]][match(eligible$ModelID, context_matrix$ModelID)],
  lineage = factor(eligible$OncotreeLineage), ky = eligible$LibraryKY, humagne = eligible$LibraryHumagneCD)
braf_lm <- coef(lm(y ~ lineage + ky + humagne + altered, data = braf_data))[["altered"]]
record("benchmark estimate agrees with plain lm", abs(braf_lm - braf$effect) < 1e-6, sprintf("BRAF in BRAF V600E %.4f against %.4f", braf_lm, braf$effect))
shifts <- read_result("power_shift_detection")
unshifted <- shifts[unit == "sd" & shift == 0]
record("no shift, no detection", max(unshifted$one_type) <= 0.01 && max(unshifted$as_run) <= 0.01 && abs(mean(unshifted$candidate) - 0.05) < 0.02,
  sprintf("nominal rate %.3f, largest genome-wide rate %.3f", mean(unshifted$candidate), max(unshifted$one_type)))
largest <- shifts[unit == "sd" & shift == max(shift)]
record("the largest shift is found at the nominal level in every cancer type", all(largest$candidate >= 0.95), sprintf("smallest rate %.3f over %d cancer types", min(largest$candidate), nrow(largest)))
# The four established dependencies sit in every Y-negative line, so the prediction set beside them must too.
positive_all <- read_result("power_positive_controls")
record("positive controls carry a prediction with the effect in every Y-negative line",
  "predicted_everywhere_one_type" %in% names(positive_all) && all(is.finite(positive_all$predicted_everywhere_one_type)) && all(positive_all$predicted_everywhere_one_type >= 0 & positive_all$predicted_everywhere_one_type <= 1),
  if ("predicted_everywhere_one_type" %in% names(positive_all)) sprintf("expected recoveries among lineages %.1f, observed %d", positive_all[level == "lineage", sum(predicted_everywhere_one_type)], positive_all[level == "lineage", sum(own_fdr < 0.05)]) else "column missing")
# Every gene that passes with a cancer type analysed alone, counted again from the stored screen.
cancer_screen <- context[class %in% c("lineage", "disease", "subtype")]
cancer_screen[, own_fdr := p.adjust(inside_p, method = "BH"), by = context_id]
alone_recount <- cancer_screen[own_fdr < 0.05, .(context_id, feature)]
alone_file <- file.path(derived_dir, "power_alone_hits.rds")
alone_hits <- if (file.exists(alone_file)) read_result("power_alone_hits") else NULL
record("genes passing with a cancer type analysed alone are listed", !is.null(alone_hits) && nrow(alone_hits) == nrow(alone_recount) && setequal(paste(alone_hits$context_id, alone_hits$feature), paste(alone_recount$context_id, alone_recount$feature)),
  sprintf("%d passes, %d outside the established four", nrow(alone_recount), sum(!alone_recount$feature %in% hits$feature)))
# Statements the final review found unsupported must be gone, and the ones it asked for present.
report_text <- gsub("\\s+", " ", gsub("<[^>]+>", " ", paste(readLines(file.path(results_dir, "report.html"), warn = FALSE), collapse = " ")))
retracted <- c("not the cause", "none of them can be changed by reanalysis", "would not have been, in any cancer type", "Most cancer types cannot be tested", "Two bars are reported")
required <- c("named in advance", "Monte Carlo", "CIB1", "in every Y-negative line", "STAG1 in STAG2-mutant lines apart")
record("report wording follows the final review", !any(vapply(retracted, grepl, logical(1), x = report_text, fixed = TRUE)) && all(vapply(required, grepl, logical(1), x = report_text, fixed = TRUE)),
  paste(c(sprintf("still present: %s", retracted[vapply(retracted, grepl, logical(1), x = report_text, fixed = TRUE)]), sprintf("missing: %s", required[!vapply(required, grepl, logical(1), x = report_text, fixed = TRUE)])), collapse = "; "))

# --- organoids, RNAi and primary AML (stages 25 to 30) --------------------------------------
ext_json <- function(name) jsonlite::fromJSON(file.path(ext_dir, paste0(name, ".json")))
ext_table <- function(name) read_result(name, dir = ext_dir)
problem <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))

# Second freezes, and the parameters each stage ran with.
freeze_problem <- problem({ assert_extension_freeze("organoid", dir = ext_dir); assert_extension_freeze("beataml", dir = ext_dir) })
record("extension freezes intact", freeze_problem == "", freeze_problem)
summary_sections <- list(organoid_prepare_summary = c("shared", "organoid"), organoid_scan_summary = c("shared", "organoid"), rnai_context_run_summary = c("shared", "rnai"),
  combined_summary = c("shared", "combine"), beataml_prepare_summary = c("shared", "beataml"), beataml_scan_summary = c("shared", "beataml"))
stale <- names(summary_sections)[!vapply(names(summary_sections), function(name) identical(ext_json(name)$parameters_sha256, ext_section_hash(summary_sections[[name]])), logical(1))]
record("every extension stage ran with the current extension parameters", length(stale) == 0, paste(stale, collapse = "; "))
trial_runs <- names(summary_sections)[vapply(names(summary_sections), function(name) isTRUE(ext_json(name)$trial), logical(1))]
record("no extension table comes from a trial run", length(trial_runs) == 0, paste(trial_runs, collapse = "; "))

# Nothing from before the extension changed: every earlier table and both frozen parameter files keep their hash.
baseline_lines <- readLines(file.path(ext_dir, "extension_baseline_sha256.txt"))
baseline <- data.table(sha256 = substr(baseline_lines, 1, 64), path = substr(baseline_lines, 67, nchar(baseline_lines)))
earlier <- baseline[!grepl("^results/", path)]
current <- vapply(file.path(root_dir, earlier$path), function(path) if (file.exists(path)) file_sha256(path) else NA_character_, character(1))
changed <- earlier$path[is.na(current) | current != earlier$sha256]
record("every table and parameter file from before the extension is unchanged", nrow(earlier) >= 450 && length(changed) == 0,
  sprintf("%d files compared%s", nrow(earlier), if (length(changed) > 0) paste0("; changed: ", paste(head(changed, 5), collapse = ", ")) else ""))

# The report outside the new section reads as it did before the extension.
baseline_report <- file.path(ext_dir, "extension_baseline_report.html")
current_report <- file.path(results_dir, "report.html")
before <- report_sections(baseline_report)
after <- report_sections(current_report)
untouched <- setdiff(names(before), c(names(before)[1], "Summary", "Limits", "Files"))
added <- setdiff(names(after), names(before))
record("report sections outside the extension read as before", length(untouched) >= 8 && all(untouched %in% names(after)) && all(before[untouched] == after[untouched]) && identical(added, c("Organoids, RNAi and primary AML", "DepMap NextGen models", "Driver status inside cancer types")),
  sprintf("%d sections compared; added: %s", length(untouched), paste(added, collapse = ", ")))
lost <- setdiff(report_blocks(baseline_report), report_blocks(current_report))
record("every earlier paragraph is still in the report, apart from the two edited on purpose",
  length(lost) == 2 && any(startsWith(lost, "In short")) && any(grepl("parameters.json holds every threshold", lost, fixed = TRUE)), sprintf("%d earlier paragraphs differ", length(lost)))
typed <- prose_numbers(rmd_section(readLines(file.path(script_dir, "22_report.Rmd")), "Organoids, RNAi and primary AML"))
record("every number in the new report section is computed from the tables", all(typed %in% c("2026", "10", "04", "0.05", "0.5", "0", "20,000")), sprintf("typed constants: %s", paste(unique(typed), collapse = " ")))

# Organoids: cohort, scan, controls, calibration.
org_summary <- ext_json("organoid_prepare_summary")
org_scan <- ext_json("organoid_scan_summary")
org_calls <- ext_table("organoid_y_calls")
org_cohort <- ext_table("organoid_cohort")
org_frame <- org_cohort[in_primary == TRUE]
org_pan <- ext_table("organoid_pan_results")
org_annotation <- ext_table("organoid_gene_annotation")
found <- unlist(org_summary$counts_found)
expected <- unlist(org_summary$counts_expected)[names(found)]
record("organoid counts are those fixed before the screen", all(abs(found - expected) <= ext_params$shared$count_tolerance) && nrow(org_frame) == org_summary$primary,
  sprintf("%d screened, %d male, %d in the cohort (%d Y-negative, %d Y-positive); largest gap %d", org_summary$screened, org_summary$males, nrow(org_frame), sum(org_frame$y_negative == 1), sum(org_frame$y_negative == 0), max(abs(found - expected))))
record("only concordant male cancers enter the organoid cohort, one per patient",
  all(is.na(org_calls[male == FALSE, y_call])) && org_frame[, all(male & malignant & !sex_discordant & !is.na(y_call))] && !anyDuplicated(org_frame$PatientID) && !anyDuplicated(org_frame$ModelID),
  sprintf("%d organoids from %d patients", nrow(org_frame), uniqueN(org_frame$PatientID)))
model_table <- fread(input_path("model"))
normalise_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
organoid_ids <- setdiff(unique(normalise_id(c(org_cohort$ModelID, org_cohort$individual_ID, org_cohort$PatientID))), "")
text_columns <- names(model_table)[vapply(model_table, is.character, logical(1))]
in_depmap <- unlist(lapply(text_columns, function(column) intersect(normalise_id(model_table[[column]]), organoid_ids)))
record("no organoid or donor identifier appears in the DepMap model table", length(in_depmap) == 0, sprintf("%d identifiers against %d models", length(organoid_ids), nrow(model_table)))
noise <- org_pan[, abs(depleted_fraction_difference) > 0 & abs(depleted_fraction_difference) < 1e-12]
same_call <- org_pan[, fraction_depleted_negative == fraction_depleted_positive & fraction_depleted_negative %in% c(0, 1)]
record("a gene whose depletion calls leave no Y contrast carries none", !any(noise, na.rm = TRUE) && all(is.na(org_pan$depleted_fraction_p[same_call])) && all(org_pan$depleted_fraction_difference[same_call] == 0),
  sprintf("%d genes called the same way in every organoid; %d estimates at rounding level", sum(same_call), sum(noise, na.rm = TRUE)))
record("no chrY gene is in the organoid scan", !any(org_pan$symbol %in% org_annotation[is_y_gene == TRUE, symbol]) && !any(org_pan$chromosome %in% c("Y", "chrY")), sprintf("%d genes scanned", nrow(org_pan)))
# The two leading genes of the scan, refitted with weighted lm outside the pipeline's fitting code.
org_effect <- readRDS(file.path(ext_dir, "organoid_effect.rds"))
org_weight <- ext_table("organoid_weights")
org_top <- org_pan[order(p)][1:2]
org_refit <- rbindlist(lapply(org_top$symbol, function(gene) {
  data <- data.frame(y = org_effect[org_frame$ModelID, gene], y_negative = org_frame$y_negative, type = factor(org_frame$OncotreeLineage), library = org_frame$LibraryV1, msi = org_frame$MSI,
    weight = org_weight$weight[match(org_frame$ModelID, org_weight$ModelID)])
  fit <- summary(lm(y ~ type + library + msi + y_negative, data = data, weights = weight))$coefficients["y_negative", ]
  data.table(symbol = gene, lm = fit[["Estimate"]], lm_p = fit[["Pr(>|t|)"]])
}))
org_refit[, pipeline := org_top$estimate]
write_result_stable(org_refit, "extension_validation_organoid_refit", dir = ext_dir)
record("organoid estimates agree with weighted lm", all(abs(org_refit$pipeline - org_refit$lm) < 1e-6), paste(sprintf("%s %.4f against %.4f", org_refit$symbol, org_refit$lm, org_refit$pipeline), collapse = "; "))
org_controls <- ext_table("organoid_controls")
record("organoid controls recover known biology", nrow(org_controls) == 3 && all(org_controls$pass) && all(sign(org_controls$estimate) == org_controls$expected_sign),
  paste(sprintf("%s by %s %.2f (p %.1e)", org_controls$gene, org_controls$biomarker, org_controls$estimate, org_controls$p), collapse = "; "))
org_shuffled <- ext_table("organoid_calibration")[run != "observed"]
record("organoid tests are calibrated under shuffled Y status",
  nrow(org_shuffled) == ext_params$shared$null_runs && all(org_shuffled$pan_lambda > 0.85 & org_shuffled$pan_lambda < 1.15) && all(org_shuffled$inside_lambda > 0.85 & org_shuffled$inside_lambda < 1.15) &&
    all(org_shuffled$pan_fdr_hits <= 1) && all(org_shuffled$inside_fdr_hits <= 1) && all(org_shuffled$type_specific == 0),
  sprintf("lambda %.2f to %.2f with all types together, %.2f to %.2f inside a type", min(org_shuffled$pan_lambda), max(org_shuffled$pan_lambda), min(org_shuffled$inside_lambda), max(org_shuffled$inside_lambda)))

# RNAi: the per-gene gate, and one estimate with plain lm.
rnai_summary <- ext_json("rnai_context_run_summary")
rnai_screen <- readRDS(file.path(ext_dir, "rnai_context_screen.rds"))
setDT(rnai_screen)
rnai_smallest <- pmin(rnai_screen$n_neg_in, rnai_screen$n_pos_in, rnai_screen$n_neg_out, rnai_screen$n_pos_out)
record("every RNAi test meets the cell gate on the lines with data for that gene", all(rnai_smallest >= ext_params$shared$min_cell) && nrow(rnai_screen) == rnai_summary$screen_tests,
  sprintf("%s tests, smallest cell %d", format(nrow(rnai_screen), big.mark = ","), min(rnai_smallest)))
rnai <- readRDS(file.path(matrix_dir, "rnai.rds"))
rnai_top <- rnai_screen[which.min(inside_p)]
rnai_lines <- cohort[ModelID %in% rownames(rnai) & !is.na(y_call)]
rnai_data <- data.frame(y = rnai[rnai_lines$ModelID, rnai_top$feature], y_negative = as.numeric(rnai_lines$y_call == "Y_negative"), lineage = factor(rnai_lines$OncotreeLineage),
  inside = context_matrix[[rnai_top$context_id]][match(rnai_lines$ModelID, context_matrix$ModelID)])
rnai_lm <- coef(lm(y ~ lineage + y_negative, data = rnai_data[rnai_data$inside %in% 1, ]))[["y_negative"]]
record("RNAi estimate agrees with plain lm", nrow(rnai_lines) == rnai_summary$lines && abs(rnai_lm - rnai_top$inside) < 1e-6, sprintf("%s inside %s %.4f against %.4f", rnai_top$symbol, rnai_top$label, rnai_lm, rnai_top$inside))

# Combination: the stored p-values enter unchanged, and the pooled statistic is recomputed by hand.
comb_summary <- ext_json("combined_summary")
comb_pan <- ext_table("combined_pan")
comb_within <- ext_table("combined_within_type")
record("the combination uses the stored p-values of both scans", all(comb_pan$cell_p == genes$p[match(comb_pan$feature, genes$feature)]) && all(comb_pan$organoid_p == org_pan$p[match(comb_pan$organoid_symbol, org_pan$symbol)]),
  sprintf("%s genes", format(nrow(comb_pan), big.mark = ",")))
pairs <- as.data.table(ext_params$combine$pairs)
comb_within[, context_id := pairs$cell_line_context[match(pair, pairs$label)]]
stored_inside <- context[comb_within, on = c("feature", "context_id"), inside_p]
unadjusted <- comb_within[msi_adjusted_cell_lines == FALSE]
record("cell-line tests inside a cancer type reproduce the stored screen", !anyNA(stored_inside) && all(abs(log10(comb_within$cell_inside_p_stored_model) - log10(stored_inside)) <= ext_params$combine$reproduce_log10_max) &&
  all(unadjusted$cell_inside_p == unadjusted$cell_inside_p_stored_model) && all(unlist(comb_summary$reproduce_gap_log10) <= ext_params$combine$reproduce_log10_max),
  sprintf("%s tests, largest gap %.1e in log10 p", format(nrow(comb_within), big.mark = ","), max(abs(log10(comb_within$cell_inside_p_stored_model) - log10(stored_inside)))))
by_hand <- comb_pan[order(combined_p)][1:20]
by_hand[, cell_unit := cell_smd / cell_z]
by_hand[, organoid_unit := organoid_smd / organoid_z]
by_hand[, z := (cell_z / cell_unit + organoid_z / organoid_unit) / sqrt(1 / cell_unit^2 + 1 / organoid_unit^2)]
record("the pooled statistic equals the weighted sum recomputed by hand",
  all(abs(by_hand$cell_z - sign(by_hand$cell_estimate) * qnorm(by_hand$cell_p / 2, lower.tail = FALSE)) < 1e-8) && all(abs(by_hand$z - by_hand$combined_z) < 1e-8) &&
    all(abs(log10(2 * pnorm(-abs(by_hand$z))) - log10(by_hand$combined_p)) < 1e-8), sprintf("20 leading genes, largest gap %.1e", max(abs(by_hand$z - by_hand$combined_z))))

# BeatAML2: one specimen per patient, the control, and the leading slope with plain lm.
aml_cohort <- ext_table("beataml_cohort")
aml_frame <- aml_cohort[in_primary == TRUE & !is.na(age)][order(SpecimenID)]
aml_drugs <- ext_table("beataml_drug_results")
aml_scan <- ext_json("beataml_scan_summary")
record("one BeatAML2 specimen per patient, eligible males only", !anyDuplicated(aml_cohort[in_primary == TRUE, PatientID]) && aml_cohort[in_primary == TRUE, all(Sex == "Male" & eligible & has_drug)] && nrow(aml_frame) == aml_scan$specimens,
  sprintf("%d patients, %d with age in the screen", aml_cohort[in_primary == TRUE, .N], nrow(aml_frame)))
aml_controls <- ext_table("beataml_controls")
record("FLT3 inhibitors are more active in FLT3-ITD specimens", nrow(aml_controls) > 0 && all(aml_controls$estimate < 0) && all(aml_controls$pass) && min(aml_controls$p) < 0.001,
  paste(sprintf("%s %.2f (p %.1e)", aml_controls$inhibitor, aml_controls$estimate, aml_controls$p), collapse = "; "))
aml_auc <- readRDS(file.path(ext_dir, "beataml_auc.rds"))
aml_top <- aml_drugs[which.min(p)]
aml_data <- data.frame(y = rank_normalise(aml_auc[aml_frame$SpecimenID, aml_top$inhibitor]), f = aml_frame$f, specimen = factor(aml_frame$specimen_type), stage = factor(aml_frame$stage), wave = factor(aml_frame$wave),
  age = aml_frame$age, fusion = aml_frame$cbf_fusion)
aml_lm <- summary(lm(y ~ specimen + stage + wave + age + fusion + f, data = aml_data))$coefficients["f", ]
record("BeatAML2 leading slope agrees with plain lm", abs(aml_lm[["Estimate"]] - aml_top$estimate) < 1e-6 && abs(log10(aml_lm[["Pr(>|t|)"]]) - log10(aml_top$p_classical)) < 1e-6,
  sprintf("%s %.4f against %.4f", aml_top$inhibitor, aml_lm[["Estimate"]], aml_top$estimate))
aml_shuffled <- ext_table("beataml_calibration")[run != "observed"]
record("BeatAML2 gives no hit under shuffled Y loss", nrow(aml_shuffled) == ext_params$shared$null_runs && all(aml_shuffled$fdr_hits == 0), paste(aml_shuffled$fdr_hits, collapse = ", "))

# The public inputs added for these stages.
inputs <- as.data.table(ext_params$inputs)
download_problem <- problem(for (i in seq_len(nrow(inputs))) verify_download(file.path(root_dir, inputs$file[i]), inputs$checksum[i], inputs$algorithm[i]))
record("the added public inputs match their recorded checksums", download_problem == "", if (download_problem == "") sprintf("%d files", nrow(inputs)) else download_problem)

# Statements the report makes in words, checked against the tables they come from.
org_hits <- org_pan[associated == TRUE]
rnai_inside <- ext_table("rnai_context_inside_only")
comb_hits <- comb_pan[combined_hit == TRUE]
aml_hits <- aml_drugs[associated == TRUE]
aml_follow <- ext_table("beataml_hit_followup")
aml_families <- ext_table("beataml_family_results")
ras <- ext_table("beataml_ras_by_y_loss")
leads <- ext_table("combined_cell_line_leads_in_rnai")
seen <- ext_table("combined_seen_before")
mek <- aml_drugs[grepl("Trametinib|CI-1040", inhibitor)]
oes_lines <- cohort[ModelID %in% rownames(effect) & !is.na(y_negative) & OncotreePrimaryDisease == "Esophagogastric Adenocarcinoma"]
rectal_lines <- cohort[ModelID %in% rownames(effect) & !is.na(y_negative) & OncotreeSubtype == "Rectal Adenocarcinoma"]
colorectal_organoids <- org_frame[OncotreeLineage == "Colorectal"]
statements <- c(
  "the organoid hits are DDX3X and RPS4X, more needed in Y-negative organoids and LOY-specific" = setequal(org_hits$symbol, c("DDX3X", "RPS4X")) && all(org_hits$estimate < 0) && all(org_hits$label == "LOY-specific") && all(org_hits$seen_before),
  "both organoid hits hold in every follow-up check" = ext_table("organoid_hit_followup")[is_hit == TRUE, .N > 0 && all(robust) && all(sign(estimate[!is.na(estimate)]) == -1)],
  "ZFX does not pass in organoids and EIF1AX is not in their gene set" = org_pan[symbol == "ZFX", !associated && p > 0.05] && !"EIF1AX" %in% org_pan$symbol && comb_pan[symbol == "ZFX", smd_difference_p < 0.05],
  "CORO1A has no support in any source" = seen[symbol == "CORO1A", organoid_fdr > 0.05 && cell_p > 0.05 && rnai_p > 0.05],
  "nothing in organoids is specific to a cancer type or a genotype" = org_scan$type_specific == 0 && org_scan$inside_only == 0 && org_scan$genotype_specific == 0,
  "the genotype detections in organoids are DDX3X and RPS4X" = all(ext_table("organoid_genotype_screen")[detected_inside_only == TRUE, symbol] %in% c("DDX3X", "RPS4X")),
  "Y-negative organoids have the better screens, and common-essential depletion does not differ" = org_frame[y_negative == 1, median(TPR)] > org_frame[y_negative == 0, median(TPR)] && org_scan$essential_shift$p > 0.05,
  "RNAi makes no context-specific call" = rnai_summary$screen_context_specific == 0 && rnai_summary$focused_context_specific == 0,
  "RNAi detections inside a context are RPS4X and DDX3X, outside cancer types" = all(rnai_inside$symbol %in% c("RPS4X", "DDX3X")) && !any(rnai_inside$class %in% c("lineage", "disease", "subtype")),
  "the leads RNAi confirms are paired genes, and CIB1 in bladder cannot be tested there" = leads[tested_in_rnai == TRUE, all(symbol %in% c("DDX3X", "RPS4X", "ZFX"))] && leads[symbol == "CIB1", .N == 1 && !tested_in_rnai] && !"EIF1AX" %in% colnames(rnai),
  "the combined hits are DDX3X, RPS4X and ZFX, the last carried by the cell lines" = setequal(comb_hits$symbol, c("DDX3X", "RPS4X", "ZFX")) && comb_hits[symbol == "ZFX", evidence == "carried by one source"] && all(comb_hits[symbol != "ZFX", evidence] == "supported by both"),
  "the standardised effects differ between sources for DDX3X and not for RPS4X" = comb_pan[symbol == "DDX3X", smd_difference_p < 0.05 && abs(organoid_smd) < abs(cell_smd)] && comb_pan[symbol == "RPS4X", smd_difference_p > 0.05],
  "no combined hit inside a cancer type" = !any(comb_within$combined_hit),
  "MSI is more common among Y-positive colorectal organoids" = colorectal_organoids[y_negative == 0, mean(MSI)] > colorectal_organoids[y_negative == 1, mean(MSI)],
  "all oesophageal adenocarcinoma and rectal cell lines are Y-negative" = oes_lines[OncotreeSubtype == "Esophageal Adenocarcinoma", .N > 0 && all(y_negative == 1)] && rectal_lines[, .N > 0 && all(y_negative == 1)] &&
    !any(grepl("Esophageal", oes_lines[y_negative == 0, OncotreeSubtype])),
  "one inhibitor passes in BeatAML2, selumetinib, less active with Y loss and LOY-specific" = nrow(aml_hits) == 1 && grepl("^Selumetinib", aml_hits$inhibitor) && aml_hits$estimate > 0 && aml_hits$label == "LOY-specific" && all(aml_follow[is_hit == TRUE, robust]),
  "the other MEK inhibitors point the same way without passing, and no family passes" = nrow(mek) == 2 && all(mek$estimate > 0) && all(mek$p > 0.05) && !any(aml_families$fdr < 0.05),
  "no inhibitor is more active with Y loss" = !any(aml_drugs[estimate < 0, fdr] < 0.05) && aml_scan$general_response$p > 0.05,
  "no patient with Y loss carries a RAS-pathway mutation, and the adjusted slope is smaller and would not pass" = ras[y_loss != "no Y loss", sum(ras_pathway_mutant)] == 0 &&
    aml_follow[grepl("^plus RAS-pathway", variant), estimate] < aml_follow[grepl("^main", variant), estimate] && all(aml_follow[check == "added after the result", p] > ext_table("beataml_power")$family_critical_p),
  "selumetinib has no single-agent match in the cell-line screens" = !any(grepl("^Selumetinib", ext_table("beataml_cell_line_comparison")$inhibitor))
)
record("statements the new section makes in words hold in the tables", !anyNA(statements) && all(statements), paste(names(statements)[is.na(statements) | !statements], collapse = "; "))
new_text <- after[["Organoids, RNAi and primary AML"]]
new_required <- c("added after the hit was seen", "both detections are withdrawn", "not independent replication", "best-case figures")
record("the new section states its two changes after the fact and its limits", all(vapply(new_required, grepl, logical(1), x = new_text, fixed = TRUE)), paste(new_required[!vapply(new_required, grepl, logical(1), x = new_text, fixed = TRUE)], collapse = "; "))

# Findings of the fresh review of 2026-10-04: the facts behind each corrected sentence, and the wording itself.
aml_lead <- aml_drugs[which.min(p)]
ras_high <- ras[y_loss == "f at least 0.5"]
ras_rest <- ras[y_loss != "f at least 0.5"]
ras_fisher <- fisher.test(matrix(c(ras_high$ras_pathway_mutant, ras_high$ras_pathway_wild_type, sum(ras_rest$ras_pathway_mutant), sum(ras_rest$ras_pathway_wild_type)), 2))$p.value
led <- org_cohort[male == TRUE & y_tier == "rna_led_dna_ambiguous"]
continuous <- ext_table("organoid_continuous")
review_facts <- c(
  "the hit rests on fewer patients with Y loss than the cohort holds" = aml_lead$high_loss < aml_frame[f >= ext_params$beataml$high_loss_fraction, .N] && aml_lead$any_loss >= aml_lead$high_loss,
  "the stored RAS p-value is Fisher's test of a fraction of at least 0.5 against the rest" = abs(ras_fisher - ras$fisher_p[1]) < 1e-12,
  "Y loss does not go with age in the cohort and does go with t(8;21)" = wilcox.test(age ~ (f >= ext_params$beataml$high_loss_fraction), data = aml_frame, exact = FALSE)$p.value > 0.05 &&
    fisher.test(table(aml_frame$f >= ext_params$beataml$high_loss_fraction, aml_frame$t8_21))$p.value < 0.05,
  "RNA-led organoid calls disagree with DNA in both directions" = led[rna_call == "depleted" & dna_y_retained > 0.5, .N] > 0 && led[rna_call == "expressed" & dna_y_ratio < params$dna_call$retained_ratio_min, .N] > 0,
  "the continuous organoid scan gives the same genes as the two-group scan" = setequal(continuous[fdr < ext_params$shared$fdr, symbol], org_hits$symbol)
)
record("facts behind the sentences corrected after the fresh review hold", !anyNA(review_facts) && all(review_facts), paste(names(review_facts)[is.na(review_facts) | !review_facts], collapse = "; "))
review_required <- c(sprintf("rests on %d patients", aml_lead$high_loss), "a noisier screen lowers a standardised effect", "share most of their lines", "for the first group against the rest",
  "t(8;21)", "little chrY left in DNA", "rank-sum", "continuous score", "each controlled on its own")
review_retracted <- c("is weaker in organoids than in cell lines", "where unrelated assays would give", "travels with age", sprintf("rests on %d patients", aml_frame[f >= ext_params$beataml$high_loss_fraction, .N]))
has <- function(phrases) vapply(phrases, grepl, logical(1), x = report_text, fixed = TRUE)
record("report wording follows the fresh review of the extension", all(has(review_required)) && !any(has(review_retracted)),
  paste(c(sprintf("missing: %s", review_required[!has(review_required)]), sprintf("still present: %s", review_retracted[has(review_retracted)])), collapse = "; "))

# --- DepMap NextGen, a third CRISPR source (stages 31 to 33) -----------------------------------
np <- ext_params$nextgen
ng_freeze <- problem(assert_extension_freeze("nextgen", dir = ext_dir))
record("NextGen freeze intact", ng_freeze == "", ng_freeze)
ng_summaries <- c("nextgen_prepare_summary", "nextgen_scan_summary", "nextgen_combined_summary")
ng_stale <- ng_summaries[!vapply(ng_summaries, function(name) identical(ext_json(name)$parameters_sha256, ext_section_hash(c("shared", "nextgen"))), logical(1))]
ng_trial <- ng_summaries[vapply(ng_summaries, function(name) isTRUE(ext_json(name)$trial), logical(1))]
record("every NextGen stage ran in full with the current parameters", length(ng_stale) == 0 && length(ng_trial) == 0, paste(c(ng_stale, ng_trial), collapse = "; "))

# Nothing from before the NextGen stages changed. The extension parameter file gained a section; its
# earlier sections are held by the two freezes checked above.
ng_baseline_lines <- readLines(file.path(ext_dir, "nextgen_baseline_sha256.txt"))
ng_baseline <- data.table(sha256 = substr(ng_baseline_lines, 1, 64), path = substr(ng_baseline_lines, 67, nchar(ng_baseline_lines)))
ng_earlier <- ng_baseline[!grepl("^results/", path) & path != "scripts/loy_reanalysis_v1/extension_parameters.json"]
ng_current <- vapply(file.path(root_dir, ng_earlier$path), function(path) if (file.exists(path)) file_sha256(path) else NA_character_, character(1))
ng_changed <- ng_earlier$path[is.na(ng_current) | ng_current != ng_earlier$sha256]
record("every table and parameter file from before the NextGen stages is unchanged", nrow(ng_earlier) >= 650 && length(ng_changed) == 0,
  sprintf("%d files compared%s", nrow(ng_earlier), if (length(ng_changed) > 0) paste0("; changed: ", paste(head(ng_changed, 5), collapse = ", ")) else ""))
ng_baseline_report <- file.path(ext_dir, "nextgen_baseline_report.html")
ng_before <- report_sections(ng_baseline_report)
ng_untouched <- setdiff(names(ng_before), c(names(ng_before)[1], "Summary", "Limits", "Files"))
ng_added <- setdiff(names(after), names(ng_before))
record("report sections outside the NextGen section read as before", length(ng_untouched) >= 9 && all(ng_untouched %in% names(after)) && all(ng_before[ng_untouched] == after[ng_untouched]) && identical(ng_added, c("DepMap NextGen models", "Driver status inside cancer types")),
  sprintf("%d sections compared; added: %s", length(ng_untouched), paste(ng_added, collapse = ", ")))
ng_lost <- setdiff(report_blocks(ng_baseline_report), report_blocks(current_report))
record("every paragraph from before the NextGen stages is still in the report, apart from the three edited on purpose",
  length(ng_lost) == 3 && any(startsWith(ng_lost, "Other models.")) && any(startsWith(ng_lost, "In short")) && any(grepl("parameters.json holds every threshold", ng_lost, fixed = TRUE)), sprintf("%d earlier paragraphs differ", length(ng_lost)))
ng_typed <- prose_numbers(rmd_section(readLines(file.path(script_dir, "22_report.Rmd")), "DepMap NextGen models"))
record("every number in the NextGen section is computed from the tables", all(ng_typed %in% c("2026", "10", "05", "0.05", "0.10", "20,000")), sprintf("typed constants: %s", paste(unique(ng_typed), collapse = " ")))

# The cohort.
ng_prep <- ext_json("nextgen_prepare_summary")
ng_scan <- ext_json("nextgen_scan_summary")
ng_comb <- ext_json("nextgen_combined_summary")
ng_cohort <- ext_table("nextgen_cohort")
ng_frame <- ng_cohort[in_primary == TRUE][order(ModelID)]
ng_genes <- ext_table("nextgen_gene_annotation")
ng_pan <- ext_table("nextgen_pan_results")
ng_effect <- readRDS(file.path(ext_dir, "nextgen_effect.rds"))
ng_found <- unlist(ng_prep$counts_found)
ng_expected <- unlist(ng_prep$counts_expected)[names(ng_found)]
record("NextGen counts are those fixed before the screen", all(abs(ng_found - ng_expected) <= ext_params$shared$count_tolerance) && nrow(ng_frame) == ng_prep$primary && identical(unlist(np$expected_counts)[names(ng_found)], ng_expected),
  sprintf("%d screened, %d male, %d in the cohort (%d Y-negative, %d Y-positive); largest gap %d", ng_prep$screened, ng_prep$males, nrow(ng_frame), sum(ng_frame$y_negative == 1), sum(ng_frame$y_negative == 0), max(abs(ng_found - ng_expected))))
sanger_donors <- hcmi_donor(org_cohort$ModelID)
ng_damaging <- read_model_matrix(file.path(raw_dir, "studies", "depmap_nextgen", "next_gen_damaging.csv"))
ng_similarity <- genotype_similarity(ng_damaging, np$genotype_rare_fraction)[ng_frame$ModelID, ng_frame$ModelID]
ng_called <- ext_table("nextgen_y_calls")
record("only recorded-male cancers enter the NextGen cohort, one per donor, none shared with the Sanger cohort or the cell lines",
  all(ng_frame$Sex == "Male") && !any(ng_frame$OncotreePrimaryDisease %in% np$non_malignant) && !anyDuplicated(ng_frame$PatientID) && !anyNA(ng_frame$y_negative) &&
    !any(ng_frame$hcmi_donor %in% sanger_donors[!is.na(sanger_donors)]) && !any(ng_cohort$ModelID %in% cohort$ModelID[cohort$ModelID %in% rownames(effect)]) &&
    max(ng_similarity[upper.tri(ng_similarity)]) < np$genotype_same_donor && all(ng_frame$y_call == ng_called$y_call[match(ng_frame$ModelID, ng_called$ModelID)]),
  sprintf("%d models from %d donors; largest genotype similarity between two of them %.2f", nrow(ng_frame), uniqueN(ng_frame$PatientID), max(ng_similarity[upper.tri(ng_similarity)])))
record("no chrY gene and no gene with a missing score is in the NextGen scan", nrow(ng_pan) == sum(ng_genes$in_scan) && !any(ng_pan$chromosome %in% "Y") && !any(ng_pan$feature %in% ng_genes[is_y_gene == TRUE, feature]) && !anyNA(ng_effect[, ng_pan$feature]),
  sprintf("%d genes", nrow(ng_pan)))

# The scan: estimates with weighted lm, the controls, the bar from the shuffles, the rule for the four.
ng_weights <- ext_table("nextgen_weights")[, setNames(weight, ModelID)]
ng_paired <- ext_table("nextgen_paired_genes")
ng_refit <- rbindlist(lapply(seq_len(nrow(ng_paired)), function(i) {
  fit <- lm(ng_effect[ng_frame$ModelID, ng_paired$feature[i]] ~ factor(ng_frame$stratum) + ng_frame$y_negative, weights = ng_weights[ng_frame$ModelID])
  data.table(symbol = ng_paired$symbol[i], pipeline = ng_paired$estimate[i], lm = unname(tail(coef(fit), 1)))
}))
record("NextGen estimates agree with weighted lm", all(abs(ng_refit$pipeline - ng_refit$lm) < 1e-6), paste(sprintf("%s %.4f against %.4f", ng_refit$symbol, ng_refit$pipeline, ng_refit$lm), collapse = "; "))
ng_controls <- ext_table("nextgen_controls")
record("NextGen controls recover known biology", nrow(ng_controls) == 2 && all(ng_controls$pass) && all(sign(ng_controls$estimate) == ng_controls$expected_sign) && all(ng_controls$p < ng_controls$p_max),
  paste(sprintf("%s %.2f (p = %.2g)", ng_controls$gene, ng_controls$estimate, ng_controls$p), collapse = "; "))
ng_null <- ext_table("nextgen_calibration")[run != "observed"]
ng_bar <- sort(ng_null$pan_min_p)[floor(np$family_wise_level * np$null_runs + 1e-9)]
ng_discoveries <- ng_pan[discovery == TRUE]
record("the bar for a NextGen discovery is the one the shuffles give, and every call respects it",
  nrow(ng_null) == np$null_runs && abs(log10(ng_bar) - log10(ng_scan$null$pan_bar)) < 1e-5 && mean(ng_null$pan_min_p < ng_bar) < np$family_wise_level && isTRUE(all.equal(ng_pan$fdr, bh(ng_pan$p))) &&
    identical(ng_pan$associated, ng_pan$fdr < ext_params$shared$fdr & abs(ng_pan$estimate) >= np$min_effect_chronos) && identical(ng_pan$discovery, ng_pan$associated & ng_pan$p < ng_bar & !ng_pan$named) &&
    all(ng_discoveries$p < ng_bar),
  sprintf("%d shuffles; bar p < %.3g; %d runs below it; %d discoveries", nrow(ng_null), ng_bar, sum(ng_null$pan_min_p < ng_bar), nrow(ng_discoveries)))
ng_cell <- genes[match(ng_paired$symbol, symbol)]
ng_rule <- ifelse(sign(ng_paired$estimate) == sign(ng_cell$estimate) & p.adjust(ng_paired$p, "BH") < ext_params$shared$fdr, "replicated",
  ifelse(sign(ng_cell$estimate) * (ng_cell$estimate - ng_paired$estimate) > 0 & ng_paired$difference_p < ext_params$shared$fdr, "not replicated", "inconclusive"))
ng_chance <- ext_table("nextgen_replication_chance")
record("the paired-gene outcomes follow the rule fixed in advance, from the stored cell-line effects",
  identical(ng_rule, ng_paired$outcome) && all(ng_paired$cell_line_estimate == ng_cell$estimate) && all(ng_paired$estimate == ng_pan$estimate[match(ng_paired$feature, ng_pan$feature)]) &&
    identical(ng_paired$holds_under_permutation, ng_paired$outcome == "replicated" & p.adjust(ng_paired$permutation_p, "BH") < ext_params$shared$fdr) &&
    all(abs(ng_paired$chance_at_strictest - ng_chance$chance_at_strictest[match(ng_paired$symbol, ng_chance$symbol)]) < 1e-12),
  paste(sprintf("%s %s", ng_paired$symbol, ng_paired$outcome), collapse = "; "))

# Three sources: the stored p-values of every scan, and the pooled statistic by hand.
ng_three <- ext_table("nextgen_combined_pan")
ng_units <- c(ext_table("extension_power")[analysis == "Pan-cancer" & source == "Cell lines", unit], org_pan$unit[1], ng_pan$unit[1])
ng_by_hand <- function(gene) {
  row <- ng_three[symbol == gene]
  z <- c(row$cell_z, row$organoid_z, row$nextgen_z)
  c(by_hand = sum((z / ng_units)[!is.na(z)]) / sqrt(sum(1 / ng_units[!is.na(z)]^2)), stored = row$combined_z)
}
ng_hand <- vapply(c("DDX3X", "EIF1AX", "POT1"), ng_by_hand, numeric(2))
record("the three-source combination uses the stored p-values of every scan, and its statistic equals the weighted sum by hand",
  all(ng_three$cell_p == genes$p[match(ng_three$feature, genes$feature)]) && all(ng_three$organoid_p == org_pan$p[match(ng_three$feature, org_pan$feature)], na.rm = TRUE) &&
    all(ng_three$nextgen_p == ng_pan$p[match(ng_three$feature, ng_pan$feature)], na.rm = TRUE) && nrow(ng_three) == ext_table("nextgen_test_counts")$combination_pan_genes[1] &&
    all(ng_three$two_source_p == comb_pan$combined_p[match(ng_three$feature, comb_pan$feature)], na.rm = TRUE) && sum(!is.na(ng_three$two_source_p)) == nrow(comb_pan) &&
    max(abs(ng_hand["by_hand", ] - ng_hand["stored", ])) < 1e-8 && isTRUE(all.equal(ng_three$combined_fdr, bh(ng_three$combined_p))),
  sprintf("%d genes; largest gap of the statistic recomputed for three genes %.2g", nrow(ng_three), max(abs(ng_hand["by_hand", ] - ng_hand["stored", ]))))

# The pooled test, rebuilt in base R: ranks within source, lm with one term per stratum, and every
# arrangement of Y status counted for the NextGen models alone.
ng_pooled_all <- ext_table("nextgen_pooled_type")
ng_pooled_named <- ext_table("nextgen_pooled_named")
ng_sensitivity <- ext_table("nextgen_pooled_sensitivity")
pool_cell <- cohort[ModelID %in% rownames(effect) & !is.na(y_call) & OncotreeLineage == "Pancreas"]
pool_ng <- ng_frame[OncotreeLineage == "Pancreas"]
normal_scores <- function(v) qnorm((rank(v) - 0.5) / length(v))
pool_check <- rbindlist(lapply(ng_paired$feature, function(target) {
  scores <- c(normal_scores(effect[pool_cell$ModelID, target]), normal_scores(ng_effect[pool_ng$ModelID, target]))
  negative <- c(as.numeric(pool_cell$y_call == "Y_negative"), pool_ng$y_negative)
  stratum <- c(paste("cell", pool_cell$LibraryKY, pool_cell$LibraryHumagneCD), paste("nextgen", pool_ng$ScreenType))
  fit <- lm(scores ~ factor(stratum) + negative)
  # NextGen alone: every way to place its Y-negative models, and the two-sided p-value by counting.
  own <- normal_scores(ng_effect[pool_ng$ModelID, target])
  placed <- combn(nrow(pool_ng), sum(pool_ng$y_negative == 1), function(index) sum(own[index]))
  seen <- sum(own[pool_ng$y_negative == 1])
  data.table(feature = target, lm = unname(coef(fit)[["negative"]]), counted_p = min(1, 2 * min(mean(placed <= seen + 1e-9), mean(placed >= seen - 1e-9))), arrangements = length(placed))
}))
pool_stored <- ng_pooled_named[source == "Pooled"][match(pool_check$feature, feature)]
pool_alone <- ng_pooled_named[source == "NextGen"][match(pool_check$feature, feature)]
record("the pooled pancreas estimates agree with plain lm, and the NextGen p-values with a count of every arrangement",
  uniqueN(pool_ng$ScreenType) == 1 && all(abs(pool_check$lm - pool_stored$estimate_rank) < 1e-8) && all(abs(pool_check$counted_p - pool_alone$p) < 1e-12) && all(pool_check$arrangements == pool_alone$arrangements),
  paste(sprintf("%s %.3f against %.3f", gene_symbol(pool_check$feature), pool_check$lm, pool_stored$estimate_rank), collapse = "; "))
pool_sizes <- ng_pooled_named[source != "Pooled", .(n_negative = n_negative[1], n_positive = n_positive[1]), by = source]
record("no pooled p-value is below the floor of its test, and the floor is two over the number of arrangements",
  all(ng_pooled_all$p >= ng_pooled_all$floor - 1e-12) && all(ng_pooled_named$p >= ng_pooled_named$floor - 1e-12) && all(ng_sensitivity$p >= ng_sensitivity$floor - 1e-12) &&
    all(abs(ng_pooled_named[exact == TRUE, floor - 2 / arrangements]) < 1e-15) && all(ng_pooled_all$exact) && ng_pooled_all$arrangements[1] == prod(choose(pool_sizes$n_negative + pool_sizes$n_positive, pool_sizes$n_positive)),
  sprintf("%s arrangements, floor %.3g, smallest p %.3g over %d genes", format(ng_pooled_all$arrangements[1], big.mark = ","), ng_pooled_all$floor[1], min(ng_pooled_all$p), nrow(ng_pooled_all)))
# The cancer types to pool, recounted from the three cohorts.
recount <- rbind(
  cohort[ModelID %in% rownames(effect) & !is.na(y_call), .(negative = sum(y_call == "Y_negative"), positive = sum(y_call == "Y_positive")), by = .(type = OncotreeLineage)],
  org_frame[, .(negative = sum(y_negative == 1), positive = sum(y_negative == 0)), by = .(type = np$sanger_lineage$lineage[match(OncotreeLineage, np$sanger_lineage$type)])],
  ng_frame[, .(negative = sum(y_negative == 1), positive = sum(y_negative == 0)), by = .(type = OncotreeLineage)]
)
by_rule <- recount[, .(alone = any(negative >= np$pooled$min_group & positive >= np$pooled$min_group), together = sum(negative) >= np$pooled$min_group && sum(positive) >= np$pooled$min_group), by = type][!alone & together, sort(type)]
record("the cancer types pooled are those the rule gives when recounted", identical(by_rule, np$pooled$expected_types) && identical(sort(unique(ng_pooled_all$type)), by_rule), paste(by_rule, collapse = ", "))

# The inputs supplied on 2026-10-05.
ng_inputs <- as.data.table(ext_params$nextgen_inputs)
ng_download <- problem(for (i in seq_len(nrow(ng_inputs))) verify_download(file.path(root_dir, ng_inputs$file[i]), ng_inputs$checksum[i], ng_inputs$algorithm[i]))
record("the NextGen inputs match their recorded checksums", ng_download == "", if (ng_download == "") sprintf("%d files", nrow(ng_inputs)) else ng_download)

# Statements the NextGen section makes in words.
ng_types <- ext_table("nextgen_type_screen")
ng_follow <- ext_table("nextgen_hit_followup")[context == "All models"]
ng_three_hits <- ng_three[combined_hit == TRUE]
ng_new <- ng_three_hits[!symbol %in% ng_paired$symbol]
ng_other <- ng_pan[associated == TRUE & named == FALSE]
ng_three_follow <- ext_table("nextgen_combined_followup")
ng_within <- ext_table("nextgen_combined_within_type")
ng_shared <- ext_table("nextgen_shared_donors")
ng_rnai <- read_result("gene_pan_rnai")
of <- function(gene) ng_paired[symbol == gene]
ng_combined_hits <- ext_table("nextgen_combined_hits")
hit <- function(gene) ng_combined_hits[symbol == gene]
ng_extra <- ng_combined_hits[named == FALSE]
ng_stomach <- ng_frame[OncotreeLineage == "Esophagus/Stomach"]
replicated_share <- mean(ng_null$named_replicated > 0)
ng_statements <- c(
  "on the model's p-values all four pass the rule; on permutation p-values from a million shuffles DDX3X, EIF1AX and RPS4X replicate and ZFX does not" = all(ng_paired$outcome == "replicated") && all(ng_paired$estimate < 0) &&
    setequal(ng_paired[outcome_permutation == "replicated", symbol], c("DDX3X", "EIF1AX", "RPS4X")) && setequal(ng_paired[holds_under_permutation == TRUE, symbol], c("DDX3X", "EIF1AX", "RPS4X")) &&
    of("ZFX")$p < 0.05 && of("ZFX")$permutation_p > 0.05 && all(ng_paired$permutation_shuffles >= 1e6) && ng_follow[symbol == "ZFX" & check == "Y-status definition", all(estimate < 0)],
  "DDX3X also passes genome-wide, holds in every follow-up check, and is the only one of the four to meet the robustness rule" = ng_pan[symbol == "DDX3X", associated && passes_family_wise_bar] && identical(ng_paired[robust == TRUE, symbol], "DDX3X") &&
    ng_follow[symbol == "DDX3X" & !is.na(estimate) & check != "other source", all(estimate < 0)],
  "for EIF1AX and ZFX the Y-negative models score alike in the two sources and the Y-positive NextGen models are the more dependent" = all(ng_paired[symbol %in% c("EIF1AX", "ZFX"), abs(mean_negative - cell_line_mean_negative) < 0.05 & cell_line_mean_positive - mean_positive > 0.1]),
  "the KRAS control clears its gate with the weights and not without them" = ng_controls[gene == "KRAS", p < p_max && p_unweighted > p_max],
  "ZFX is not significant with screen quality or the genomic covariates in the model" = ng_follow[symbol == "ZFX" & variant %in% c("plus screen NNMD", "plus WGD, TP53, aneuploidy"), .N == 2 && all(p > 0.05)],
  "EIF1AX and ZFX are weaker than in the cell lines; DDX3X and RPS4X do not differ" = all(ng_paired[symbol %in% c("EIF1AX", "ZFX"), difference_p < 0.05 & abs(estimate) < abs(cell_line_estimate)]) && all(ng_paired[symbol %in% c("DDX3X", "RPS4X"), difference_p > 0.05]),
  "RPS4X had the lowest chance of replication" = ng_paired[which.min(chance_at_strictest), symbol] == "RPS4X",
  "shuffles call a paired gene replicated about twice as often as calibrated p-values would (about 2.5% of runs)" = replicated_share > 1.5 * 0.0255 && replicated_share < 3 * 0.0255,
  "two other genes pass FDR, both less needed in Y-negative models, neither clears the bar and nothing is called" = nrow(ng_other) == 2 && all(ng_other$estimate > 0) && !any(ng_other$passes_family_wise_bar) && nrow(ng_discoveries) == 0 && ng_other[which.min(p), p / ng_bar] < 1.1,
  "neither of the two shows anything in the cell lines, the Sanger organoids or RNAi, which hold more models" = all(genes[symbol %in% ng_other$symbol, p] > 0.05) && all(org_pan[symbol %in% ng_other$symbol, p] > 0.05) && all(ng_rnai[symbol %in% ng_other$symbol, p] > 0.05) &&
    nrow(org_frame) > nrow(ng_frame) && min(ng_rnai[symbol %in% ng_other$symbol, n]) > nrow(ng_frame),
  "shuffled Y status gives an FDR hit in more than 5% of runs and more genes below 1e-5 than expected" = mean(ng_null$pan_fdr_hits > 0) > 0.05 && mean(ng_null$pan_below_1e5) > 1e-5 * nrow(ng_pan),
  "few shuffles reach the observed inflation factor; the excess sits in the common-essential genes, leans positive, and falls with screen quality in the model" = ng_scan$null$runs_at_or_above_lambda < 0.05 * nrow(ng_null) &&
    ng_scan$essential$lambda > 1.5 && ng_scan$essential$positive > 0.6 && ng_scan$essential$runs_at_or_above_lambda < 0.05 * nrow(ng_null) && ng_scan$essential$lambda > ng_scan$essential$lambda_other_genes &&
    ng_scan$lambda_quality_adjusted < 1.1 && ng_pan[symbol == "POLR1G", common_essential && estimate > 0],
  "DDX3X survives the adjustments, its dependency contrast points the same way without being significant, and no gene is LOY-specific" = ng_pan[symbol == "DDX3X", survives_adjustment && survives_quality && survives_essential_shift && !dependency_contrast &&
    dependent_fraction_difference > 0 && dependent_fraction_p > 0.05 && dependent_fraction_p_without_quality < 0.05 && label == "LOY-associated"] && ng_scan$loy_specific == 0,
  "the continuous scan gives DDX3X, EIF1AX and OCEL1" = setequal(ext_table("nextgen_continuous")[fdr < ext_params$shared$fdr, symbol], c("DDX3X", "EIF1AX", "OCEL1")),
  "only oesophagus/stomach can be tested by cancer type, and nothing passes there, observed or shuffled" = identical(ng_scan$types, "Esophagus/Stomach") && ng_scan$type_specific == 0 && ng_scan$inside_only == 0 && all(ng_null$inside_fdr_hits == 0) && all(ng_null$type_specific == 0),
  "Y-negative screens are noisier, common-essential depletion does not differ, and Y loss goes with genome doubling" = ng_frame[y_negative == 1, median(ScreenNNMD)] > ng_frame[y_negative == 0, median(ScreenNNMD)] && ng_scan$quality_shift$permutation_p < 0.05 &&
    ng_scan$essential_shift$permutation_p > 0.05 && ng_frame[y_negative == 1, mean(WGD, na.rm = TRUE)] > ng_frame[y_negative == 0, mean(WGD, na.rm = TRUE)],
  "Y status is tied to screen format, and brain models are mostly Y-positive" = ng_frame[y_negative == 1, mean(dome)] > ng_frame[y_negative == 0, mean(dome)] + 0.3 && ng_frame[OncotreeLineage == "CNS/Brain", mean(y_negative)] < 0.5 &&
    all(ng_frame[OncotreeLineage == "CNS/Brain", ScreenType] %in% c("2DN", "3DN")),
  "recorded females show almost no Y signal next to the cell lines" = ng_prep$rna_thresholds$female_mean_max < 0.1 * jsonlite::fromJSON(file.path(derived_dir, "call_thresholds.json"))$female_mean_max,
  "six genes pass with three sources, three with two: the four paired genes, POT1 and ARHGAP26" = setequal(ng_three_hits$symbol, c("DDX3X", "EIF1AX", "ZFX", "RPS4X", "POT1", "ARHGAP26")) && sum(ng_three$two_source_hit, na.rm = TRUE) == 3,
  "DDX3X, EIF1AX and RPS4X are supported by every source that tests them and ZFX by two of three, the organoids being the one that does not" = all(ng_three_hits[symbol %in% c("DDX3X", "EIF1AX", "RPS4X"), evidence] == "supported by every source") &&
    ng_three_hits[symbol == "ZFX", evidence == "supported by 2 of 3 sources" && abs(organoid_z) < 1.96 && nextgen_z < 0] && ng_three_hits[symbol == "EIF1AX", is.na(organoid_p) && n_sources == 2] && !"EIF1AX" %in% org_pan$symbol,
  "POT1 is less needed in all three sources and ARHGAP26 more needed; neither is a hit in any single source, and both are below the single-source size in the cell lines" = ng_new[symbol == "POT1", cell_z > 0 && organoid_z > 0 && nextgen_z > 0] &&
    ng_new[symbol == "ARHGAP26", cell_z < 0 && nextgen_z < 0] && !any(ng_new$cell_hit) && !any(ng_new$organoid_hit) && !any(ng_new$nextgen_hit) && all(abs(ng_new$cell_estimate) < np$min_effect_chronos),
  "neither extra gene is a vulnerability of Y loss: POT1 is less needed with it, and ARHGAP26 is not a dependency" = hit("POT1")[, cell_estimate > 0 && nextgen_estimate > 0 && cell_dependent_negative < cell_dependent_positive] &&
    hit("ARHGAP26")[, cell_dependent_negative < 0.01 && cell_dependent_positive < 0.01 && nextgen_dependent_negative == 0 && nextgen_dependent_positive == 0 && cell_mean_negative > 0 && nextgen_mean_negative > 0],
  "neither extra gene was a hit with two sources, neither stays one with the genomic covariates in NextGen's model or with every source adjusted, and RNAi shows neither" = all(ng_extra$two_source_fdr > 0.05) &&
    all(ng_extra$fdr_nextgen_with_genomic_covariates > 0.05) && all(ng_extra$fdr_every_source_adjusted > 0.05) && all(ng_extra$rnai_p > 0.05) && all(ng_pan[symbol %in% ng_new$symbol, p_adjusted] > 0.05),
  "the follow-up of the combined hits uses a million joint shuffles, and every hit keeps its sign with RNA-only and DNA-confirmed calls" = all(ng_combined_hits$permutation_shuffles >= 1e6) &&
    ng_three_follow[check == "Y-status definition", all(sign(z) == sign(ng_three_follow[check == "all sources, least squares", z][match(symbol, ng_three_follow[check == "all sources, least squares", symbol])]))],
  "the four paired genes stay below the p-value they need at their rank, under the shuffles and with any one model left out" = all(ng_combined_hits[named == TRUE, permutation_below_cut & leave_one_out_below_cut]),
  "no NextGen model is among the cell lines, and NextGen carries the smallest share of the information" = !any(ng_cohort$ModelID %in% rownames(effect)) && ng_comb$weight_share$nextgen < ng_comb$weight_share$organoids && ng_comb$weight_share$organoids < ng_comb$weight_share$cell_lines,
  "no gene passes inside oesophagogastric cancer with three sources; DDX3X comes closest and gains on two sources" = !any(ng_within$combined_hit) && ng_within[which.min(combined_p), symbol == "DDX3X" && combined_fdr < two_source_fdr],
  "NextGen has too few bowel models for the gate" = ng_frame[OncotreeLineage == "Bowel", min(sum(y_negative == 1), sum(y_negative == 0))] < ext_params$shared$min_cell,
  "inside pancreas DDX3X, EIF1AX and ZFX pass among the four and RPS4X does not; all Y-positive pancreatic cell lines sit in one library" = setequal(ng_pooled_named[source == "Pooled" & named_fdr < ext_params$shared$fdr, symbol], c("DDX3X", "EIF1AX", "ZFX")) &&
    all(ng_pooled_named[source == "Pooled" & symbol != "RPS4X", estimate_chronos] < 0) && ng_pooled_named[source == "Pooled" & symbol == "RPS4X", p > 0.05] && pool_cell[y_call == "Y_positive", uniqueN(paste(LibraryKY, LibraryHumagneCD))] == 1,
  "a single gene cannot be called genome-wide in pancreas, far fewer genes sit at the smallest p-value than BH would need, and none is called" = ng_comb$pooled$Pancreas$floor > ng_comb$pooled$Pancreas$bh_critical_p &&
    ng_comb$pooled$Pancreas$min_family_wise_p > 0.05 && is.na(ext_table("nextgen_power")[analysis == "Pancreas, pooled", genome_wide]) &&
    sum(ng_pooled_all$p <= ng_pooled_all$floor * (1 + 1e-9)) < ceiling(ng_pooled_all$floor[1] * nrow(ng_pooled_all) / ext_params$shared$fdr) && min(ng_pooled_all$fdr) > 0.05,
  "the cell lines alone had been tested inside pancreas before without a significant result, and NextGen alone gives p below 0.05 for the three genes" = read_result("robustness_gene_hits")[check == "within lineage" & variant == "Pancreas" & symbol %in% ng_paired$symbol, .N == 4 && all(p > 0.05)] &&
    all(ng_pooled_named[source == "NextGen" & symbol %in% c("DDX3X", "EIF1AX", "ZFX"), p] < 0.05) && all(abs(ng_pooled_named[source == "Cell lines"][match(ng_paired$symbol, symbol), estimate_chronos] -
      read_result("robustness_gene_hits")[check == "within lineage" & variant == "Pancreas"][match(ng_paired$symbol, symbol), estimate]) < 1e-6),
  "all four paired genes have p below 0.05 in each lineage of the secondary table" = ng_sensitivity[, .N == 12 && all(p < 0.05) && all(estimate_rank < 0)],
  "every shared model carries the same name at both institutes; Y calls agree in all but one of those called at both, and the exception is one model, Y-positive at Sanger and Y-negative at the Broad" = all(ng_shared$same_name) &&
    sum(!ng_shared$calls_agree, na.rm = TRUE) == 1 && ng_shared[calls_agree %in% FALSE, same_name && nextgen_call == "Y_negative" && sanger_call == "Y_positive" && sanger_tier == "concordant" && nextgen_tier == "concordant" && sanger_in_cohort],
  "two screens of one model share little of their profile, though more than unrelated models do" = ng_comb$shared_donors$profile_median > ng_comb$shared_donors$unrelated_median && ng_comb$shared_donors$profile_median < 0.1 &&
    ng_comb$shared_donors$above_unrelated_q95 > 0.05 * ng_comb$shared_donors$compared,
  "in oesophagus/stomach every Y-positive model comes from a treated patient and most Y-negative ones do not; the Y-positive pancreatic models come from one laboratory" = all(ng_stomach[y_negative == 0, treatment_status] == "Treated") &&
    ng_stomach[y_negative == 1, mean(treatment_status %in% c("Treated", "Post-treatment"))] < 0.5 && pool_ng[y_negative == 0, uniqueN(source_lab)] == 1,
  "the four paired genes keep their sign without the Y-negative models whose X copy number is female-like" = ng_follow[symbol %in% ng_paired$symbol & variant == "x_checked", .N == 4 && all(estimate < 0)] && sum(ng_frame$possible_female) > 0,
  "the thresholds typed in the text are the parameters" = np$min_effect_chronos == 0.10 && ext_params$shared$fdr == 0.05 && ext_params$shared$followup_permutations == 20000,
  "with three sources the detectable effect falls a little, and NextGen alone finds only large effects" = ext_table("nextgen_power")[, genome_wide[analysis == "All cancer types" & source == "Three sources"] < genome_wide[analysis == "All cancer types" & source == "Two sources"] &&
    genome_wide[analysis == "Oesophagogastric adenocarcinoma" & source == "Three sources"] < genome_wide[analysis == "Oesophagogastric adenocarcinoma" & source == "Two sources"] && genome_wide[analysis == "All cancer types" & source == "NextGen"] > 1.5]
)
record("statements the NextGen section makes in words hold in the tables", !anyNA(ng_statements) && all(ng_statements), paste(names(ng_statements)[is.na(ng_statements) | !ng_statements], collapse = "; "))
ng_text <- after[["DepMap NextGen models"]]
ng_required <- c("named in advance", "frozen before", "not that it is specific to pancreas", "no gene is called from it", "about twice as often")
record("the NextGen section states what was fixed in advance and what its results cannot claim", all(vapply(ng_required, grepl, logical(1), x = ng_text, fixed = TRUE)), paste(ng_required[!vapply(ng_required, grepl, logical(1), x = ng_text, fixed = TRUE)], collapse = "; "))

# Findings of the fresh review of 2026-10-05: the tables behind each corrected sentence, and the wording.
ng_donor_table <- ext_table("nextgen_donor_groups")
at_sanger <- ng_donor_table[kind == "donor screened at Sanger"]
record("a model shared with the Sanger cohort is paired with the Sanger organoid of the same name, and the model kept is one of the Sanger cohort",
  nrow(ng_shared) > 0 && all(ng_shared$name == ng_shared$sanger_model) && ng_cohort[shared_with_sanger == TRUE, all(name == sanger_model)] && all(at_sanger$other == at_sanger$name) &&
    all(unlist(strsplit(at_sanger[!is.na(kept), kept], "; ", fixed = TRUE)) %in% org_frame$ModelID) && ng_cohort[shared_with_sanger == FALSE, all(is.na(sanger_model))],
  sprintf("%d shared models; %d with a Sanger organoid in its cohort", nrow(at_sanger), sum(!is.na(at_sanger$kept))))
# The first NextGen freeze was replaced after the review. What it held is proven unchanged: every frozen file
# but the cohort has its first hash; the cohort has it again once the pairing column is put back as the first
# join made it; and the parameters have theirs once the two settings added with the fix are taken out.
archived_freezes <- sort(list.files(ext_dir, pattern = "^extension_freeze_nextgen_archived_.*\\.json$", full.names = TRUE))
freeze_proof <- problem({
  stopifnot(length(archived_freezes) == 1)
  first <- jsonlite::fromJSON(archived_freezes[1])
  now <- jsonlite::fromJSON(file.path(ext_dir, "extension_freeze_nextgen.json"))
  others <- setdiff(names(first$sha256), "cohort")
  stopifnot(setequal(names(first$sha256), names(now$sha256)), all(unlist(first$sha256[others]) == unlist(now$sha256[others])), first$sha256$cohort != now$sha256$cohort)
  as_first <- copy(ng_cohort)
  first_join <- sanger_donors
  as_first[, sanger_model := org_cohort$ModelID[match(hcmi_donor, first_join)]]
  as_first[shared_with_sanger == FALSE, sanger_model := NA_character_]
  rewritten <- tempfile(fileext = ".csv")
  fwrite(as_first, rewritten)
  stopifnot(identical(file_sha256(rewritten), first$sha256$cohort))
  unlink(rewritten)
  kept_settings <- ext_params[c("shared", "nextgen")]
  kept_settings$nextgen[c("paired_permutations", "combined_permutations")] <- NULL
  as_text <- as.character(jsonlite::toJSON(kept_settings, auto_unbox = TRUE, digits = NA, null = "null"))
  stopifnot(identical(digest::digest(as_text, algo = "sha256", serialize = FALSE), first$parameters_sha256))
})
record("the replaced NextGen freeze differs from the first only in the Sanger pairing column and two added settings", freeze_proof == "", freeze_proof)
record("pooled estimates in Chronos units are given only where every pooled source is on that scale",
  ng_sensitivity[grepl("Sanger", sources), .N > 0 && all(is.na(estimate_chronos))] && ng_sensitivity[!grepl("Sanger", sources), .N > 0 && !anyNA(estimate_chronos)] && !anyNA(ng_pooled_named$estimate_chronos),
  sprintf("%d rows with a Sanger source, %d without", ng_sensitivity[grepl("Sanger", sources), .N], ng_sensitivity[!grepl("Sanger", sources), .N]))
# A dependency call that varies only inside strata holding one Y group says nothing about Y once screen quality
# is out of the model; a p-value there would be rounding noise and differ between processors.
ng_dependency <- readRDS(file.path(ext_dir, "nextgen_dependency.rds"))
ng_calls_made <- (ng_dependency[ng_frame$ModelID, ng_pan$feature] >= params$gene_scan$dependency_probability_cut) * 1
ng_strata <- model.matrix(~ factor(ng_frame$stratum))
ng_exposure_left <- qr.resid(qr(ng_strata), ng_frame$y_negative)
ng_without_contrast <- colSums(abs(qr.resid(qr(ng_strata), ng_calls_made) * ng_exposure_left) > 1e-8) == 0
record("a gene whose dependency calls leave no Y contrast without screen quality carries no p-value there",
  "dependent_fraction_p_without_quality" %in% names(ng_pan) && all(is.na(ng_pan$dependent_fraction_p_without_quality[ng_without_contrast])) && !anyNA(ng_pan$dependent_fraction_p_without_quality[!ng_without_contrast]),
  sprintf("%d of %d genes without contrast; %d of them with a p-value", sum(ng_without_contrast), length(ng_without_contrast), sum(!is.na(ng_pan$dependent_fraction_p_without_quality[ng_without_contrast]))))
ng_retracted <- c("ll four paired genes replicate", "confirm all four", "All four paired genes are now seen", "first test of the paired genes", "tested for the first time", "Both pass the follow-up",
  "two different tumour samples", "inside the range of the shuffles", "No genome-wide call is possible", "candidates found by pooling", "not established vulnerabilities", "sizes compare directly")
ng_review_required <- c("ZFX is not confirmed", "added after the first full run", "designed after the two genes were seen", "is not a dependency", "same model name at both institutes", "belongs to the culture",
  "had been tested inside pancreas before", "common-essential genes", "A Cas12a guide that also cuts the Y copy", "second look")
record("report wording follows the fresh review of the NextGen stages", all(has(ng_review_required)) && !any(has(ng_retracted)),
  paste(c(sprintf("missing: %s", ng_review_required[!has(ng_review_required)]), sprintf("still present: %s", ng_retracted[has(ng_retracted)])), collapse = "; "))

# --- driver status inside cancer types (stages 34 to 36) ----------------------------------------
dp <- ext_params$driver_type
sp <- ext_params$shared
dt_freeze <- problem(assert_extension_freeze("driver_type", dir = ext_dir))
record("driver-by-type freeze intact", dt_freeze == "", dt_freeze)
dt_summaries <- c("driver_type_prepare_summary", "driver_type_tests_summary", "driver_type_power_summary")
dt_stale <- dt_summaries[!vapply(dt_summaries, function(name) identical(ext_json(name)$parameters_sha256, ext_section_hash(c("shared", "driver_type"))), logical(1))]
dt_trial <- dt_summaries[vapply(dt_summaries, function(name) isTRUE(ext_json(name)$trial), logical(1))]
record("every driver-by-type stage ran in full with the current parameters", length(dt_stale) == 0 && length(dt_trial) == 0, paste(c(dt_stale, dt_trial), collapse = "; "))
record("the power stage of the driver-by-type layer used the frozen power parameters", identical(ext_json("driver_type_power_summary")$power_parameters_sha256, file_sha256(file.path(script_dir, "power_parameters.json"))), "")

# Nothing from before these stages changed. The extension parameter file gained a section; its
# earlier sections are held by the three freezes checked above.
dt_baseline_lines <- readLines(file.path(ext_dir, "driver_type_baseline_sha256.txt"))
dt_baseline <- data.table(sha256 = substr(dt_baseline_lines, 1, 64), path = substr(dt_baseline_lines, 67, nchar(dt_baseline_lines)))
dt_earlier <- dt_baseline[!grepl("^results/", path) & path != "scripts/loy_reanalysis_v1/extension_parameters.json"]
dt_current <- vapply(file.path(root_dir, dt_earlier$path), function(path) if (file.exists(path)) file_sha256(path) else NA_character_, character(1))
dt_changed <- dt_earlier$path[is.na(dt_current) | dt_current != dt_earlier$sha256]
record("every table and parameter file from before the driver-by-type stages is unchanged", nrow(dt_earlier) >= 775 && length(dt_changed) == 0,
  sprintf("%d files compared%s", nrow(dt_earlier), if (length(dt_changed) > 0) paste0("; changed: ", paste(head(dt_changed, 5), collapse = ", ")) else ""))
dt_baseline_report <- file.path(ext_dir, "driver_type_baseline_report.html")
dt_before <- report_sections(dt_baseline_report)
dt_untouched <- setdiff(names(dt_before), c(names(dt_before)[1], "Summary", "Limits", "Files"))
dt_added <- setdiff(names(after), names(dt_before))
record("report sections outside the driver-by-type section read as before", length(dt_untouched) >= 10 && all(dt_untouched %in% names(after)) && all(dt_before[dt_untouched] == after[dt_untouched]) && identical(dt_added, "Driver status inside cancer types"),
  sprintf("%d sections compared; added: %s", length(dt_untouched), paste(dt_added, collapse = ", ")))
dt_lost <- setdiff(report_blocks(dt_baseline_report), report_blocks(current_report))
record("every paragraph from before the driver-by-type stages is still in the report, apart from the two edited on purpose",
  length(dt_lost) == 2 && any(startsWith(dt_lost, "In short")) && any(grepl("parameters.json holds every threshold", dt_lost, fixed = TRUE)), sprintf("%d earlier paragraphs differ", length(dt_lost)))
dt_typed <- prose_numbers(rmd_section(readLines(file.path(script_dir, "22_report.Rmd")), "Driver status inside cancer types"))
record("every number in the driver-by-type section is computed from the tables", all(dt_typed %in% c("2026", "10", "06", "0.05")), sprintf("typed constants: %s", paste(unique(dt_typed), collapse = " ")))

# Statuses, census and splits.
dt_prep <- ext_json("driver_type_prepare_summary")
dt_tests <- ext_json("driver_type_tests_summary")
dt_power <- ext_json("driver_type_power_summary")
dt_status <- ext_table("driver_type_status")
dt_events <- ext_table("driver_type_events")
dt_census <- ext_table("driver_type_census")
dt_splits <- ext_table("driver_type_splits")
dt_lines <- ext_table("driver_type_split_lines")
dt_frame <- cohort[match(dt_status$ModelID, ModelID)]
dt_frame[, exposed := as.numeric(y_call == "Y_negative")]
dt_matrix <- as.matrix(dt_status[, !"ModelID"])
record("driver-by-type lines are the male CRISPR lines of the frozen cohort with a Y call", nrow(dt_status) == dp$expected_counts$lines && !anyNA(dt_frame$ModelID) && !anyNA(dt_frame$exposed) && sum(dt_frame$exposed) == dp$expected_counts$y_negative,
  sprintf("%d lines, %d Y-negative", nrow(dt_status), sum(dt_frame$exposed)))
# One status per gene: a line is altered exactly where the event table lists it.
dt_altered <- which(dt_matrix == 1, arr.ind = TRUE)
record("a line is altered for a gene exactly where an event is listed", nrow(dt_altered) == nrow(dt_events) &&
  setequal(paste(dt_status$ModelID[dt_altered[, 1]], colnames(dt_matrix)[dt_altered[, 2]]), paste(dt_events$ModelID, dt_events$gene)) && dt_events[, all(mutation | amplification | deletion)],
  sprintf("%d altered pairs, %d events", nrow(dt_altered), nrow(dt_events)))
# The same definitions as the frozen context matrix, kind by kind.
dt_context <- read_result("context_matrix")
dt_stored <- read_result("context_definitions")[kept == TRUE & class %in% c("driver", "amplification", "deletion")]
dt_kind <- c(driver = "mutation", amplification = "amplification", deletion = "deletion")
dt_same <- vapply(seq_len(nrow(dt_stored)), function(i) {
  listed <- dt_events[gene == dt_stored$label[i] & get(dt_kind[[dt_stored$class[i]]]) == TRUE, ModelID]
  identical(as.numeric(dt_status$ModelID %in% listed), as.numeric(dt_context[[dt_stored$context_id[i]]][match(dt_status$ModelID, dt_context$ModelID)]))
}, logical(1))
record("the rebuilt statuses equal every stored alteration context, kind by kind", nrow(dt_stored) >= 85 && all(dt_same), sprintf("%d of %d contexts equal", sum(dt_same), nrow(dt_stored)))
dt_levels <- unlist(dp$levels)
# Counted here from the long table, without the helper that stage 34 used.
dt_long <- melt(cbind(dt_frame[, c("exposed", unname(dt_levels)), with = FALSE], as.data.table(dt_matrix)), id.vars = c("exposed", unname(dt_levels)), variable.name = "gene", value.name = "status", variable.factor = FALSE)
dt_recount <- rbindlist(lapply(names(dt_levels), function(level) {
  out <- dt_long[!is.na(status) & !is.na(get(dt_levels[[level]])) & get(dt_levels[[level]]) != "", .(
    neg_alt = sum(exposed == 1 & status == 1), pos_alt = sum(exposed == 0 & status == 1), neg_wt = sum(exposed == 1 & status == 0), pos_wt = sum(exposed == 0 & status == 0)
  ), by = c(dt_levels[[level]], "gene")]
  setnames(out, dt_levels[[level]], "type")
  out <- out[neg_alt + pos_alt > 0]
  out[, label := fifelse(pmin(neg_alt, pos_alt, neg_wt, pos_wt) >= sp$min_cell, "specificity testable", fifelse(pmin(neg_alt, pos_alt) >= sp$min_cell, "Y effect among altered only", fifelse(pmin(neg_alt, pos_alt) == 0, "one Y group missing", "below the gate")))]
  out[, level := level]
}))
rm(dt_long)
dt_cells <- dt_census[source == "Cell lines"]
dt_joined <- merge(dt_cells, dt_recount, by = c("level", "type", "gene"), suffixes = c("", "_again"))
record("the census recounted from the status table equals the stored one", nrow(dt_joined) == nrow(dt_cells) && nrow(dt_recount) == nrow(dt_cells) &&
  dt_joined[, all(neg_alt == neg_alt_again & pos_alt == pos_alt_again & neg_wt == neg_wt_again & pos_wt == pos_wt_again & label == label_again)], sprintf("%d pairs", nrow(dt_cells)))
dt_found <- dt_prep$counts_found
record("driver-by-type counts match those fixed before the stages ran", all(vapply(names(dp$expected_counts), function(key) identical(as.numeric(dt_found[[key]]), as.numeric(dp$expected_counts[[key]])), logical(1))),
  sprintf("%d splits, %d with a specificity test", dt_found$distinct_splits, dt_found$distinct_specificity))
dt_fitted <- dt_splits[fitted == TRUE]
dt_membership <- dt_lines[, .(neg_alt = sum(exposed == 1 & ctx == 1), pos_alt = sum(exposed == 0 & ctx == 1), neg_wt = sum(exposed == 1 & ctx == 0), pos_wt = sum(exposed == 0 & ctx == 0)), by = split_id]
dt_membership <- merge(dt_fitted[, .(split_id, gene, neg_alt, pos_alt, neg_wt, pos_wt, label)], dt_membership, by = "split_id", suffixes = c("", "_lines"))
dt_line_status <- dt_matrix[cbind(match(dt_lines$ModelID, dt_status$ModelID), match(dt_fitted$gene[match(dt_lines$split_id, dt_fitted$split_id)], colnames(dt_matrix)))]
record("every fitted split reaches the gate, and its stored lines carry the frozen status and Y call", nrow(dt_membership) == nrow(dt_fitted) &&
  dt_membership[, all(pmin(neg_alt, pos_alt) >= sp$min_cell & neg_alt == neg_alt_lines & pos_alt == pos_alt_lines & neg_wt == neg_wt_lines & pos_wt == pos_wt_lines)] &&
  identical(as.numeric(dt_line_status), as.numeric(dt_lines$ctx)) && identical(as.numeric(dt_lines$exposed), dt_frame$exposed[match(dt_lines$ModelID, dt_frame$ModelID)]),
  sprintf("%d splits, %d lines listed", nrow(dt_fitted), nrow(dt_lines)))

# The tests.
dt_screen <- as.data.table(readRDS(file.path(ext_dir, "driver_type_screen.rds")))
dt_controls <- ext_table("driver_type_controls")
record("driver-by-type controls point the expected way and are significant", nrow(dt_controls) == length(dp$controls$gene) && dt_controls[, all(pass & sign(estimate) == expected_direction & p < dp$controls$p_max)],
  paste(sprintf("%s %+.2f (p %.1e)", dt_controls$control, dt_controls$estimate, dt_controls$p), collapse = "; "))
dt_four <- dt_fitted[label == "specificity testable", split_id]
record("no specificity result exists for a split without five lines in all four cells", dt_screen[!split_id %in% dt_four, all(is.na(interaction_p))] && dt_screen[split_id %in% dt_four, all(is.finite(interaction_p))] &&
  setequal(unique(dt_screen$split_id), dt_fitted$split_id), sprintf("%d splits with the test, %d without", length(dt_four), nrow(dt_fitted) - length(dt_four)))
dt_adjust <- function(p) { out <- rep(NA_real_, length(p)); out[is.finite(p)] <- p.adjust(p[is.finite(p)], method = "BH"); out }
dt_rule <- dt_screen[, .(inside_fdr_again = dt_adjust(inside_p), interaction_fdr_again = dt_adjust(interaction_p), inside_fdr, interaction_fdr, specific, inside_only)]
dt_rule[, specific_again := !is.na(interaction_fdr_again) & interaction_fdr_again < sp$fdr & inside_fdr_again < sp$fdr]
dt_hits <- ext_table("driver_type_hits")
record("the rule of the driver-by-type layer recomputed from the stored p-values gives the stored calls", dt_rule[, all(abs(inside_fdr - inside_fdr_again) < 1e-12) && all(abs(interaction_fdr - interaction_fdr_again) < 1e-12, na.rm = TRUE) &&
  identical(specific, specific_again) && identical(inside_only, inside_fdr_again < sp$fdr & !specific_again)] && nrow(dt_hits) == dt_rule[, sum(specific | inside_only)] && dt_tests$calls == sum(dt_rule$specific) && dt_tests$inside_only == sum(dt_rule$inside_only),
  sprintf("%d calls, %d among altered lines only, %d tests", sum(dt_rule$specific), sum(dt_rule$inside_only), nrow(dt_rule)))
# One estimate of each quantity against plain lm, on the raw Chronos scale.
dt_check_split <- dt_fitted[label == "specificity testable"][order(-pmin(neg_alt, pos_alt, neg_wt, pos_wt))][1, split_id]
dt_check_lines <- dt_lines[split_id == dt_check_split]
dt_check_frame <- cbind(dt_frame[match(dt_check_lines$ModelID, ModelID), .(ModelID, LibraryKY, LibraryHumagneCD)], dt_check_lines[, .(ctx, exposed_in = exposed * ctx, exposed_out = exposed * (1 - ctx))])
dt_check_row <- dt_screen[split_id == dt_check_split & symbol == "DDX3X"]
dt_effect <- readRDS(file.path(matrix_dir, "crispr_effect.rds"))[dt_check_frame$ModelID, dt_check_row$feature]
dt_lm <- coef(lm(dt_effect ~ LibraryKY + LibraryHumagneCD + ctx + exposed_in + exposed_out, data = dt_check_frame))
dt_rank <- rank_normalise(readRDS(file.path(matrix_dir, "crispr_effect.rds"))[dt_status$ModelID, dt_check_row$feature])[match(dt_check_frame$ModelID, dt_status$ModelID)]
dt_lm_rank <- summary(lm(dt_rank ~ LibraryKY + LibraryHumagneCD + ctx + exposed_in + exposed_out, data = dt_check_frame))$coefficients
record("the least-squares p-value stored beside a driver-by-type test is that of plain lm on the ranked scores", abs(dt_check_row$inside_p_ols / dt_lm_rank["exposed_in", "Pr(>|t|)"] - 1) < 1e-8 &&
  abs(dt_check_row$inside_rank - dt_lm_rank["exposed_in", "Estimate"]) < 1e-8, sprintf("%s, DDX3X: p %.2e by least squares, %.2e by the layer's test", dt_check_split, dt_check_row$inside_p_ols, dt_check_row$inside_p))
record("a driver-by-type estimate equals plain lm for each of its three quantities", nrow(dt_check_row) == 1 && abs(dt_check_row$inside - dt_lm[["exposed_in"]]) < 1e-8 && abs(dt_check_row$outside - dt_lm[["exposed_out"]]) < 1e-8 &&
  abs(dt_check_row$interaction - (dt_lm[["exposed_in"]] - dt_lm[["exposed_out"]])) < 1e-8, sprintf("%s, DDX3X: among altered %+.3f, among unaltered %+.3f", dt_check_split, dt_lm[["exposed_in"]], dt_lm[["exposed_out"]]))
dt_null <- ext_table("driver_type_calibration")[run != "observed"]
dt_null_cells <- dt_null[source == "Cell lines"]
# Bounds written after the shuffles had been read (4.8% below 0.05, no run with a call); they guard a rerun, they do not judge the first run.
record("shuffled Y status gives nominal p-values and no more calls than chance in the driver-by-type layer", nrow(dt_null_cells) == dp$null_runs && dt_null[source == "Sanger organoids", .N] == dp$null_runs &&
  dt_null_cells[, mean(inside_below_0.05) > 0.04 && mean(inside_below_0.05) < 0.06 && mean(interaction_below_0.05) > 0.04 && mean(interaction_below_0.05) < 0.06 && mean(specific > 0) <= 0.05 && mean(inside_only > 0) <= 0.10],
  sprintf("%d runs; %.3f below 0.05; %d with a call, %d with a hit among altered lines", nrow(dt_null_cells), mean(dt_null_cells$inside_below_0.05), sum(dt_null_cells$specific > 0), sum(dt_null_cells$inside_only > 0)))
dt_paired <- ext_table("driver_type_paired_genes")
record("the paired genes are judged as a family of their own inside the splits", nrow(dt_paired) == 4 * nrow(dt_fitted) && dt_paired[, all(abs(named_fdr - bh(inside_p)) < 1e-12) && identical(shown, named_fdr < sp$fdr & inside_rank < 0)],
  sprintf("%d of %d shown", sum(dt_paired$shown), nrow(dt_paired)))
dt_follow <- ext_table("driver_type_hit_followup")
dt_followed <- unique(dt_follow[, .(split_id, feature, robust)])
record("every hit and every paired gene shown inside a split has its follow-up", nrow(dt_followed) == nrow(dt_hits) + dt_paired[shown == TRUE, .N] && dt_tests$followed_up == nrow(dt_followed) &&
  dt_follow[check == "permutation within status", .N] == nrow(dt_followed), sprintf("%d followed up, %d robust", nrow(dt_followed), sum(dt_followed$robust)))
dt_org_splits <- ext_table("driver_type_organoid_splits")
dt_org <- ext_table("driver_type_organoid_screen")
dt_comb <- ext_table("driver_type_combined")
dt_pairs <- ext_table("driver_type_pairs")
record("organoid splits and the combined pairs are those frozen", setequal(unique(dt_org$split_id), dt_org_splits$split_id) && dt_org[, all(abs(inside_fdr - bh(inside_p)) < 1e-12)] &&
  setequal(unique(paste(dt_comb$pair, dt_comb$gene)), paste(dt_pairs$pair, dt_pairs$gene)) && dt_comb[, all(abs(combined_fdr - bh(combined_p)) < 1e-12)],
  sprintf("%d organoid splits, %d pairs, %d combined tests", nrow(dt_org_splits), nrow(dt_pairs), nrow(dt_comb)))
dt_z <- dt_comb[1, c(signed_z(cell_inside_p, cell_inside_rank), signed_z(organoid_inside_p, organoid_inside_rank))]
dt_by_hand <- dt_comb[1, sum(dt_z / c(cell_unit, organoid_unit)) / sqrt(1 / cell_unit^2 + 1 / organoid_unit^2)]
record("a combined driver-by-type statistic equals the weighted sum computed by hand", abs(dt_by_hand - dt_comb$combined_z[1]) < 1e-10 && abs(2 * pnorm(-abs(dt_by_hand)) / dt_comb$combined_p[1] - 1) < 1e-8, sprintf("z = %.3f", dt_by_hand))

# What each split could detect.
dt_info <- ext_table("driver_type_power_splits")
dt_limits <- ext_table("driver_type_limits")
dt_enlarged <- ext_table("driver_type_enlarged_limits")
dt_needed <- ext_table("driver_type_lines_needed")
record("the power stage reproduced the stored p-values of every split, by refitting and through its simulation path", setequal(dt_info$split_id, dt_fitted$split_id) && max(dt_info$reproduce_gap) <= 1e-8 &&
  "path_gap" %in% names(dt_info) && all(is.finite(dt_info$path_gap)) && max(dt_info$path_gap) <= 0.02 && abs(dt_power$bars$family_inside / dt_tests$bars$family_inside - 1) < 1e-10,
  sprintf("%d splits; largest gap in log10 p %.1e by refitting, %.1e through the simulation path", nrow(dt_info), max(dt_info$reproduce_gap), if ("path_gap" %in% names(dt_info)) max(dt_info$path_gap) else NA_real_))
dt_floor_again <- dt_info[, (1 + (neg_alt == pos_alt)) / choose(neg_alt + pos_alt, neg_alt) + 1 / (sp$followup_permutations + 1)]
dt_robust_rows <- unique(dt_follow[robust == TRUE, split_id])
# The floor is an average over the shuffles, so a single p-value may fall below it; none may fall below half of the
# most favourable case (the observed arrangement alone), which no draw of 20,000 shuffles does for these counts.
dt_plain_floor <- dt_info[, 1 / choose(neg_alt + pos_alt, neg_alt) + 1 / (sp$followup_permutations + 1)]
record("a robust result exists only in splits whose shuffles can give a p-value below the robustness threshold", "robust_reachable" %in% names(dt_info) && identical(dt_info$robust_reachable, dt_floor_again < sp$robust_permutation_p) &&
  identical(dt_info$robust_reachable, dt_plain_floor < sp$robust_permutation_p) && !any(dt_robust_rows %in% dt_info[robust_reachable == FALSE, split_id]) &&
  dt_follow[check == "permutation within status", all(p >= 0.5 * dt_plain_floor[match(split_id, dt_info$split_id)])],
  sprintf("%d splits where no robust result is possible; %d splits with a robust result", sum(!dt_info$robust_reachable), length(dt_robust_rows)))
dt_named_limits <- dt_limits[unit == "sd" & target == 0.8 & bar == "candidate" & reachable == TRUE]
# Bounds fixed before the power stage was read: at the bar of a named gene the simulation and normal theory should agree within a factor of 1.5.
dt_ratio <- dt_named_limits[, limit / best_case]
record("at the bar of a named gene the simulated limits agree with normal theory on the design", nrow(dt_named_limits) >= 0.9 * nrow(dt_fitted) && mean(dt_ratio > 1 / 1.5 & dt_ratio < 1.5) >= 0.9,
  sprintf("%d splits; ratio simulated to normal theory: median %.2f, range %.2f to %.2f", nrow(dt_named_limits), median(dt_ratio), min(dt_ratio), max(dt_ratio)))
dt_four_bars <- c("specificity", "named_difference")
record("no limit that needs the difference is given for a split that cannot test it", dt_limits[bar %in% dt_four_bars, all(split_id %in% dt_four)] && dt_needed[bar %in% dt_four_bars, all(split_id %in% dt_four)] &&
  all(vapply(dt_four_bars, function(b) uniqueN(dt_limits[bar == b, split_id]) == length(dt_four), logical(1))), sprintf("%d splits", uniqueN(dt_limits[bar == "named_difference", split_id])))
record("limits exist for a least-squares test on the Chronos scale at the layer's bar", dt_limits[bar == "family_raw" & unit == "sd" & target == 0.8, .N == nrow(dt_fitted)], sprintf("%d of %d splits reach 80%%", dt_limits[bar == "family_raw" & unit == "sd" & target == 0.8, sum(reachable)], nrow(dt_fitted)))
dt_on_grid <- dt_needed[below_grid == FALSE & above_grid == FALSE]
dt_needed_again <- dt_on_grid[, .(multiple, at_most, again = {
  curve <- dt_enlarged[split_id == .BY$split_id & bar == .BY$bar][order(multiple)]
  lines_needed(curve$multiple, curve$limit, effect_sd)
}, bound_again = {
  curve <- dt_enlarged[split_id == .BY$split_id & bar == .BY$bar][order(multiple)]
  lines_needed_bound(curve$multiple, curve$limit, effect_sd)
}), by = .(split_id, bar)]
dt_today <- merge(dt_on_grid, dt_enlarged[multiple == 1, .(split_id, bar, limit_today = limit)], by = c("split_id", "bar"))
dt_host_sd <- dt_power$host_sd
record("lines needed recomputed from the enlarged limits equals the stored table, under both ways of sizing a known effect", nrow(dt_needed_again) == nrow(dt_on_grid) && dt_needed_again[, all(multiple == again | abs(multiple / again - 1) < 1e-10) && identical(at_most, bound_again)] &&
  dt_today[, all((multiple == 1) == (!is.na(limit_today) & limit_today <= effect_sd))] && dt_needed[below_grid == TRUE | above_grid == TRUE, all(is.na(multiple))] &&
  setequal(unique(dt_needed$convention), c("own spread", "typical gene", "size")) && dt_needed[convention == "typical gene", all(abs(effect_sd - abs(effect_chronos) / dt_host_sd) < 1e-10)],
  sprintf("%d rows; %d found today, %d beyond the largest multiple, %d upper bounds", nrow(dt_needed_again), sum(dt_needed_again$multiple == 1), sum(is.infinite(dt_needed_again$multiple)), sum(dt_needed_again$at_most)))
dt_monotone <- dt_on_grid[order(split_id, bar, -effect_sd), .(ok = !is.unsorted(multiple)), by = .(split_id, bar)]
record("a larger reference effect never needs more lines than a smaller one in the same split", all(dt_monotone$ok), sprintf("%d split-and-bar curves", nrow(dt_monotone)))
dt_cut <- ext_table("driver_type_cut_check")
# The first form of this check set halving a real split against doubling by copies. The review replaced it: a split cut
# to 1/k and copied k times is set against the real split. Criterion written before the new run was read: every
# comparison that the simulated shifts allow lies within the tolerance fixed in the parameters.
dt_compared <- if ("copied_against_real_log2" %in% names(dt_cut)) dt_cut[is.finite(copied_against_real_log2) & at_grid_floor == FALSE] else dt_cut[0]
record("a split cut down and copied back to full size finds about what the real split finds", nrow(dt_compared) >= 10 && dt_compared[bar == "candidate", .N >= 5] && dt_compared[bar == "family", .N >= 3] &&
  all(abs(dt_compared$copied_against_real_log2) <= dp$power$cut_tolerance_log2),
  if (nrow(dt_compared) > 0) paste(dt_compared[, sprintf("%s: %d comparisons, log2 of copied over real %.2f (%.2f to %.2f)", bar[1], .N, mean(copied_against_real_log2), min(copied_against_real_log2), max(copied_against_real_log2)), by = bar]$V1, collapse = "; ") else "no comparison")

# The statements the driver-by-type section makes in words.
dt_sd <- dt_limits[unit == "sd" & target == 0.8]
dt_layer <- dt_sd[bar == "family"]
dt_layer[, smaller := pmin(neg_alt, pos_alt)]
dt_nine <- dt_fitted[gene %in% c("CDKN2A", "CDKN2B", "MTAP")]
dt_altered_of <- function(id) dt_lines[split_id == id & ctx == 1, ModelID]
dt_nine_pairs <- dt_nine[, CJ(a = split_id, b = split_id)[a < b], by = .(level, type)]
dt_nine_pairs[, shared := mapply(function(a, b) length(intersect(dt_altered_of(a), dt_altered_of(b))) / length(union(dt_altered_of(a), dt_altered_of(b))), a, b)]
dt_other_small <- dt_screen[inside_p < 1e-5 & !symbol %in% hits$symbol]
dt_paired_difference <- dt_paired[is.finite(interaction_p)]
dt_raw_hits <- dt_screen[dt_adjust(inside_p_raw) < sp$fdr]
dt_distinct_small <- unique(dt_screen[inside_p < 1e-5, .(symbol, altered_set, p = signif(inside_p, 10))])
dt_statements <- c(
  "every split finds a shift for a gene named in advance" = dt_sd[bar == "candidate", .N == nrow(dt_fitted) && all(reachable)],
  "no split with five lines in the smaller altered cell reaches 80% at the layer's bar, and every split with ten or more does" = dt_layer[smaller == sp$min_cell, .N > 0 && !any(reachable)] && dt_layer[smaller >= 10, .N > 0 && all(reachable)],
  "a split that does not reach 80% still finds the largest shift some of the time" = dt_layer[reachable == FALSE, .N > 0 && all(ceiling > 0 & ceiling < 0.8)],
  "KRAS-altered lung is among the splits that can test a difference" = dt_fitted[split_id == "lineage__Lung__KRAS", .N == 1 && label == "specificity testable"],
  "CDKN2A, CDKN2B and MTAP often mark the same lines" = nrow(dt_nine_pairs) >= 5 && median(dt_nine_pairs$shared) >= 0.75,
  "the tests below 1e-5 outside the paired genes are two, point to less dependence and do not pass" = nrow(dt_other_small) == 2 && all(dt_other_small$inside > 0) && all(dt_other_small$inside_fdr > sp$fdr),
  "two splits share their altered lines, so fewer distinct tests than tests lie below 1e-5" = nrow(dt_distinct_small) < dt_screen[inside_p < 1e-5, .N] && dt_tests$distinct_tests < nrow(dt_screen),
  "the ten smallest p-values among altered lines belong to the paired genes" = all(dt_screen[order(inside_p)][1:10, symbol] %in% hits$symbol),
  "by least squares only paired genes pass, in cell lines and organoids, and no difference does" = dt_screen[inside_fdr_ols < sp$fdr, .N > 0 && all(symbol %in% hits$symbol)] && dt_org[inside_fdr_ols < sp$fdr, all(symbol %in% hits$symbol)] &&
    !any(dt_adjust(dt_screen$interaction_p_ols) < sp$fdr, na.rm = TRUE),
  "the raw-scale hits outside the paired genes sit in splits with few lines" = dt_raw_hits[!symbol %in% hits$symbol, .N > 0 && max(pmin(n_neg_in, n_pos_in)) <= 7],
  "no paired gene's effect differs by driver status in any split" = dt_paired_difference[, .N == 4 * length(dt_four) && min(p.adjust(interaction_p, method = "BH")) > sp$fdr],
  "no organoid or combined test passes the rule" = dt_org[, !any(specific) && !any(inside_only)] && !any(dt_comb$combined_hit),
  "at the layer's bar the copies' limit is on average below the real split's" = dt_compared[bar == "family", .N >= 3 && mean(copied_against_real_log2) < 0],
  "the same Chronos shift is more SDs of a typical gene than of a paired gene's own spread" = dt_needed[class == "Y paralog", .(typical = effect_sd[convention == "typical gene"][1], own = effect_sd[convention == "own spread"][1]), by = reference][, all(typical > own)]
)
record("statements of the driver-by-type section hold in the tables", all(dt_statements), paste(names(dt_statements)[!dt_statements], collapse = "; "))
dt_section_text <- paste(rmd_section(readLines(file.path(script_dir, "22_report.Rmd")), "Driver status inside cancer types"), collapse = " ")
record("the driver-by-type section says \"differs by\" and never \"only in\"", grepl("differ", dt_section_text, fixed = TRUE) && !grepl("only in ", dt_section_text, fixed = TRUE) && !grepl("specific to", dt_section_text, fixed = TRUE), "")
# Wording after the fresh review of these stages: what was withdrawn must be gone, what was added must be there.
dt_withdrawn <- c("no effect size is enough", "of any size", "cannot be detected in any split", "No split would find", "with the same code", "could not have shown")
dt_required <- c("with one exception", "describes the rank-based test", "depends on the gene", "somewhat too low", "These p-values describe; the rule fixed in advance makes the calls", "not calibrated in cells this small", "counted once")
dt_in_section <- function(phrases) vapply(phrases, grepl, logical(1), x = dt_section_text, fixed = TRUE)
dt_marked <- lengths(regmatches(dt_section_text, gregexpr("added after the review", dt_section_text, fixed = TRUE)))
record("driver-by-type wording follows the fresh review of these stages", !any(dt_in_section(dt_withdrawn)) && all(dt_in_section(dt_required)) && dt_marked >= 4 &&
  all(has(c("would have seen only shifts of several SD", "Four parts of the section were added after an independent review"))) && !any(has(c("could not have shown one of any size", "and a difference by driver status in none"))),
  paste(c(sprintf("still present: %s", dt_withdrawn[dt_in_section(dt_withdrawn)]), sprintf("missing: %s", dt_required[!dt_in_section(dt_required)]), sprintf("%d additions marked", dt_marked)), collapse = "; "))

# The extension tables in the final folder are the ones in the extension folder.
ext_final <- list.files(file.path(derived_dir, "final"), pattern = paste0("^(", paste(ext_prefixes, collapse = "|"), ")"))
record("extension tables are in the final folder, identical to their source", length(ext_final) >= 90 && all(file.exists(file.path(ext_dir, ext_final))) &&
  all(tools::md5sum(file.path(derived_dir, "final", ext_final)) == tools::md5sum(file.path(ext_dir, ext_final))), sprintf("%d files", length(ext_final)))

# --- deliverables ------------------------------------------------------------------------
final_dir <- file.path(derived_dir, "final")
final_files <- list.files(final_dir, full.names = TRUE)
record("final tables present and non-empty", length(final_files) >= 100 && all(file.size(final_files) > 0), sprintf("%d files", length(final_files)))
figures <- list.files(file.path(results_dir, "figures"), pattern = "^fig([1-9]|1[0-4])_.*\\.(pdf|svg|png)$", full.names = TRUE)
record("fourteen figures in three formats", length(figures) == 42 && all(file.size(figures) > 5000), sprintf("%d files", length(figures)))
sources <- list.files(file.path(results_dir, "figure_source_data"), pattern = "\\.csv$")
record("source data for every figure panel", length(sources) >= 20, sprintf("%d files", length(sources)))
report <- file.path(results_dir, c("report.html", "report.pdf"))
record("report rendered", all(file.exists(report)) && all(file.size(report) > 50000), paste(basename(report), round(file.size(report) / 1e6, 2), "MB", collapse = "; "))
ledger <- fread(file.path(results_dir, "submission_ledger.tsv"), header = FALSE, sep = "\t")
snapshots <- unique(ledger$V4)
record("every submission has its code snapshot and hash", all(dir.exists(file.path(root_dir, snapshots))) && all(file.exists(file.path(root_dir, paste0(snapshots, ".code.sha256")))), sprintf("%d submissions", nrow(ledger)))

passed <- all(vapply(checks, function(x) x$pass, logical(1)))
jsonlite::write_json(list(passed = passed, checked = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), checks = checks), file.path(results_dir, "validation_summary.json"), auto_unbox = TRUE, pretty = TRUE)
if (!passed) stop("validation failed: ", paste(names(checks)[!vapply(checks, function(x) x$pass, logical(1))], collapse = "; "))
log_msg("validation passed: %d checks", length(checks))
