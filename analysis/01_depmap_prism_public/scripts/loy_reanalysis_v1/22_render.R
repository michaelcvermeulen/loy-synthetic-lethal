# Stage 22: export the deliverable tables and render the report.
# The report only assembles pre-computed tables and figures.
# "collect" is run first for the expression forests (17_forest_expression.R collect).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("organoid")
assert_extension_freeze("beataml")
assert_extension_freeze("nextgen")
assert_extension_freeze("driver_type")

# --- deliverable tables ----------------------------------------------------------
final_dir <- file.path(derived_dir, "final")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
tables <- c(
  "y_status_calls", "call_summary", "cohort", "cohort_sizes", "context_definitions", "gene_annotation", "drug_feature_annotation",
  "pipeline_controls", "gene_pan_results", "gene_pan_hits", "gene_pan_continuous", "gene_pan_rnai", "gene_pan_hits_rnai", "gene_pan_pathways", "gene_pan_xy_partners",
  "drug_pan_results", "drug_pan_consensus", "drug_pan_mechanism",
  "context_gene_summary", "context_gene_focused", "context_gene_hits", "context_gene_calibration", "context_gene_lineage_heterogeneity", "context_gene_focused_by_lineage",
  "context_drug_summary", "context_drug_hits", "context_drug_calibration",
  "expression_de_results", "expression_de_hits", "expression_de_pathways", "expression_de_positional", "expression_de_dosage_control", "expression_de_xy_partners",
  "expression_dependency_set", "expression_dependency_summary", "expression_dependency_top", "expression_dependency_y_genes", "expression_dependency_de_overlap", "expression_dependency_hallmark",
  "forest_expression_accuracy", "forest_expression_features", "forest_expression_partial_dependence", "forest_expression_null", "forest_confounders", "forest_confounders_summary",
  "robustness_gene_hits", "robustness_gene_scan_by_definition", "robustness_drug_scan_by_definition", "robustness_context_gene_hits",
  "comparison_call_agreement", "comparison_call_crosstab", "comparison_gene_scan", "comparison_xy_partners", "comparison_prior_candidates",
  "power_census", "power_cancer_types", "power_shift_detection", "power_limits", "power_benchmarks", "power_benchmark_detection", "power_benchmark_summary",
  "power_benchmark_null", "power_benchmark_levels", "power_positive_controls", "power_alone_hits", "power_genotype_context_limits", "power_drug_limits"
)
for (name in tables) {
  for (extension in c(".csv", ".rds")) {
    stopifnot(file.exists(file.path(derived_dir, paste0(name, extension))))
    file.copy(file.path(derived_dir, paste0(name, extension)), file.path(final_dir, paste0(name, extension)), overwrite = TRUE)
  }
}
large <- c("context_gene_screen", "context_drug_screen", "expression_dependency_associations")
for (name in large) {
  file.copy(file.path(derived_dir, paste0(name, ".csv.gz")), file.path(final_dir, paste0(name, ".csv.gz")), overwrite = TRUE)
  file.copy(file.path(derived_dir, paste0(name, ".rds")), file.path(final_dir, paste0(name, ".rds")), overwrite = TRUE)
}
# Tables of stages 25 to 36, copied from the extension folder. Their prefixes keep them apart from the tables above.
extension_tables <- c(
  "organoid_y_calls", "organoid_cohort", "organoid_gene_annotation", "organoid_test_counts", "organoid_input_manifest", "organoid_weights", "organoid_controls",
  "organoid_pan_results", "organoid_continuous", "organoid_xy_partners", "organoid_type_screen", "organoid_genotype_screen", "organoid_hits", "organoid_calibration", "organoid_hit_followup",
  "rnai_context_summary", "rnai_context_focused", "rnai_context_hits", "rnai_context_inside_only", "rnai_context_calibration", "rnai_hit_followup",
  "combined_pan", "combined_seen_before", "combined_within_type", "combined_cross_assay_summary", "combined_cell_line_leads_in_rnai", "combined_xy_partners", "extension_power",
  "beataml_cohort", "beataml_karyotypes", "beataml_test_counts", "beataml_families", "beataml_compound_match", "beataml_input_manifest", "beataml_controls",
  "beataml_drug_results", "beataml_calibration", "beataml_family_results", "beataml_cell_line_comparison", "beataml_cell_line_comparison_summary", "beataml_ras_by_y_loss",
  "beataml_hit_followup", "beataml_leading_inhibitor", "beataml_power",
  # Stages 31 to 33.
  "nextgen_input_manifest", "nextgen_y_calls", "nextgen_donor_groups", "nextgen_cohort", "nextgen_gene_annotation", "nextgen_test_counts", "nextgen_pooled_types", "nextgen_replication_chance",
  "nextgen_weights", "nextgen_controls", "nextgen_calibration", "nextgen_pan_results", "nextgen_paired_genes", "nextgen_xy_partner_scan", "nextgen_continuous", "nextgen_type_screen",
  "nextgen_hits", "nextgen_hit_followup", "nextgen_combined_pan", "nextgen_combined_followup", "nextgen_combined_hits", "nextgen_combined_within_type", "nextgen_pooled_type", "nextgen_pooled_named",
  "nextgen_pooled_sensitivity", "nextgen_shared_donors", "nextgen_xy_partners", "nextgen_power",
  # Stages 34 to 36.
  "driver_type_events", "driver_type_status", "driver_type_organoid_events", "driver_type_organoid_status", "driver_type_census", "driver_type_census_summary",
  "driver_type_splits", "driver_type_split_lines", "driver_type_organoid_splits", "driver_type_organoid_split_lines", "driver_type_pairs", "driver_type_genes",
  "driver_type_controls", "driver_type_paired_genes", "driver_type_calibration", "driver_type_hits", "driver_type_split_results", "driver_type_organoid_screen",
  "driver_type_combined", "driver_type_hit_followup",
  "driver_type_shift_detection", "driver_type_benchmark_detection", "driver_type_enlarged_detection", "driver_type_limits", "driver_type_power_splits",
  "driver_type_enlarged_limits", "driver_type_lines_needed", "driver_type_cut_check", "driver_type_organoid_limits", "driver_type_gate_distance"
)
extension_large <- c("rnai_context_screen", "combined_cross_assay", "driver_type_screen")
stopifnot(!any(c(extension_tables, extension_large) %in% c(tables, large)), all(grepl(paste0("^(", paste(ext_prefixes, collapse = "|"), ")"), c(extension_tables, extension_large))))
for (name in extension_tables) {
  for (extension in c(".csv", ".rds")) {
    stopifnot(file.exists(file.path(ext_dir, paste0(name, extension))))
    file.copy(file.path(ext_dir, paste0(name, extension)), file.path(final_dir, paste0(name, extension)), overwrite = TRUE)
  }
}
for (name in extension_large) {
  for (extension in c(".csv.gz", ".rds")) {
    stopifnot(file.exists(file.path(ext_dir, paste0(name, extension))))
    file.copy(file.path(ext_dir, paste0(name, extension)), file.path(final_dir, paste0(name, extension)), overwrite = TRUE)
  }
}
forest_files <- list.files(file.path(derived_dir, "forest_y_status"), full.names = TRUE)
forest_y <- rbindlist(lapply(grep("forest_y_status_summary_.*\\.rds$", forest_files, value = TRUE), readRDS))
forest_features <- rbindlist(lapply(grep("forest_y_status_(?!summary).*\\.rds$", forest_files, value = TRUE, perl = TRUE), readRDS), fill = TRUE)
write_result(forest_y, "forest_y_status_summary", dir = final_dir)
write_result(forest_features, "forest_y_status_features", dir = final_dir)
log_msg("final tables: %d files in %s", length(list.files(final_dir)), final_dir)

# --- report -------------------------------------------------------------------------
file.copy(file.path(script_dir, "22_report.Rmd"), file.path(results_dir, "report.Rmd"), overwrite = TRUE)
Sys.setenv(LOY_REPORT_DERIVED = derived_dir)
old <- setwd(results_dir)
outputs <- rmarkdown::render("report.Rmd", output_format = "all", quiet = TRUE, envir = new.env(parent = globalenv()))
setwd(old)
stopifnot(all(file.exists(outputs)), all(file.size(outputs) > 50000))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info_report.txt"))
log_msg("report rendered: %s", paste(basename(outputs), collapse = ", "))
