# Stage 20 (plan Stage 8): comparison with the earlier analyses in this collection.
# This is the only stage that reads their outputs; nothing here feeds back into a result.
# 1. The new Y-status calls against each earlier call table.
# 2. Pan-cancer gene effects against loy_dependencies_v1.
# 3. Every candidate the two later analyses prioritised, looked up in the new results.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
assert_freeze()

calls <- read_result("y_status_calls")[Sex == "Male"]
calls[, new_call := fifelse(is.na(y_call), "set aside", y_call)]
prior_path <- function(...) file.path(root_dir, ...)

# --- 1. calls -------------------------------------------------------------------
reviewed <- fread(prior_path("derived", "sex_chromosome_review_v1", "sex_chromosome_calls_reviewed.csv"), select = c("ModelID", "review_LOY_binary", "review_LOY_recommended_call"))
dosage <- fread(prior_path("derived", "y_dosage_vulnerabilities_v1", "final", "y_state_model_evidence.csv"), select = c("ModelID", "y_state"))
rna_prior <- fread(prior_path("derived", "rna_y_vulnerabilities_v1", "final", "rna_y_calls.csv"), select = c("ModelID", "y_state"))
legacy <- as.data.table(readRDS(prior_path("depmap_loy_classifications.RDS")))[, .(ModelID, legacy = as.character(LOY_binary))]
sources <- list(
  "sex_chromosome_review_v1 (review_LOY_binary)" = reviewed[, .(ModelID, old = fifelse(is.na(review_LOY_binary), "unresolved", fifelse(review_LOY_binary == 1, "LOY", "not LOY")), negative = "LOY", positive = "not LOY")],
  "y_dosage_vulnerabilities_v1 (y_state)" = dosage[, .(ModelID, old = y_state, negative = "LOY", positive = "RET or GAIN")][, old := fifelse(old %in% c("RET", "GAIN"), "RET or GAIN", old)],
  "rna_y_vulnerabilities_v1 (y_state)" = rna_prior[, .(ModelID, old = y_state, negative = "DEPLETED", positive = "EXPRESSED")],
  "depmap_loy_classifications.RDS" = legacy[, .(ModelID, old = legacy, negative = "LOY", positive = "XY")]
)
cross <- rbindlist(lapply(names(sources), function(name) {
  joined <- merge(calls[, .(ModelID, new_call)], sources[[name]], by = "ModelID", all.x = TRUE)
  joined[is.na(old), old := "not in table"]
  joined[, .N, by = .(old, new_call)][, source := name]
}))
write_result(cross[order(source, old, new_call)], "comparison_call_crosstab")
agreement <- rbindlist(lapply(names(sources), function(name) {
  src <- sources[[name]]
  joined <- merge(calls[, .(ModelID, new_call)], src, by = "ModelID")
  old_definite <- joined$old == joined$negative | joined$old == joined$positive
  new_definite <- joined$new_call != "set aside"
  both <- joined[old_definite & new_definite]
  data.table(
    source = name, males_in_table = nrow(joined), old_definite = sum(old_definite), new_definite_in_table = sum(new_definite),
    both_definite = nrow(both), agree = sum((both$old == both$negative) == (both$new_call == "Y_negative")),
    newly_classified = sum(new_definite & !old_definite), old_only = sum(old_definite & !new_definite)
  )
}))
agreement[, agreement := agree / both_definite]
write_result(agreement, "comparison_call_agreement")
print(agreement)

# --- 2. pan-cancer gene effects --------------------------------------------------
new_genes <- read_result("gene_pan_results")
old_genes <- fread(prior_path("derived", "loy_dependencies_v1", "loy_dependency_results.csv"), select = c("gene", "symbol", "tested", "primary_delta", "primary_fdr", "primary_n_LOY", "primary_n_Y_present"))
genes <- merge(new_genes[, .(feature, symbol, estimate, fdr, label, n_negative, n_positive)], old_genes[tested == TRUE], by.x = "feature", by.y = "gene", suffixes = c("", "_old"))
gene_summary <- data.table(
  genes_in_both = nrow(genes), effect_correlation = cor(genes$estimate, genes$primary_delta, use = "complete.obs"),
  old_lines = paste(genes$primary_n_LOY[1], genes$primary_n_Y_present[1], sep = " vs "), new_lines = paste(max(genes$n_negative), max(genes$n_positive), sep = " vs "),
  old_fdr_hits = sum(genes$primary_fdr < 0.05, na.rm = TRUE), new_fdr_hits = sum(genes$fdr < 0.05),
  hits_in_both = sum(genes$primary_fdr < 0.05 & genes$fdr < 0.05, na.rm = TRUE),
  old_hits = paste(genes[primary_fdr < 0.05, symbol], collapse = ";"), new_hits = paste(genes[fdr < 0.05, symbol], collapse = ";")
)
write_result(gene_summary, "comparison_gene_scan")
print(gene_summary)
partners <- genes[symbol %in% params$xy_paralogs$x_partner, .(symbol, old_estimate = primary_delta, old_fdr = primary_fdr, new_estimate = estimate, new_fdr = fdr, new_label = label)]
write_result(partners[order(new_fdr)], "comparison_xy_partners")
print(partners[order(new_fdr)])

# --- 3. earlier candidates in the new results -------------------------------------
candidates <- rbind(
  fread(prior_path("derived", "y_dosage_vulnerabilities_v1", "final", "prioritized_vulnerabilities.csv"))[, analysis := "y_dosage_vulnerabilities_v1"],
  fread(prior_path("derived", "rna_y_vulnerabilities_v1", "final", "candidate_associations.csv"))[, analysis := "rna_y_vulnerabilities_v1"],
  fill = TRUE
)[, .(analysis, outcome, dataset, contrast, test_type, context, old_effect = effect, old_p = p_value, old_fdr = fdr, symbol, compound_name)]
context_map <- c(
  genomic_MSIScore = "genomic__MSI", genomic_CIN = "genomic_continuous__CIN", WGD = "genomic__WGD", driver_TP53 = "driver__TP53",
  deletion_CDKN2A = "deletion__CDKN2A", Lung = "lineage__Lung", "CNS/Brain" = "lineage__CNS_Brain"
)
dataset_map <- c(prism_extended_REP_PRIMARY = "prism_REPPRIMARY", prism24_REP300 = "prism_REP300", prism24_REP1M = "prism_REP1M", prism20_secondary = "prism_secondary", gdsc1 = "gdsc1", gdsc2 = "gdsc2", ctrp_v2 = "ctrp")
gene_context <- readRDS(file.path(derived_dir, "context_gene_screen.rds"))
gene_focused <- read_result("context_gene_focused")
drug_pan <- read_result("drug_pan_results")
drug_context <- readRDS(file.path(derived_dir, "context_drug_screen.rds"))
new_feature <- function(outcome, dataset) {
  if (dataset == "crispr_screen") return(outcome)
  if (grepl("^GDSC", outcome)) return(sub("_", ":", outcome))
  if (dataset == "prism20_secondary") return(substr(outcome, 1, 13))
  paste0("BRD:", sub("::.*$", "", outcome))
}
# Local names are prefixed with "want_" so they cannot be mistaken for a column inside [.data.table.
lookup <- rbindlist(lapply(seq_len(nrow(candidates)), function(i) {
  row <- candidates[i]
  want_feature <- new_feature(row$outcome, row$dataset)
  want_context <- unname(context_map[row$context])
  out <- data.table(new_dataset = NA_character_, new_context = NA_character_, new_quantity = NA_character_, new_estimate = NA_real_, new_p = NA_real_, new_fdr = NA_real_, new_note = "not tested in the new analysis")
  if (row$dataset == "crispr_screen" && row$test_type == "pan") {
    hit <- new_genes[which(new_genes$feature == want_feature)]
    if (nrow(hit) == 1) out <- data.table(new_dataset = "crispr_effect", new_context = "pan-cancer", new_quantity = "main effect", new_estimate = hit$estimate, new_p = hit$p, new_fdr = hit$fdr, new_note = hit$label)
  } else if (row$dataset == "crispr_screen") {
    pool <- rbind(gene_context, gene_focused, fill = TRUE)
    hit <- pool[which(pool$feature == want_feature & pool$context_id %in% want_context)]
    if (nrow(hit) > 0) {
      hit <- hit[1]
      inside <- grepl("^within", row$test_type)
      out <- data.table(new_dataset = "crispr_effect", new_context = hit$context_id, new_quantity = if (inside) "effect inside context" else "interaction",
        new_estimate = if (inside) hit$inside else hit$interaction, new_p = if (inside) hit$inside_p else hit$interaction_p,
        new_fdr = if (inside) NA_real_ else hit$interaction_fdr_class, new_note = if (inside) "effect inside context" else if (isTRUE(hit$context_specific)) "context-specific" else "not context-specific")
    }
  } else {
    want_dataset <- unname(dataset_map[row$dataset])
    hit <- drug_context[which(drug_context$dataset %in% want_dataset & drug_context$feature == want_feature & drug_context$context_id %in% want_context)]
    if (nrow(hit) > 0) {
      hit <- hit[1]
      out <- data.table(new_dataset = want_dataset, new_context = hit$context_id, new_quantity = "interaction", new_estimate = hit$interaction_standardized, new_p = hit$interaction_p,
        new_fdr = hit$interaction_fdr_class, new_note = if (isTRUE(hit$context_specific)) "context-specific" else "not context-specific")
    }
  }
  cbind(row, out)
}))
write_result(lookup, "comparison_prior_candidates")
print(lookup[, .(analysis = sub("_vulnerabilities_v1", "", analysis), name = fifelse(is.na(symbol) | symbol == "", compound_name, symbol), test_type = substr(test_type, 1, 22), context = substr(context, 1, 16), old_fdr = signif(old_fdr, 2), new_quantity, new_p = signif(new_p, 2), new_fdr = signif(new_fdr, 2), new_note)])
log_msg("comparison done: %d earlier candidates looked up, %d with a counterpart", nrow(lookup), sum(lookup$new_note != "not tested in the new analysis"))
