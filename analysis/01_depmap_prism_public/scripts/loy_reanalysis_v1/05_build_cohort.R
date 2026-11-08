# Stage 05: eligible cohort, covariates and the context panel.
# Eligible = annotated male, cancer, not engineered. All models of a donor are
# kept; PatientID is the blocking factor. No outcome table is read here.
# Outputs: cohort.csv, context_matrix.csv, context_definitions.csv.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))

ctx_par <- params$context
calls <- read_result("y_status_calls")
model <- fread(input_path("model"), select = c("ModelID", "GrowthPattern", "AgeCategory", "PrimaryOrMetastasis"))

# --- eligibility -------------------------------------------------------------
call_cols <- grep("^y_call", names(calls), value = TRUE)
cohort <- calls[, c(
  "ModelID", "PatientID", "StrippedCellLineName", "Sex", "OncotreeLineage", "OncotreePrimaryDisease",
  "OncotreeSubtype", "EngineeredModel", "ModelType", call_cols, "y_tier", "rna_call", "dna_call",
  "dna_source", "dna_partial_type", "rna_y_mean", "rna_y_score", "rna_n_genes_on", "dna_y_ratio", "y_copies",
  "rna_library_size", "xist_high", "x_two_copy", "possible_female", "call_stable",
  paste0("rna_", params$y_genes_core), paste0("rna_on_", params$y_genes_core)
), with = FALSE]
exclusions <- cohort[, .(
  all_models = .N,
  not_male = sum(Sex != "Male"),
  male_non_cancerous = sum(Sex == "Male" & OncotreePrimaryDisease == "Non-Cancerous"),
  male_engineered = sum(Sex == "Male" & OncotreePrimaryDisease != "Non-Cancerous" & EngineeredModel != "")
)]
cohort <- cohort[Sex == "Male" & OncotreePrimaryDisease != "Non-Cancerous" & (is.na(EngineeredModel) | EngineeredModel == "")]
cohort <- merge(cohort, model, by = "ModelID", all.x = TRUE, sort = FALSE)
cohort[, in_primary := !is.na(y_call)]
cohort[, y_negative := fifelse(y_call == "Y_negative", 1, fifelse(y_call == "Y_positive", 0, NA_real_))]
exclusions[, eligible_males := nrow(cohort)]
exclusions[, eligible_in_primary := sum(cohort$in_primary)]
write_result(exclusions, "cohort_exclusions")
print(exclusions)
stopifnot(!anyDuplicated(cohort$ModelID), all(cohort$PatientID != ""))

# --- covariates --------------------------------------------------------------
signatures <- fread(input_path("signatures"))[IsDefaultEntryForModel == "Yes", .(ModelID, MSIScore, LoHFraction, WGD, CIN, Ploidy, Aneuploidy)]
confounders <- as.data.table(read_model_matrix(input_path("crispr_confounders")), keep.rownames = "ModelID")
cohort <- merge(cohort, signatures, by = "ModelID", all.x = TRUE, sort = FALSE)
cohort <- merge(cohort, confounders, by = "ModelID", all.x = TRUE, sort = FALSE)

# Driver alterations: any variant carrying one of DepMap's driver flags, in the
# model's default profile. A model without mutation profiling is NA, never wild-type.
as_flag <- function(x) toupper(as.character(x)) %in% c("TRUE", "Y", "YES", "1")
mut <- fread(input_path("mutations"), select = c(
  "ModelID", "IsDefaultEntryForModel", "HugoSymbol", "Hotspot", "OncogeneHighImpact", "TumorSuppressorHighImpact", "HessDriver"
), showProgress = FALSE)
mut <- mut[IsDefaultEntryForModel == "Yes"]
mut[, driver := as_flag(Hotspot) | as_flag(OncogeneHighImpact) | as_flag(TumorSuppressorHighImpact) | as_flag(HessDriver)]
driver_pairs <- unique(mut[driver == TRUE & HugoSymbol != "", .(ModelID, HugoSymbol)])
profiled <- fread(input_path("mutation_hotspot"), select = c("ModelID", "IsDefaultEntryForModel"))[IsDefaultEntryForModel == "Yes", unique(ModelID)]
audit_ids("mutation-profiled models", profiled, calls$ModelID)
rm(mut)
driver_status <- function(symbol) {
  altered <- driver_pairs[HugoSymbol == symbol, ModelID]
  fifelse(cohort$ModelID %in% altered, 1, fifelse(cohort$ModelID %in% profiled, 0, NA_real_))
}
cohort[, tp53_altered := driver_status("TP53")]
log_msg("eligible males: %d (%d in the primary call); mutation-profiled %d; TP53 altered %d", nrow(cohort), sum(cohort$in_primary), sum(!is.na(cohort$tp53_altered)), sum(cohort$tp53_altered, na.rm = TRUE))

# --- context panel -----------------------------------------------------------
primary <- cohort$in_primary
context_list <- list()
definitions <- list()
make_id <- function(class, label) paste0(class, "__", gsub("[^A-Za-z0-9]+", "_", label))
add_context <- function(class, label, values, type, source) {
  id <- make_id(class, label)
  context_list[[id]] <<- as.numeric(values)
  definitions[[id]] <<- data.table(context_id = id, class = class, label = label, type = type, source = source)
}

for (level in sort(unique(cohort$OncotreeLineage))) add_context("lineage", level, cohort$OncotreeLineage == level, "binary", "Model.csv OncotreeLineage")
for (level in sort(unique(cohort$OncotreePrimaryDisease))) add_context("disease", level, cohort$OncotreePrimaryDisease == level, "binary", "Model.csv OncotreePrimaryDisease")
for (level in sort(unique(cohort$OncotreeSubtype[cohort$OncotreeSubtype != ""]))) add_context("subtype", level, cohort$OncotreeSubtype == level, "binary", "Model.csv OncotreeSubtype")

# DepMap's inferred molecular subtypes keep NA for models that could not be assessed.
subtypes <- fread(input_path("molecular_subtypes"))
subtype_cols <- setdiff(names(subtypes), c("ModelID", "MSI"))
subtypes <- subtypes[match(cohort$ModelID, ModelID)]
for (col in subtype_cols) add_context("molecular_subtype", col, as.logical(subtypes[[col]]), "binary", "OmicsInferredMolecularSubtypes.csv")

# Driver genes chosen by prevalence only, among eligible males in the primary call.
driver_counts <- driver_pairs[ModelID %in% cohort[in_primary == TRUE, ModelID], .N, by = HugoSymbol][N >= ctx_par$driver_min_altered][order(-N)]
log_msg("driver genes with at least %d altered eligible males: %d", ctx_par$driver_min_altered, nrow(driver_counts))
for (symbol in driver_counts$HugoSymbol) add_context("driver", symbol, driver_status(symbol), "binary", "OmicsSomaticMutations.csv driver flags")

# Gene-level amplification and deep deletion on the relative copy ratio (WGS, else WES).
cn_genes <- union(ctx_par$amp_genes, ctx_par$del_genes)
wgs <- read_default_profiles(input_path("cn_wgs"), features = cn_genes)
colnames(wgs) <- gene_symbol(colnames(wgs))
wes_header <- read_header(input_path("cn_wes"))
wes_keep <- wes_header[gene_symbol(wes_header) %in% cn_genes]
wes <- fread(input_path("cn_wes"), select = c("ModelConditionID", "IsDefaultEntryForMC", wes_keep), showProgress = FALSE)[IsDefaultEntryForMC == "Yes"]
wes <- merge(wes, fread(input_path("model_condition"), select = c("ModelConditionID", "ModelID")), by = "ModelConditionID")
wes <- wes[, lapply(.SD, median, na.rm = TRUE), by = ModelID, .SDcols = wes_keep]
setnames(wes, wes_keep, gene_symbol(wes_keep))
copy_ratio <- function(symbol) {
  from_wgs <- if (symbol %in% colnames(wgs)) wgs[match(cohort$ModelID, rownames(wgs)), symbol] else rep(NA_real_, nrow(cohort))
  from_wes <- if (symbol %in% names(wes)) wes[[symbol]][match(cohort$ModelID, wes$ModelID)] else rep(NA_real_, nrow(cohort))
  fifelse(!is.na(from_wgs), from_wgs, from_wes)
}
for (symbol in ctx_par$amp_genes) add_context("amplification", symbol, copy_ratio(symbol) >= ctx_par$amp_ratio_min, "binary", sprintf("relative copy ratio >= %s", ctx_par$amp_ratio_min))
for (symbol in ctx_par$del_genes) add_context("deletion", symbol, copy_ratio(symbol) <= ctx_par$del_ratio_max, "binary", sprintf("relative copy ratio <= %s", ctx_par$del_ratio_max))

add_context("genomic", "WGD", cohort$WGD, "binary", "OmicsGlobalSignatures.csv WGD")
add_context("genomic", "MSI", cohort$MSIScore >= ctx_par$msi_score_min, "binary", sprintf("OmicsGlobalSignatures.csv MSIScore >= %s", ctx_par$msi_score_min))
for (col in c("CIN", "Ploidy", "Aneuploidy", "LoHFraction")) add_context("genomic_continuous", col, cohort[[col]], "continuous", paste("OmicsGlobalSignatures.csv", col))

context <- as.data.table(context_list)
definitions <- rbindlist(definitions)

# Keep a context only if the design can support it: every cell of Y status by
# context has the minimum number of lines (per-assay gates are applied again later).
y_neg <- cohort$y_negative
cell_counts <- function(values, type) {
  if (type == "continuous") {
    return(c(sum(primary & y_neg %in% 1 & !is.na(values)), sum(primary & y_neg %in% 0 & !is.na(values)), NA, NA))
  }
  c(sum(y_neg %in% 1 & values %in% 1), sum(y_neg %in% 0 & values %in% 1), sum(y_neg %in% 1 & values %in% 0), sum(y_neg %in% 0 & values %in% 0))
}
cells <- t(vapply(seq_len(nrow(definitions)), function(i) cell_counts(context[[i]], definitions$type[i]), numeric(4)))
definitions[, c("n_neg_in", "n_pos_in", "n_neg_out", "n_pos_out") := as.data.table(cells)]
definitions[, n_assessed := vapply(seq_len(.N), function(i) sum(primary & !is.na(context[[i]])), numeric(1))]
definitions[, supported := fifelse(
  type == "continuous", n_neg_in >= ctx_par$min_group_continuous & n_pos_in >= ctx_par$min_group_continuous,
  pmin(n_neg_in, n_pos_in, n_neg_out, n_pos_out) >= ctx_par$min_cell
)]
# Drop binary contexts that are identical to an earlier one (a disease that is its whole lineage).
signature <- vapply(seq_len(nrow(definitions)), function(i) paste(context[[i]], collapse = ""), character(1))
definitions[, duplicate_of := context_id[match(signature, signature)]]
definitions[duplicate_of == context_id, duplicate_of := NA_character_]
definitions[, kept := supported & is.na(duplicate_of)]
print(definitions[, .(candidates = .N, supported = sum(supported), kept = sum(kept)), by = class])

kept_ids <- definitions[kept == TRUE, context_id]
context_out <- cbind(cohort[, .(ModelID)], context[, kept_ids, with = FALSE])
write_result(cohort, "cohort")
write_result(context_out, "context_matrix")
write_result(definitions, "context_definitions")
log_msg("cohort %d eligible males (%d Y-negative, %d Y-positive, %d set aside); %d contexts kept of %d candidates",
  nrow(cohort), sum(cohort$y_negative %in% 1), sum(cohort$y_negative %in% 0), sum(!cohort$in_primary), length(kept_ids), nrow(definitions))
print(cohort[in_primary == TRUE, .(
  n = .N, wgd = mean(WGD, na.rm = TRUE), tp53 = mean(tp53_altered, na.rm = TRUE), loh = mean(LoHFraction, na.rm = TRUE),
  aneuploidy = mean(Aneuploidy, na.rm = TRUE), library_ky = mean(LibraryKY, na.rm = TRUE), donors = uniqueN(PatientID)
), by = y_call])
