# Stage 06: outcome matrices and their annotations.
# Every matrix is models (rows) by features (columns), saved as a plain numeric
# .rds so both R environments can read it. Filters applied here never look at
# Y status. Outputs go to derived/loy_reanalysis_v1/matrices/ plus
# gene_annotation.csv, drug_feature_annotation.csv and cohort_sizes.csv.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))

dir.create(matrix_dir, showWarnings = FALSE, recursive = TRUE)
cohort <- read_result("cohort")
eligible <- cohort$ModelID
registry <- fread(input_path("model"), select = "ModelID")$ModelID
gene_par <- params$gene_scan
drug_par <- params$drug_scan

save_matrix <- function(mat, name) {
  stopifnot(is.matrix(mat), !anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)))
  saveRDS(mat, file.path(matrix_dir, paste0(name, ".rds")))
  log_msg("matrix %s: %d models x %d features (%d eligible males)", name, nrow(mat), ncol(mat), sum(rownames(mat) %in% eligible))
  invisible(mat)
}

# --- CRISPR ------------------------------------------------------------------
effect <- read_model_matrix(input_path("crispr_effect"))
dependency <- read_model_matrix(input_path("crispr_dependency"))
audit_ids("CRISPR models", rownames(effect), registry)
# The two tables hold the same models and genes in a different order.
stopifnot(setequal(colnames(effect), colnames(dependency)), setequal(rownames(effect), rownames(dependency)))
dependency <- dependency[rownames(effect), colnames(effect)]
save_matrix(effect, "crispr_effect")
save_matrix(dependency, "crispr_dependency")

gene <- fread(input_path("gene"), select = c("symbol", "entrez_id", "location", "essentiality", "selectivity"))
gene[, chromosome := fifelse(grepl("^([0-9]{1,2}|X|Y)[pq]", location), sub("^([0-9]{1,2}|X|Y)[pq].*$", "\\1", location), NA_character_)]
gene[, par := grepl("^X", location) & grepl(" and Y", location)]
gene[, entrez := as.character(as.integer(entrez_id))]
annotate_genes <- function(features) {
  ann <- data.table(feature = features, symbol = gene_symbol(features), entrez = gene_entrez(features))
  by_entrez <- match(ann$entrez, gene$entrez)
  by_symbol <- match(ann$symbol, gene$symbol)
  idx <- fifelse(!is.na(ann$entrez) & !is.na(by_entrez), by_entrez, by_symbol)
  ann[, c("location", "chromosome", "par", "common_essential", "strongly_selective") := gene[idx, .(
    location, chromosome, par, essentiality == "common essential", selectivity == "strongly selective"
  )]]
  pairs <- params$xy_paralogs
  ann[, y_partner := pairs$y_gene[match(symbol, pairs$x_partner)]]
  ann[, is_y_gene := chromosome %in% "Y" & !(par %in% TRUE)]
  ann[]
}
gene_ann <- annotate_genes(colnames(effect))
in_cohort <- intersect(eligible, rownames(effect))
gene_ann[, coverage := colMeans(!is.na(effect[in_cohort, , drop = FALSE]))]
gene_ann[, n_dependent := colSums(dependency[in_cohort, , drop = FALSE] >= gene_par$dependency_probability_cut, na.rm = TRUE)]
gene_ann[, n_not_dependent := colSums(dependency[in_cohort, , drop = FALSE] < gene_par$dependency_probability_cut, na.rm = TRUE)]
# Scan universe: everything except chrY genes, with enough coverage. Selective
# dependencies (dependent in some lines, not in others) form the context universe.
gene_ann[, in_scan := !is_y_gene & coverage >= gene_par$min_coverage]
gene_ann[, selective := in_scan & n_dependent >= gene_par$selective_min_lines & n_not_dependent >= gene_par$selective_min_lines]
write_result(gene_ann, "gene_annotation")
log_msg("CRISPR genes: %d total, %d chrY removed, %d in scan, %d selective, %d common essential, %d PAR, %d X partners",
  nrow(gene_ann), sum(gene_ann$is_y_gene), sum(gene_ann$in_scan), sum(gene_ann$selective),
  sum(gene_ann$common_essential, na.rm = TRUE), sum(gene_ann$par, na.rm = TRUE), sum(!is.na(gene_ann$y_partner)))

# --- RNAi --------------------------------------------------------------------
rnai <- read_model_matrix(input_path("rnai"))
audit_ids("RNAi models", rownames(rnai), registry)
save_matrix(rnai, "rnai")

# --- expression --------------------------------------------------------------
expr <- read_default_profiles(input_path("expr_tpm_pc"))
save_matrix(expr, "expr_tpm_pc")
# chrY genes and XIST from the all-genes table, which also holds the non-coding ones.
y_symbols <- gene[chromosome %in% "Y" & par == FALSE, symbol]
expr_y <- read_default_profiles(input_path("expr_tpm_all"), features = c(y_symbols, "XIST"))
expr_y <- expr_y[, !duplicated(colnames(expr_y)), drop = FALSE]
save_matrix(expr_y, "expr_y_genes")
rm(expr, expr_y)

# --- compound key ------------------------------------------------------------
portal <- fread(input_path("portal_compounds"))
sample_map <- portal[, .(sample_id = unlist(strsplit(SampleIDs, ";", fixed = TRUE))), by = .(CompoundID, CompoundName, TargetOrMechanism, GeneSymbolOfTargets)]
sample_map <- sample_map[!duplicated(sample_id)]
sample_map[, stem := fifelse(grepl("^BRD:", sample_id), substr(sub("^BRD:", "", sample_id), 1, 13), NA_character_)]
annotate_compounds <- function(dataset, feature, sample_id, name, mechanism = NA_character_, target = NA_character_, screen = NA_character_) {
  out <- data.table(dataset = dataset, feature = feature, sample_id = sample_id, source_name = name, source_mechanism = mechanism, source_target = target, screen = screen)
  idx <- match(out$sample_id, sample_map$sample_id)
  # A PRISM batch missing from the portal list falls back to the compound stem.
  stem_idx <- match(substr(sub("^BRD:", "", out$sample_id), 1, 13), sample_map$stem)
  idx <- fifelse(is.na(idx) & grepl("^BRD:", out$sample_id), stem_idx, idx)
  out[, compound_id := sample_map$CompoundID[idx]]
  out[, compound_name := fifelse(is.na(idx), toupper(source_name), sample_map$CompoundName[idx])]
  portal_mechanism <- sample_map$TargetOrMechanism[idx]
  out[, mechanism := toupper(fifelse(!is.na(portal_mechanism) & portal_mechanism != "", portal_mechanism, source_mechanism))]
  portal_target <- sample_map$GeneSymbolOfTargets[idx]
  out[, targets := fifelse(!is.na(portal_target) & portal_target != "", portal_target, source_target)]
  out[]
}
feature_ann <- list()

# --- PRISM 24Q2 extended primary --------------------------------------------
prism <- fread(input_path("prism_primary"), showProgress = FALSE)
prism_ids <- prism[[1]]
assert_unique(prism_ids, "PRISM primary treatments")
prism_mat <- t(as.matrix(prism[, -1]))
colnames(prism_mat) <- prism_ids
audit_ids("PRISM primary models", rownames(prism_mat), registry)
compound_list <- fread(input_path("prism_compounds"))
# One ID (AZ-628) is listed under two screens and one matrix row is not listed at all.
compound_list[, n_listed := .N, by = IDs]
compound_list <- compound_list[!duplicated(IDs)]
listed <- compound_list[match(prism_ids, IDs)]
log_msg("PRISM primary: %d treatments, %d without annotation, %d listed under more than one screen", length(prism_ids), sum(is.na(listed$screen)), sum(listed$n_listed > 1, na.rm = TRUE))
prism_cohort <- intersect(eligible, rownames(prism_mat))
n_active <- colSums(prism_mat[prism_cohort, , drop = FALSE] <= drug_par$prism_active_lfc, na.rm = TRUE)
for (screen_name in sort(unique(na.omit(listed$screen)))) {
  cols <- which(listed$screen %in% screen_name)
  name <- paste0("prism_", gsub("[^A-Za-z0-9]+", "", screen_name))
  mat <- prism_mat[, cols, drop = FALSE]
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  save_matrix(mat, name)
  ann <- annotate_compounds(name, prism_ids[cols], prism_ids[cols], listed$Drug.Name[cols], listed$MOA[cols], listed$repurposing_target[cols], screen_name)
  ann[, active := n_active[cols] >= drug_par$prism_active_min_lines]
  ann[, screen_ambiguous := listed$n_listed[cols] > 1]
  feature_ann[[name]] <- ann
}
rm(prism, prism_mat)

# --- PRISM 20Q2 secondary ----------------------------------------------------
secondary <- fread(input_path("prism_secondary"), select = c(
  "broad_id", "depmap_id", "screen_id", "auc", "passed_str_profiling", "name", "moa", "target"
), showProgress = FALSE)
n_all <- nrow(secondary)
secondary <- secondary[depmap_id != "" & !is.na(depmap_id) & passed_str_profiling == TRUE & is.finite(auc)]
log_msg("PRISM secondary: %d of %d curves kept (model ID present, STR profiling passed)", nrow(secondary), n_all)
audit_ids("PRISM secondary models", secondary$depmap_id, registry)
secondary[, stem := substr(broad_id, 1, 13)]
# One curve per compound and line: the newest screen available, as the provider recommends.
secondary[, priority := match(screen_id, c("MTS010", "MTS006", "MTS005", "HTS002"))]
stopifnot(!anyNA(secondary$priority))
secondary <- secondary[, .SD[priority == min(priority)], by = .(stem, depmap_id)]
secondary_wide <- secondary[, .(auc = mean(auc), screen_id = screen_id[1]), by = .(stem, depmap_id)]
to_matrix <- function(dt, row, col, value) {
  wide <- dcast(dt, as.formula(paste(row, "~", col)), value.var = value)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide[[1]]
  mat
}
save_matrix(to_matrix(secondary_wide, "depmap_id", "stem", "auc"), "prism_secondary")
secondary_screen <- dcast(secondary_wide, depmap_id ~ stem, value.var = "screen_id")
saveRDS(as.data.frame(secondary_screen), file.path(matrix_dir, "prism_secondary_screen.rds"))
secondary_names <- secondary[, .(broad_id = broad_id[1], name = name[1], moa = moa[1], target = target[1]), by = stem]
ann <- annotate_compounds("prism_secondary", secondary_names$stem, paste0("BRD:", secondary_names$broad_id), secondary_names$name, secondary_names$moa, secondary_names$target, "secondary")
ann[, active := TRUE]
ann[, screen_ambiguous := FALSE]
feature_ann[["prism_secondary"]] <- ann
rm(secondary)

# --- GDSC --------------------------------------------------------------------
gdsc <- fread(input_path("gdsc"), select = c("DATASET", "DRUG_ID", "ARXSPAN_ID", "AUC_PUBLISHED", "DRUG_NAME"), showProgress = FALSE)
n_all <- nrow(gdsc)
gdsc <- gdsc[ARXSPAN_ID != "" & !is.na(ARXSPAN_ID) & is.finite(AUC_PUBLISHED)]
log_msg("GDSC: %d of %d rows kept (model ID present)", nrow(gdsc), n_all)
audit_ids("GDSC models", gdsc$ARXSPAN_ID, registry)
gdsc[, feature := paste0(DATASET, ":", DRUG_ID)]
stopifnot(!anyDuplicated(gdsc[, .(feature, ARXSPAN_ID)]))
for (dataset in c("GDSC1", "GDSC2")) {
  name <- tolower(dataset)
  part <- gdsc[DATASET == dataset]
  save_matrix(to_matrix(part, "ARXSPAN_ID", "feature", "AUC_PUBLISHED"), name)
  names_tab <- unique(part[, .(feature, DRUG_NAME)])
  ann <- annotate_compounds(name, names_tab$feature, names_tab$feature, names_tab$DRUG_NAME, screen = dataset)
  ann[, active := TRUE]
  ann[, screen_ambiguous := FALSE]
  feature_ann[[name]] <- ann
}
rm(gdsc)

# --- CTRP v2 -----------------------------------------------------------------
ctrp <- fread(input_path("ctrp"), select = c(
  "ModelID", "master_ccl_id", "master_cpd_id", "cpd_name", "area_under_curve", "n_experiments",
  "gene_symbol_of_protein_target", "target_or_activity_of_compound"
), showProgress = FALSE)
# Spot-check the derived table against the raw CTRP curves it was built from.
raw_curves <- fread(input_path("ctrp_raw_curves"), select = c("experiment_id", "master_cpd_id", "area_under_curve"), showProgress = FALSE)
raw_experiments <- unique(fread(input_path("ctrp_raw_experiments"), select = c("experiment_id", "master_ccl_id")))
raw_curves <- merge(raw_curves, raw_experiments, by = "experiment_id", allow.cartesian = TRUE)
raw_auc <- raw_curves[, .(raw_auc = mean(area_under_curve), raw_n = .N), by = .(master_ccl_id, master_cpd_id)]
check <- merge(ctrp, raw_auc, by = c("master_ccl_id", "master_cpd_id"), all.x = TRUE)
ctrp_agreement <- check[, .(
  pairs = .N, in_raw = sum(!is.na(raw_auc)),
  equal = sum(abs(area_under_curve - raw_auc) < 1e-3, na.rm = TRUE),
  max_abs_difference = max(abs(area_under_curve - raw_auc), na.rm = TRUE)
)]
write_result(ctrp_agreement, "ctrp_derived_vs_raw")
print(ctrp_agreement)
stopifnot(ctrp_agreement$in_raw == ctrp_agreement$pairs, ctrp_agreement$equal / ctrp_agreement$pairs > 0.99)
rm(raw_curves, check)
cell_map <- fread(input_path("ctrp_cell_map"), select = c("master_ccl_id", "ModelID", "match_tier", "tissue_check"))
trusted <- cell_map[match_tier %in% c("1_CCLEName", "2_StrippedName") & tissue_check == "agree", master_ccl_id]
n_all <- nrow(ctrp)
ctrp <- ctrp[ModelID != "" & !is.na(ModelID) & master_ccl_id %in% trusted & is.finite(area_under_curve)]
log_msg("CTRP: %d of %d rows kept (mapped by name with agreeing tissue)", nrow(ctrp), n_all)
audit_ids("CTRP models", ctrp$ModelID, registry)
ctrp[, feature := paste0("CTRP:", master_cpd_id)]
ctrp_wide <- ctrp[, .(area_under_curve = mean(area_under_curve)), by = .(ModelID, feature)]
save_matrix(to_matrix(ctrp_wide, "ModelID", "feature", "area_under_curve"), "ctrp")
names_tab <- unique(ctrp[, .(feature, cpd_name, target_or_activity_of_compound, gene_symbol_of_protein_target)])
names_tab <- names_tab[!duplicated(feature)]
ann <- annotate_compounds("ctrp", names_tab$feature, names_tab$feature, names_tab$cpd_name, names_tab$target_or_activity_of_compound, names_tab$gene_symbol_of_protein_target, "CTRPv2")
ann[, active := TRUE]
ann[, screen_ambiguous := FALSE]
feature_ann[["ctrp"]] <- ann

feature_ann <- rbindlist(feature_ann, use.names = TRUE)
write_result(feature_ann, "drug_feature_annotation")
print(feature_ann[, .(features = .N, mapped = sum(!is.na(compound_id)), with_mechanism = sum(!is.na(mechanism) & mechanism != ""), active = sum(active)), by = dataset])
shared <- feature_ann[!is.na(compound_id), .(datasets = uniqueN(sub("_.*$", "", dataset))), by = compound_id]
log_msg("compounds mapped to the portal key: %d; in two or more dataset families: %d", nrow(shared), sum(shared$datasets >= 2))

# --- cohort sizes per assay ---------------------------------------------------
datasets <- c("crispr_effect", "rnai", "expr_tpm_pc", unique(feature_ann$dataset))
sizes <- rbindlist(lapply(datasets, function(name) {
  mat <- readRDS(file.path(matrix_dir, paste0(name, ".rds")))
  lines <- cohort[ModelID %in% rownames(mat)]
  lines[, .(
    dataset = name, eligible_males = .N,
    y_negative = sum(y_negative %in% 1), y_positive = sum(y_negative %in% 0), set_aside = sum(is.na(y_negative)),
    donors_negative = uniqueN(PatientID[y_negative %in% 1]), donors_positive = uniqueN(PatientID[y_negative %in% 0]),
    concordant_negative = sum(y_call_concordant %in% "Y_negative"), concordant_positive = sum(y_call_concordant %in% "Y_positive")
  )]
}))
write_result(sizes, "cohort_sizes")
print(sizes)
