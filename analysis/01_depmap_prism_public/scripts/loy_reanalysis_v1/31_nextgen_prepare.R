# Stage 31: DepMap NextGen models (Neiswender et al. 2026), preparation.
# Screened organoid and neurosphere models: recorded sex, Y status (RNA-led with a
# DNA check), donors, the cohort, the matrices, the gene annotation, the tests to
# come, the cancer types to pool across sources and the chance of replicating the
# four paired genes, then a freeze that stages 32 and 33 assert. No dependency
# value is compared with Y status here. Outputs: extension/nextgen_input_manifest,
# nextgen_y_calls, nextgen_donor_groups, nextgen_cohort, nextgen_gene_annotation,
# nextgen_test_counts, nextgen_pooled_types, nextgen_replication_chance,
# nextgen_effect.rds, nextgen_dependency.rds, nextgen_prepare_summary.json,
# extension_freeze_nextgen.json.
# Pass "--replace" to archive an existing freeze and write a new one.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
suppressPackageStartupMessages(library(readxl))
assert_freeze()
assert_extension_freeze("organoid", dir = ext_dir)
dir.create(ext_out_dir, showWarnings = FALSE, recursive = TRUE)

np <- ext_params$nextgen
sp <- ext_params$shared
dna_par <- params$dna_call
core <- sp$y_genes
study_dir <- file.path(raw_dir, "studies", "depmap_nextgen")
study_path <- function(...) file.path(study_dir, ...)
within_tolerance <- function(found, expected) abs(found - expected) <= sp$count_tolerance
unit_of <- function(design, coef) sqrt(solve(crossprod(design))[coef, coef])

# --- inputs ----------------------------------------------------------------------------
# The six files supplied on 2026-10-05 against the checksums recorded with them; the
# files that were already in the collection against the collection's own list.
inputs <- as.data.table(ext_params$nextgen_inputs)
supplied <- rbindlist(lapply(seq_len(nrow(inputs)), function(i) verify_download(file.path(root_dir, inputs$file[i]), inputs$checksum[i], inputs$algorithm[i])))
supplied[, source := inputs$source]
listed <- tstrsplit(readLines(file.path(root_dir, "docs", "SHA256SUMS")), "  ", fixed = TRUE)
listed <- setNames(listed[[1]], listed[[2]])
earlier <- file.path("raw", "studies", "depmap_nextgen", c(
  "screen_metadata.csv", "crispr_screen_qc.csv", "screen_gene_effect.csv", "screen_gene_dependency.csv",
  "next_gen_damaging.csv", "next_gen_hotspot.csv", "next_gen_omics_signatures.csv"
))
stopifnot(all(earlier %in% names(listed)))
in_collection <- rbindlist(lapply(earlier, function(path) verify_download(file.path(root_dir, path), listed[[path]], "sha256")))
in_collection[, source := "DepMap portal file set 'NextGen Model Manuscript 2026', in the collection since September 2026; SHA-256 from docs/SHA256SUMS"]
write_ext(rbind(supplied, in_collection), "nextgen_input_manifest")

# --- screens: NextGen models that pass QC, one screen each -----------------------------------
screens <- fread(study_path("screen_metadata.csv"))[IsNextGen == TRUE & PassesQC == TRUE]
qc <- fread(study_path("crispr_screen_qc.csv"))
screens <- merge(screens, qc[, .(ScreenID, CasActivity, ScreenDoublingTime, ScreenNNMD, ScreenROCAUC)], by = "ScreenID", all.x = TRUE, sort = FALSE)
assert_unique(screens$ModelID, "NextGen screens per model")
stopifnot(uniqueN(screens$Library) == 1, !anyNA(screens$ScreenNNMD), !anyNA(screens$ScreenType))
log_msg("screens: %d NextGen models pass QC, library %s; formats %s", nrow(screens), unique(screens$Library),
  paste(sprintf("%s %d", names(table(screens$ScreenType)), table(screens$ScreenType)), collapse = ", "))

# --- model annotation: the paper's Supplementary Table 1 ----------------------------------------
models <- as.data.table(suppressMessages(read_excel(study_path("supplementary_tables", "41586_2026_10843_MOESM3_ESM.xlsx"), sheet = "NextGen Model Metadata")))
models <- models[!is.na(`Model ID`), .(
  ModelID = `Model ID`, name = `Cell Line Name`, alias = Alias, source_institution = `Source Institution`, source_lab = `Source Lab`,
  culture_media = `Culture Media`, culture_format = `Culture Format`, plate_coating = `Plate Coating`, age = suppressWarnings(as.numeric(Age)),
  Sex, tumour_status = `Tumor Status`, treatment_status = `Treatment Status`
)]
assert_unique(models$ModelID, "NextGen model annotation")
stopifnot(all(screens$ModelID %in% models$ModelID), !anyNA(models$name), all(models$Sex %in% c("Male", "Female", "Unknown")))
models[, recorded_male := Sex == "Male"]

# --- RNA evidence and call --------------------------------------------------------------------
log_expr <- read_model_matrix(study_path("next_gen_expression.csv"), features = core)
colnames(log_expr) <- gene_symbol(colnames(log_expr))
stopifnot(setequal(colnames(log_expr), core), !anyDuplicated(colnames(log_expr)), !anyNA(log_expr), all(log_expr >= 0))
log_expr <- log_expr[, core]
with_rna <- intersect(models$ModelID, rownames(log_expr))
mean_y <- rowMeans(log_expr[with_rna, ])
recorded_female <- models[Sex == "Female" & ModelID %in% with_rna, ModelID]
rna_males <- models[Sex == "Male" & ModelID %in% with_rna, ModelID]
# A recorded female above the midpoint of the two sexes expresses Y like a male: not a reference female.
midpoint <- (median(mean_y[recorded_female]) + median(mean_y[rna_males])) / 2
discordant_female <- recorded_female[mean_y[recorded_female] > midpoint]
female_reference <- setdiff(recorded_female, discordant_female)
thresholds_cell_line <- jsonlite::fromJSON(file.path(derived_dir, "call_thresholds.json"))
# The cell-line threshold of expressed_mean_min log2(TPM + 1), as a share of the typical male level in TPM.
expressed_fraction <- (2^params$rna_call$expressed_mean_min - 1) / (2^thresholds_cell_line$expressed_male_median - 1)
rna_rule <- function(ceiling_fraction) call_rna_y(
  log_expr[with_rna, ], female_ids = female_reference, male_ids = rna_males, expressed_fraction = expressed_fraction,
  depleted_max_on = np$depleted_max_genes_on, expressed_min_on = np$expressed_min_genes_on, ceiling_fraction = ceiling_fraction
)
rna <- rna_rule(np$rna_ceiling_fraction)
rna_no_ceiling <- rna_rule(NULL)
lenient_call <- fifelse(rna$mean_y <= rna$thresholds$depleted_mean_max, "depleted", fifelse(rna$mean_y >= rna$thresholds$expressed_mean_min, "expressed", "intermediate"))
log_msg("RNA rule on the NextGen scale: %d reference females (%d recorded females express Y like males: %s), %d males; female mean max %.3f, ceiling %.3f (%.0f%% of the male level), male reference %.3f, expressed threshold %.3f (%.1f%%); %d passes",
  length(female_reference), length(discordant_female), paste(discordant_female, collapse = ", "), length(rna_males), rna$thresholds$female_mean_max, rna$thresholds$depleted_mean_max,
  100 * np$rna_ceiling_fraction, rna$thresholds$male_reference, rna$thresholds$expressed_mean_min, 100 * expressed_fraction, rna$thresholds$passes)

# --- DNA evidence and call: the seven core Y genes, copy number relative to ploidy -----------------
gene <- fread(input_path("gene"), select = c("symbol", "entrez_id", "location", "essentiality"))
gene[, chromosome := fifelse(grepl("^([0-9]{1,2}|X|Y)[pq]", location), sub("^([0-9]{1,2}|X|Y)[pq].*$", "\\1", location), NA_character_)]
gene[, par := grepl("^X", location) & grepl(" and Y", location)]
gene[, entrez := as.character(as.integer(entrez_id))]
annotate_genes <- function(features) {
  ann <- data.table(feature = features, symbol = gene_symbol(features), entrez = gene_entrez(features))
  by_entrez <- match(ann$entrez, gene$entrez)
  by_symbol <- match(ann$symbol, gene$symbol)
  idx <- fifelse(!is.na(ann$entrez) & !is.na(by_entrez), by_entrez, by_symbol)
  ann[, c("chromosome", "par", "common_essential") := gene[idx, .(chromosome, par, essentiality == "common essential")]]
  ann[, y_partner := params$xy_paralogs$y_gene[match(symbol, params$xy_paralogs$x_partner)]]
  ann[, is_y_gene := chromosome %in% "Y" & !(par %in% TRUE)]
  ann[]
}
copy_number <- read_model_matrix(study_path("next_gen_copy_number.csv"))
cn_ann <- annotate_genes(colnames(copy_number))
stopifnot(all(core %in% cn_ann$symbol), !anyDuplicated(cn_ann[symbol %in% core, symbol]))
y_copy <- copy_number[, cn_ann$feature[match(core, cn_ann$symbol)]]
colnames(y_copy) <- core
dna <- dna_gene_call(y_copy, loss_max = dna_par$loss_ratio_max, retained_min = dna_par$retained_ratio_min, min_genes = np$dna_core_genes_min)
setnames(dna, c("sample", "n_genes", "n_loss", "n_retained", "ratio", "call"), c("ModelID", "dna_y_genes", "dna_y_genes_lost", "dna_y_genes_retained", "dna_y_ratio", "dna_call"))
# X copy number relative to ploidy: about 0.5 for one X in a diploid genome. The table cannot tell a
# female from a male with a gained X, so it is used only to mark Y-negative males for one follow-up variant.
x_features <- cn_ann[chromosome %in% "X" & !(par %in% TRUE), feature]
dna[, x_ratio := apply(copy_number[, x_features], 1, median, na.rm = TRUE)]
log_msg("DNA: %d models with copy number, %d chrX genes for the X ratio; calls %s", nrow(dna), length(x_features),
  paste(sprintf("%s %d", names(table(dna$dna_call, useNA = "ifany")), table(dna$dna_call, useNA = "ifany")), collapse = ", "))

# --- genomic state and mutations ---------------------------------------------------------------------
signatures <- as.data.table(read_model_matrix(study_path("next_gen_omics_signatures.csv")), keep.rownames = "ModelID")
damaging <- read_model_matrix(study_path("next_gen_damaging.csv"))
hotspot <- read_model_matrix(study_path("next_gen_hotspot.csv"))
stopifnot(!anyNA(damaging), !anyNA(hotspot))
# 1 when the gene carries a call, 0 when the model was profiled without one, missing when not profiled.
flag_of <- function(mat, symbol, ids) {
  column <- match(symbol, gene_symbol(colnames(mat)))
  values <- if (is.na(column)) setNames(rep(0, nrow(mat)), rownames(mat)) else mat[, column]
  unname(as.numeric(values[ids] > 0))
}

# --- integration and alternative definitions ------------------------------------------------------
calls <- merge(models, dna, by = "ModelID", all.x = TRUE, sort = FALSE)
at <- match(calls$ModelID, with_rna)
calls[, has_rna := !is.na(at)]
calls[, rna_y_mean := rna$mean_y[at]]
calls[, rna_n_genes_on := rna$n_on[at]]
calls[, rna_call := rna$call[at]]
calls[, rna_call_lenient := lenient_call[at]]
calls[, rna_call_no_ceiling := rna_no_ceiling$call[at]]
for (g in core) calls[, (paste0("rna_", g)) := log_expr[with_rna, g][at]]
calls[, female_reference := ModelID %in% female_reference]
calls[, sex_discordant := ModelID %in% discordant_female]
primary_call <- integrate_y_calls(calls$rna_call, calls$dna_call, calls$recorded_male)
calls[, y_call := primary_call$call]
calls[, y_tier := primary_call$tier]
as_call <- function(rna_state) fifelse(rna_state == "depleted", "Y_negative", fifelse(rna_state == "expressed", "Y_positive", NA_character_))
calls[, y_call_rna := fifelse(recorded_male, as_call(rna_call), NA_character_)]
calls[, y_call_concordant := fifelse(y_tier %in% "concordant", y_call, NA_character_)]
calls[, y_call_lenient := integrate_y_calls(rna_call_lenient, dna_call, recorded_male)$call]
calls[, y_call_no_ceiling := integrate_y_calls(rna_call_no_ceiling, dna_call, recorded_male)$call]
# Recorded sex cannot be checked for a Y-negative model; this variant drops the ones with a female-like X.
calls[, possible_female := y_call %in% "Y_negative" & !is.na(x_ratio) & x_ratio >= np$x_ratio_flag]
calls[, y_call_x_checked := fifelse(possible_female, NA_character_, y_call)]
# Continuous score: share of the Y-positive male level that is missing, median over the core genes.
expressed_males <- calls[recorded_male == TRUE & rna_call %in% "expressed", ModelID]
fraction <- y_fraction(log2((2^log_expr[with_rna, ] - 1) + 0.1), male_ids = expressed_males, floor_sd = 0)
calls[, y_loss_fraction := fraction$f_raw[match(ModelID, fraction$sample)]]
assert_unique(calls$ModelID, "nextgen_y_calls")
write_ext(calls, "nextgen_y_calls")
all_males <- calls[recorded_male == TRUE]
print(all_males[, .N, by = .(rna_call, dna_call)][order(rna_call, dna_call)])
print(all_males[, .N, by = .(y_tier, y_call)][order(y_tier)])

# --- donors -------------------------------------------------------------------------------------------
# Shared with the Sanger organoids: the same HCMI donor ID as a screened Sanger organoid.
# Inside NextGen: the same HCMI donor ID, the same name root, or the same genotype.
sanger <- read_ext("organoid_cohort", dir = ext_dir)
sanger_donor <- hcmi_donor(sanger$ModelID)
similarity <- genotype_similarity(damaging, np$genotype_rare_fraction)
cohort <- merge(screens[, .(ModelID, ScreenID, Library, ScreenType, OncotreeLineage, OncotreePrimaryDisease, OncotreeSubtype, CasActivity, ScreenDoublingTime, ScreenNNMD, ScreenROCAUC)], calls, by = "ModelID", sort = FALSE)
stopifnot(nrow(cohort) == nrow(screens))
cohort[, hcmi_donor := fcoalesce(hcmi_donor(name), hcmi_donor(alias))]
cohort[, name_root := sub(np$name_root_suffix, "", name)]
cohort[, shared_with_sanger := !is.na(hcmi_donor) & hcmi_donor %in% sanger_donor[!is.na(sanger_donor)]]
# The Sanger organoid of the same name where there is one (a donor can have several organoids there),
# otherwise the donor's organoid in the Sanger cohort, otherwise its first by name.
cohort[, sanger_model := partner_model(name, hcmi_donor, sanger$ModelID, sanger_donor, preferred = sanger$in_primary)]
cohort[shared_with_sanger == FALSE, sanger_model := NA_character_]
genotyped <- intersect(cohort$ModelID, rownames(similarity))
cohort[, donor := donor_groups(ModelID, keys = list(hcmi = hcmi_donor, root = name_root), similarity = similarity[genotyped, genotyped], cutoff = np$genotype_same_donor)]
setnames(cohort, "donor", "PatientID")

# --- cohort of screened models -----------------------------------------------------------------------
profiled <- intersect(rownames(damaging), rownames(hotspot))
cohort[, tp53_altered := fifelse(ModelID %in% profiled, pmax(flag_of(damaging, "TP53", ModelID), flag_of(hotspot, "TP53", ModelID)), NA_real_)]
cohort[, kras_hotspot := fifelse(ModelID %in% rownames(hotspot), flag_of(hotspot, "KRAS", ModelID), NA_real_)]
cohort <- merge(cohort, signatures[, .(ModelID, WGD, Aneuploidy, Ploidy, MSIScore)], by = "ModelID", all.x = TRUE, sort = FALSE)
cohort[, malignant := !OncotreePrimaryDisease %in% np$non_malignant]
cohort[, stratum := do.call(paste, c(.SD, sep = " | ")), .SDcols = np$stratum]
cohort[, dome := as.numeric(ScreenType == np$dome_format)]
cohort[, eligible := recorded_male & malignant & !shared_with_sanger]
cohort[, called := eligible & !is.na(y_call)]
# One model per donor among the called ones: the better screen.
cohort[, first_of_donor := FALSE]
cohort[called == TRUE, first_of_donor := best_of_group(ModelID, PatientID, ScreenNNMD)]
cohort[, in_primary := called & first_of_donor]
cohort[, y_negative := fifelse(y_call == "Y_negative", 1, fifelse(y_call == "Y_positive", 0, NA_real_))]
setorder(cohort, ModelID)
write_ext(cohort, "nextgen_cohort")

# Every pair of screened models read as one donor, and every model shared with the Sanger cohort.
pair_index <- which(outer(cohort$PatientID, cohort$PatientID, "==") & upper.tri(diag(nrow(cohort))), arr.ind = TRUE)
similarity_of <- function(a, b) ifelse(a %in% genotyped & b %in% genotyped, similarity[cbind(match(a, rownames(similarity)), match(b, rownames(similarity)))], NA_real_)
donor_pairs <- data.table(kind = "same donor inside NextGen", model = cohort$ModelID[pair_index[, 1]], other = cohort$ModelID[pair_index[, 2]])
donor_pairs[, c("name", "other_name") := .(cohort$name[match(model, cohort$ModelID)], cohort$name[match(other, cohort$ModelID)])]
donor_pairs[, same_hcmi_donor := cohort$hcmi_donor[match(model, cohort$ModelID)] == cohort$hcmi_donor[match(other, cohort$ModelID)]]
donor_pairs[is.na(same_hcmi_donor), same_hcmi_donor := FALSE]
donor_pairs[, same_name_root := cohort$name_root[match(model, cohort$ModelID)] == cohort$name_root[match(other, cohort$ModelID)]]
donor_pairs[, genotype_similarity := similarity_of(model, other)]
donor_pairs[, sex := cohort$Sex[match(model, cohort$ModelID)]]
donor_pairs[, c("called", "other_called") := .(cohort$called[match(model, cohort$ModelID)], cohort$called[match(other, cohort$ModelID)])]
# The model of the pair that enters the cohort, if either does.
donor_pairs[, kept := fifelse(cohort$in_primary[match(model, cohort$ModelID)], name, fifelse(cohort$in_primary[match(other, cohort$ModelID)], other_name, NA_character_))]
# For a donor screened at Sanger: every organoid Sanger screened of that donor, and the one in its cohort.
of_donor <- function(donor, rows = rep(TRUE, nrow(sanger))) unname(vapply(donor, function(one) paste(sanger$ModelID[which(sanger_donor == one & rows)], collapse = "; "), character(1)))
with_sanger <- cohort[shared_with_sanger == TRUE, .(kind = "donor screened at Sanger", model = ModelID, other = sanger_model, name, other_name = sanger_model, same_hcmi_donor = TRUE, same_name_root = FALSE, genotype_similarity = NA_real_, sex = Sex, called = NA, other_called = NA,
  kept = of_donor(hcmi_donor, sanger$in_primary), sanger_models_of_donor = of_donor(hcmi_donor))]
with_sanger[kept == "", kept := NA_character_]
donor_pairs[, sanger_models_of_donor := NA_character_]
donor_table <- rbind(donor_pairs, with_sanger)
stopifnot(all(with_sanger$other %in% sanger$ModelID), !anyNA(with_sanger$other))
write_ext(donor_table, "nextgen_donor_groups")
print(donor_table[, .(kind, name, other_name, same_hcmi_donor, same_name_root, genotype_similarity = round(genotype_similarity, 2), sex, kept, sanger_models_of_donor)])
# How well the genotype separates: named same-donor pairs against every other pair of screened models.
screened_similarity <- similarity[genotyped, genotyped]
named_pair <- outer(cohort$hcmi_donor[match(genotyped, cohort$ModelID)], cohort$hcmi_donor[match(genotyped, cohort$ModelID)], function(a, b) !is.na(a) & !is.na(b) & a == b) |
  outer(cohort$name_root[match(genotyped, cohort$ModelID)], cohort$name_root[match(genotyped, cohort$ModelID)], "==")
upper <- upper.tri(screened_similarity)
named_values <- screened_similarity[upper & named_pair]
other_values <- screened_similarity[upper & !named_pair]
log_msg("genotype similarity: %d named same-donor pairs score %.2f to %.2f; the other %d pairs at most %.2f below the cut of %.2f and %d at or above it",
  length(named_values), min(named_values), max(named_values), length(other_values), max(other_values[other_values < np$genotype_same_donor]), np$genotype_same_donor, sum(other_values >= np$genotype_same_donor))

screened_males <- cohort[recorded_male == TRUE]
print(screened_males[, .N, by = .(rna_call, dna_call)][order(rna_call, dna_call)])
print(screened_males[eligible == TRUE, .N, by = .(y_tier, y_call)][order(y_tier)])
primary <- cohort[in_primary == TRUE]
print(dcast(primary, OncotreeLineage ~ y_call, value.var = "ModelID", fun.aggregate = length))
print(dcast(primary, ScreenType ~ y_call, value.var = "ModelID", fun.aggregate = length))
log_msg("screened %d; recorded male %d, female %d, unknown %d; malignant males %d; not shared with Sanger %d; called %d (%d donors); primary cohort %d (%d Y-negative, %d Y-positive)",
  nrow(cohort), sum(cohort$Sex == "Male"), sum(cohort$Sex == "Female"), sum(cohort$Sex == "Unknown"), sum(cohort$recorded_male & cohort$malignant), sum(cohort$eligible), sum(cohort$called),
  uniqueN(cohort[called == TRUE, PatientID]), nrow(primary), sum(primary$y_negative == 1), sum(primary$y_negative == 0))

# --- the design, read without any score ----------------------------------------------------------------
# One stratum per lineage and screen format. A stratum holding one Y group says nothing about Y.
primary[, exposed := y_negative]
stratum_frame <- function(frame) copy(frame)[, OncotreeLineage := stratum]
pan_design <- build_design(stratum_frame(primary), tested = "exposed")
stopifnot(pan_design$estimable)
pan_unit <- unit_of(pan_design$design, "exposed")
strata <- primary[, .(n_negative = sum(y_negative == 1), n_positive = sum(y_negative == 0)), by = .(stratum, OncotreeLineage, ScreenType)][order(stratum)]
strata[, informative := n_negative > 0 & n_positive > 0]
print(strata)
log_msg("design: %d strata, %d with both Y groups; %d models in strata with one group; unit %.4f, worth %.1f against %.1f balanced models",
  nrow(strata), sum(strata$informative), strata[informative == FALSE, sum(n_negative + n_positive)], pan_unit, 2 / pan_unit^2, 2 / pan_unit^2)

# --- matrices: models in rows, genes in columns ------------------------------------------------------------
by_model <- function(path) {
  mat <- read_model_matrix(path)
  stopifnot(all(cohort$ScreenID %in% rownames(mat)))
  mat <- mat[cohort$ScreenID, ]
  rownames(mat) <- cohort$ModelID
  mat
}
effect <- by_model(study_path("screen_gene_effect.csv"))
dependency <- by_model(study_path("screen_gene_dependency.csv"))
stopifnot(setequal(colnames(effect), colnames(dependency)))
dependency <- dependency[, colnames(effect)]
saveRDS(effect, file.path(ext_out_dir, "nextgen_effect.rds"))
saveRDS(dependency, file.path(ext_out_dir, "nextgen_dependency.rds"))

# --- gene annotation -------------------------------------------------------------------------------------------
cell_genes <- read_result("gene_annotation")
organoid_genes <- read_ext("organoid_gene_annotation", dir = ext_dir)
gene_ann <- annotate_genes(colnames(effect))
assert_unique(gene_ann$feature, "NextGen genes")
at <- match(gene_ann$feature, cell_genes$feature)
gene_ann[, cell_line_in_scan := cell_genes$in_scan[at] %in% TRUE]
gene_ann[, cell_line_selective := cell_genes$selective[at] %in% TRUE]
gene_ann[, organoid_in_scan := feature %in% organoid_genes[in_scan == TRUE, feature]]
gene_ann[, n_missing := colSums(is.na(effect))]
# Scan genes: a score in every screened model, and not on chrY.
gene_ann[, in_scan := !is_y_gene & n_missing == 0]
dependent <- dependency[primary$ModelID, ] >= params$gene_scan$dependency_probability_cut
gene_ann[, n_dependent := colSums(dependent, na.rm = TRUE)]
gene_ann[, n_not_dependent := colSums(!dependent, na.rm = TRUE)]
gene_ann[, selective := in_scan & n_dependent >= np$selective_min & n_not_dependent >= np$selective_min]
gene_ann[, named := symbol %in% np$named_genes]
stopifnot(all(np$named_genes %in% gene_ann[in_scan == TRUE, symbol]), all(np$controls$gene %in% gene_ann[in_scan == TRUE, symbol]), !anyDuplicated(gene_ann[in_scan == TRUE, symbol]))
write_ext(gene_ann, "nextgen_gene_annotation")
essential_median <- apply(effect[, gene_ann[common_essential %in% TRUE & in_scan == TRUE, feature]], 1, median)
log_msg("genes: %d in the matrix, %d without a score in any model, %d with a score in some, %d chrY, %d in scan, %d selective; %d also in the cell-line scan, %d in the organoid scan; common-essential median %.3f (range %.3f to %.3f)",
  nrow(gene_ann), sum(gene_ann$n_missing == nrow(effect)), sum(gene_ann$n_missing > 0 & gene_ann$n_missing < nrow(effect)), sum(gene_ann$is_y_gene), sum(gene_ann$in_scan), sum(gene_ann$selective),
  sum(gene_ann$in_scan & gene_ann$cell_line_in_scan), sum(gene_ann$in_scan & gene_ann$organoid_in_scan), median(essential_median), min(essential_median), max(essential_median))

# --- the tests to come, counted before any fit --------------------------------------------------------------------
cells_of <- function(value) c(neg_in = sum(primary$y_negative == 1 & value), pos_in = sum(primary$y_negative == 0 & value), neg_out = sum(primary$y_negative == 1 & !value), pos_out = sum(primary$y_negative == 0 & !value))
test_counts <- rbindlist(lapply(sort(unique(primary$OncotreeLineage)), function(type) {
  cbind(data.table(class = "cancer_type", context = type, diseases = paste(sort(unique(primary[OncotreeLineage == type, OncotreePrimaryDisease])), collapse = "; ")), t(cells_of(primary$OncotreeLineage == type)))
}))
test_counts[, supported := pmin(neg_in, pos_in, neg_out, pos_out) >= sp$min_cell]
# Family of the three-source combination: genes in the cell-line scan and in at least one other scan.
combination_family <- sum(cell_genes$in_scan & (cell_genes$feature %in% gene_ann[in_scan == TRUE, feature] | cell_genes$feature %in% organoid_genes[in_scan == TRUE, feature]))
test_counts <- cbind(test_counts, data.table(
  primary_models = nrow(primary), pan_tests = sum(gene_ann$in_scan), selective_genes = sum(gene_ann$selective),
  named_genes = length(np$named_genes), combination_pan_genes = combination_family, pan_unit = pan_unit
))
write_ext(test_counts, "nextgen_test_counts")
print(test_counts[, .(context, neg_in, pos_in, neg_out, pos_out, supported)])

# --- cancer types to pool across sources, by rule ---------------------------------------------------------------------
# A type qualifies when no source has min_group models per Y group inside it but the sources together do.
cell_frame <- analysis_frame(read_result("cohort"), load_matrix("crispr_effect"))
sanger_primary <- sanger[in_primary == TRUE]
stopifnot(all(sanger_primary$OncotreeLineage %in% np$sanger_lineage$type))
sanger_primary[, lineage := np$sanger_lineage$lineage[match(OncotreeLineage, np$sanger_lineage$type)]]
count_types <- function(type, negative, source) data.table(type = type, negative = negative)[, .(n_negative = sum(negative == 1), n_positive = sum(negative == 0)), by = type][, source := source]
type_counts <- rbind(
  count_types(cell_frame$OncotreeLineage, cell_frame$exposed, "Cell lines"),
  count_types(sanger_primary$lineage, sanger_primary$y_negative, "Sanger organoids"),
  count_types(primary$OncotreeLineage, primary$y_negative, "NextGen")
)
pooled_types <- pooled_gate(type_counts, np$pooled$min_group)
type_counts[, tested_alone := n_negative >= np$pooled$min_group & n_positive >= np$pooled$min_group]
type_counts[, pooled := type %in% pooled_types]
setorder(type_counts, type, source)
write_ext(type_counts, "nextgen_pooled_types")
print(type_counts[type %in% union(pooled_types, primary$OncotreeLineage)])
log_msg("cancer types to pool: %s", if (length(pooled_types)) paste(pooled_types, collapse = ", ") else "none")
stopifnot(identical(pooled_types, np$pooled$expected_types))
# What the pooled test can reach: strata are the library in the cell lines and the screen format in NextGen.
cell_frame[, pooled_stratum := paste("Cell lines", fifelse(LibraryKY == 1, "KY", "no KY"), fifelse(LibraryHumagneCD == 1, "Humagne-CD", "no Humagne-CD"), sep = " | ")]
pooled_design <- rbindlist(lapply(pooled_types, function(type) {
  rows <- rbind(
    cell_frame[OncotreeLineage == type, .(stratum = pooled_stratum, negative = exposed)],
    sanger_primary[lineage == type, .(stratum = paste("Sanger organoids", library, sep = " | "), negative = y_negative)],
    primary[OncotreeLineage == type, .(stratum = paste("NextGen", ScreenType, sep = " | "), negative = y_negative)]
  )
  out <- rows[, .(n_negative = sum(negative == 1), n_positive = sum(negative == 0)), by = stratum][order(stratum)]
  out[, informative := n_negative > 0 & n_positive > 0]
  out[, arrangements := fifelse(informative, choose(n_negative + n_positive, n_negative), 1)]
  cbind(type = type, out)
}))
print(pooled_design)
pooled_floor <- pooled_design[, .(arrangements = prod(arrangements), n_negative = sum(n_negative[informative]), n_positive = sum(n_positive[informative])), by = type][, floor := 2 / arrangements]
print(pooled_floor)

# --- chance of replicating the four paired genes, before any Y test ----------------------------------------------------
# From the cell-line effect (both are Chronos scores), the design above and each gene's spread around its
# stratum means in this cohort. The spread is read without Y status in the model.
cell_pan <- read_result("gene_pan_results")
blind_design <- pan_design$design[, colnames(pan_design$design) != "exposed", drop = FALSE]
df_pan <- nrow(pan_design$design) - ncol(pan_design$design)
replication <- rbindlist(lapply(np$named_genes, function(symbol_name) {
  feature <- gene_ann[symbol == symbol_name & in_scan == TRUE, feature]
  spread <- sqrt(sum(qr.resid(qr(blind_design), effect[primary$ModelID, feature])^2) / (nrow(blind_design) - ncol(blind_design)))
  reference <- cell_pan[symbol == symbol_name, estimate]
  stopifnot(length(feature) == 1, length(reference) == 1)
  data.table(
    symbol = symbol_name, cell_line_estimate = reference, cell_line_se = cell_pan[symbol == symbol_name, (ci_high - ci_low) / (2 * qt(0.975, df))],
    nextgen_spread = spread, unit = pan_unit, expected_se = pan_unit * spread, df = df_pan,
    chance_at_strictest = analytic_power(reference, pan_unit * spread, df_pan, sp$fdr / length(np$named_genes)),
    chance_at_most_lenient = analytic_power(reference, pan_unit * spread, df_pan, sp$fdr)
  )
}))
write_ext(replication, "nextgen_replication_chance")
print(replication)

# --- counts expected from the exploration; a different input needs a new decision ------------------------------------------
expected <- np$expected_counts
found <- list(
  screened = nrow(cohort), males = sum(cohort$recorded_male), discordant_females = length(discordant_female),
  non_malignant_males = sum(cohort$recorded_male & !cohort$malignant), shared_with_sanger = sum(cohort$recorded_male & cohort$malignant & cohort$shared_with_sanger),
  second_of_donor = sum(cohort$called & !cohort$first_of_donor), set_aside = sum(cohort$eligible & !cohort$called),
  primary_negative = sum(primary$y_negative == 1), primary_positive = sum(primary$y_negative == 0)
)
agree <- vapply(names(expected), function(key) within_tolerance(found[[key]], expected[[key]]), logical(1))
print(data.table(count = names(expected), expected = unlist(expected), found = unlist(found[names(expected)]), within_tolerance = agree))
stopifnot(found$screened == expected$screened, all(agree))
stopifnot(
  all(primary$Sex == "Male"), all(primary$malignant), !any(primary$shared_with_sanger), !anyDuplicated(primary$PatientID), !anyNA(primary$y_negative),
  !any(primary$hcmi_donor %in% sanger_donor[!is.na(sanger_donor)])
)

# --- freeze ---------------------------------------------------------------------------------------------------------
summary_list <- list(
  screened = nrow(cohort), males = sum(cohort$recorded_male), females = sum(cohort$Sex == "Female"), unknown_sex = sum(cohort$Sex == "Unknown"),
  eligible = sum(cohort$eligible), called = sum(cohort$called), primary = nrow(primary), primary_negative = sum(primary$y_negative == 1), primary_positive = sum(primary$y_negative == 0),
  counts_found = found, counts_expected = expected,
  rna_thresholds = list(
    female_reference_n = length(female_reference), male_n = length(rna_males), discordant_females = discordant_female, female_mean_max = rna$thresholds$female_mean_max,
    depleted_mean_max = rna$thresholds$depleted_mean_max, expressed_mean_min = rna$thresholds$expressed_mean_min, male_reference = rna$thresholds$male_reference,
    expressed_fraction = expressed_fraction, ceiling_fraction = np$rna_ceiling_fraction, gene_max = as.list(rna$thresholds$gene_max), female_gene_max = as.list(rna$thresholds$female_gene_max)
  ),
  all_males = list(n = nrow(all_males), tiers = as.list(table(all_males$y_tier)), rna_by_dna = all_males[, .N, by = .(rna_call, dna_call)]),
  screened_male_tiers = as.list(table(screened_males[eligible == TRUE, y_tier])),
  possible_female_in_primary = sum(primary$possible_female),
  genotype = list(named_pairs = length(named_values), named_min = min(named_values), named_max = max(named_values), other_max_below_cut = max(other_values[other_values < np$genotype_same_donor]), other_at_or_above_cut = sum(other_values >= np$genotype_same_donor)),
  design = list(strata = nrow(strata), informative_strata = sum(strata$informative), models_in_one_group_strata = strata[informative == FALSE, sum(n_negative + n_positive)], unit = pan_unit, balanced_equivalent = 2 / pan_unit^2, df = df_pan),
  genes_in_matrix = nrow(gene_ann), genes_without_score = sum(gene_ann$n_missing == nrow(effect)), genes_in_scan = sum(gene_ann$in_scan), genes_selective = sum(gene_ann$selective),
  combination_pan_genes = combination_family, pooled_types = pooled_types, pooled_floor = pooled_floor,
  parameters_sha256 = ext_section_hash(c("shared", "nextgen"))
)
jsonlite::write_json(summary_list, file.path(ext_out_dir, "nextgen_prepare_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
frozen <- function(name) file.path(ext_out_dir, name)
write_extension_freeze(
  "nextgen",
  paths = c(
    cohort = frozen("nextgen_cohort.csv"), y_calls = frozen("nextgen_y_calls.csv"), gene_annotation = frozen("nextgen_gene_annotation.csv"),
    effect = frozen("nextgen_effect.rds"), dependency = frozen("nextgen_dependency.rds"), pooled_types = frozen("nextgen_pooled_types.csv"),
    replication_chance = frozen("nextgen_replication_chance.csv"), test_counts = frozen("nextgen_test_counts.csv")
  ),
  sections = c("shared", "nextgen"),
  counts = list(primary_negative = sum(primary$y_negative == 1), primary_positive = sum(primary$y_negative == 0), pan_tests = sum(gene_ann$in_scan), selective_genes = sum(gene_ann$selective), pooled_types = pooled_types),
  replace = "--replace" %in% commandArgs(trailingOnly = TRUE)
)
log_msg("NextGen preparation done")
