# Stage 29: BeatAML2 (Bottomly et al. 2022), preparation.
# Primary AML samples are mixtures, so Y loss is measured as a fraction: the
# seven core Y genes, each relative to its X partner and to the median male
# specimen. Values inside the noise of the male bulk are set to 0. The karyotype
# gives the check (cells in clones without a Y). One specimen per patient.
# No drug response is compared with Y loss here. Outputs: extension/
# beataml_cohort, beataml_karyotypes, beataml_test_counts, beataml_families,
# beataml_compound_match, beataml_input_manifest, beataml_auc.rds,
# beataml_prepare_summary.json, extension_freeze_beataml.json.
# Pass "--replace" to archive an existing freeze and write a new one.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
suppressPackageStartupMessages(library(readxl))
assert_freeze()
dir.create(ext_out_dir, showWarnings = FALSE, recursive = TRUE)

bp <- ext_params$beataml
sp <- ext_params$shared
y_genes <- sp$y_genes
x_partners <- sp$x_partners
data_dir <- file.path(raw_dir, "drug_response", "beataml2")
within_tolerance <- function(found, expected) abs(found - expected) <= sp$count_tolerance

# --- inputs --------------------------------------------------------------------------------
inputs <- as.data.table(ext_params$inputs)[grepl("beataml2", file)]
manifest <- rbindlist(lapply(seq_len(nrow(inputs)), function(i) verify_download(file.path(root_dir, inputs$file[i]), inputs$checksum[i], inputs$algorithm[i])))
manifest[, source := inputs$source]
write_ext(manifest, "beataml_input_manifest")

clinical <- as.data.table(suppressMessages(read_excel(file.path(data_dir, "beataml_wv1to4_clinical.xlsx"), guess_max = 5000)))
mapping <- as.data.table(suppressMessages(read_excel(file.path(data_dir, "beataml_waves1to4_sample_mapping.xlsx"))))
fits <- fread(file.path(data_dir, "beataml_probit_curve_fits_v4_dbgap.txt"), showProgress = FALSE)
mutations <- fread(file.path(data_dir, "beataml_wes_wv1to4_mutations_dbgap.txt"), select = c("dbgap_sample_id", "symbol"), showProgress = FALSE)
expression <- fread(file.path(data_dir, "beataml_waves1to4_norm_exp_dbgap.txt"), showProgress = FALSE)
wanted <- c(y_genes, x_partners, "XIST")
stopifnot(all(wanted %in% expression$display_label), !anyDuplicated(expression[display_label %in% wanted, display_label]))
levels_all <- t(as.matrix(expression[match(wanted, display_label), -(1:4)]))
colnames(levels_all) <- wanted
rm(expression)

# --- specimens: tumour RNA-seq that the study kept, with clinical annotation ---------------------
kept_rna <- mapping[rna_control %in% "No" & rna_include_in_analysis %in% "yes", dbgap_rnaseq_sample]
specimens <- clinical[!is.na(dbgap_rnaseq_sample) & dbgap_rnaseq_sample %in% rownames(levels_all)]
assert_unique(specimens$dbgap_rnaseq_sample, "BeatAML2 specimens with RNA-seq")
setnames(specimens, c("dbgap_rnaseq_sample", "dbgap_subject_id", "consensus_sex"), c("SpecimenID", "PatientID", "Sex"))
specimens[, PatientID := as.character(PatientID)]
with_drug <- unique(fits[!is.na(dbgap_rnaseq_sample) & dbgap_rnaseq_sample != "", dbgap_rnaseq_sample])
found <- list(
  male_with_rna = specimens[Sex %in% "Male", .N], male_with_rna_and_drug = specimens[Sex %in% "Male" & SpecimenID %in% with_drug, .N],
  patients = uniqueN(specimens[Sex %in% "Male" & SpecimenID %in% with_drug, PatientID])
)
specimens[, rna_kept_by_study := SpecimenID %in% kept_rna]

# --- karyotype: cells in clones without a Y, female clones, donor chimerism -------------------------
read_karyotypes <- parse_karyotype_y(specimens$karyotype)
read_karyotypes[, SpecimenID := specimens$SpecimenID]
setcolorder(read_karyotypes, "SpecimenID")
write_ext(read_karyotypes, "beataml_karyotypes")
specimens[, karyotype_parsed := read_karyotypes$parsed]
specimens[, karyotype_y_fraction := read_karyotypes$fraction]
specimens[, karyotype_mentions_y_loss := read_karyotypes$mentions_y_loss %in% TRUE]
specimens[, karyotype_female_clone := read_karyotypes$female_clone %in% TRUE]
specimens[, karyotype_chimerism := read_karyotypes$chimerism %in% TRUE]
unread <- specimens[Sex %in% "Male" & !is.na(karyotype) & !karyotype_parsed, .(SpecimenID, karyotype)]
log_msg("male karyotypes: %d given, %d read, %d not read; %d mention -Y; %d with a female clone or donor chimerism",
  specimens[Sex %in% "Male" & !is.na(karyotype), .N], specimens[Sex %in% "Male" & karyotype_parsed, .N], nrow(unread),
  specimens[Sex %in% "Male" & karyotype_mentions_y_loss, .N], specimens[Sex %in% "Male" & (karyotype_female_clone | karyotype_chimerism), .N])
if (nrow(unread) > 0) print(unread[, .(SpecimenID, karyotype = substr(karyotype, 1, 70))][1:min(.N, 20)])

# --- guards: who can serve as male or as the female reference ------------------------------------------
xist <- levels_all[specimens$SpecimenID, "XIST"]
is_male <- specimens$Sex %in% "Male"
is_female <- specimens$Sex %in% "Female"
xist_midpoint <- (median(xist[is_male]) + median(xist[is_female])) / 2
specimens[, xist := xist]
specimens[, xist_female_cells := is_male & xist > xist_midpoint]
# A recorded female with a Y in the karyotype, or with male-level XIST, does not define the female level.
specimens[, karyotype_has_y := karyotype_parsed & !karyotype_female_clone & (karyotype_y_fraction < 1) %in% TRUE & grepl("Y", karyotype)]
specimens[, female_reference := is_female & rna_kept_by_study & xist > xist_midpoint & !karyotype_has_y]
specimens[, eligible := is_male & rna_kept_by_study & !xist_female_cells & !karyotype_female_clone & !karyotype_chimerism]
log_msg("males with RNA-seq %d; left out: %d not kept by the study, %d with XIST above the midpoint (%.2f), %d with a female clone or chimerism; eligible %d; female reference %d",
  sum(is_male), sum(is_male & !specimens$rna_kept_by_study), sum(specimens$xist_female_cells), xist_midpoint, sum(is_male & (specimens$karyotype_female_clone | specimens$karyotype_chimerism)),
  sum(specimens$eligible), sum(specimens$female_reference))

# --- Y-loss fraction ------------------------------------------------------------------------------------
y_level <- levels_all[specimens$SpecimenID, y_genes]
ratio_level <- y_level - levels_all[specimens$SpecimenID, x_partners]
colnames(ratio_level) <- y_genes
males <- specimens[eligible == TRUE, SpecimenID]
females <- specimens[female_reference == TRUE, SpecimenID]
primary_fraction <- y_fraction(ratio_level, male_ids = males, female_ids = females, floor_sd = bp$floor_sd)
y_only_fraction <- y_fraction(y_level, male_ids = males, female_ids = females, floor_sd = bp$floor_sd)
specimens[, y_ratio := primary_fraction$ratio]
specimens[, f_no_floor := primary_fraction$f_raw]
specimens[, f := primary_fraction$f]
specimens[, f_y_only := y_only_fraction$f]
specimens[, f_floor := primary_fraction$floor[1]]
specimens[Sex %in% "Female" | !eligible, c("f", "f_no_floor", "f_y_only") := NA_real_]
log_msg("Y-loss fraction: floor %.3f (Y relative to X partner; %.3f for the Y genes alone); eligible males with f >= %.1f: %d; with f > 0: %d",
  primary_fraction$floor[1], y_only_fraction$floor[1], bp$high_loss_fraction, specimens[eligible == TRUE & f >= bp$high_loss_fraction, .N], specimens[eligible == TRUE & f > 0, .N])
specimens[, evidence_tier := fifelse(!karyotype_parsed | is.na(karyotype_y_fraction), "no karyotype", fifelse(karyotype_y_fraction > 0, "karyotype -Y", "karyotype XY"))]
agreement <- specimens[eligible == TRUE & karyotype_mentions_y_loss & !is.na(karyotype_y_fraction)]
log_msg("f against the karyotype fraction, %d eligible males with -Y in the karyotype: rank correlation %.2f", nrow(agreement), cor(agreement$f_no_floor, agreement$karyotype_y_fraction, method = "spearman"))
print(specimens[eligible == TRUE & f > 0, .N, by = .(evidence_tier, high = f >= bp$high_loss_fraction)][order(evidence_tier, high)])

# --- covariates ------------------------------------------------------------------------------------------
number <- function(x) suppressWarnings(as.numeric(as.character(x)))
tp53_samples <- mutations[symbol == "TP53", unique(dbgap_sample_id)]
exome_samples <- unique(mutations$dbgap_sample_id)
specimens[, specimen_type := specimenType]
specimens[, stage := fifelse(is.na(diseaseStageAtSpecimenCollection), "Unknown", diseaseStageAtSpecimenCollection)]
specimens[, age := fifelse(is.na(number(ageAtSpecimenAcquisition)), number(ageAtDiagnosis), number(ageAtSpecimenAcquisition))]
specimens[, wave := cohort]
specimens[, cbf_fusion := as.numeric(consensusAMLFusions %in% bp$cbf_fusions)]
specimens[, t8_21 := as.numeric(consensusAMLFusions %in% "RUNX1-RUNX1T1")]
specimens[, flt3_itd := fifelse(`FLT3-ITD` == "positive", 1, fifelse(`FLT3-ITD` == "negative", 0, NA_real_))]
specimens[, npm1 := fifelse(NPM1 == "positive", 1, fifelse(NPM1 == "negative", 0, NA_real_))]
# TP53 from the exome calls; a specimen without exome data is unknown, never wild-type.
specimens[, tp53 := fifelse(dbgap_dnaseq_sample %in% tp53_samples, 1, fifelse(exomeSeq %in% "y" | dbgap_dnaseq_sample %in% exome_samples, 0, NA_real_))]
specimens[, blasts_bm := number(`%.Blasts.in.BM`)]
specimens[, blasts_pb := number(`%.Blasts.in.PB`)]
specimens[, blasts := fifelse(specimen_type == "Bone Marrow Aspirate", blasts_bm, blasts_pb)]
specimens[, monocytic := fifelse(is.na(fabBlastMorphology), NA_real_, as.numeric(fabBlastMorphology %in% bp$monocytic_fab))]
specimens[, collection_day := number(timeOfSampleCollectionRelativeToInclusion)]
specimens[, has_drug := SpecimenID %in% with_drug]

# --- one specimen per patient: initial diagnosis first, then the earliest --------------------------------------
candidates <- specimens[eligible == TRUE & has_drug == TRUE & !is.na(f)]
candidates[, later_stage := stage != "Initial Diagnosis"]
setorder(candidates, PatientID, later_stage, collection_day, SpecimenID, na.last = TRUE)
first_specimen <- candidates[!duplicated(PatientID), SpecimenID]
specimens[, in_primary := SpecimenID %in% first_specimen]
cohort_columns <- c(
  "SpecimenID", "PatientID", "Sex", "eligible", "has_drug", "in_primary", "female_reference", "rna_kept_by_study", "xist", "xist_female_cells",
  "y_ratio", "f", "f_no_floor", "f_y_only", "f_floor", "evidence_tier", "karyotype", "karyotype_parsed", "karyotype_y_fraction", "karyotype_mentions_y_loss", "karyotype_female_clone", "karyotype_chimerism",
  "specimen_type", "stage", "age", "wave", "cbf_fusion", "t8_21", "consensusAMLFusions", "flt3_itd", "npm1", "tp53", "blasts_bm", "blasts_pb", "blasts", "monocytic", "fabBlastMorphology", "collection_day", "dbgap_dnaseq_sample"
)
beataml_cohort <- specimens[, cohort_columns, with = FALSE]
write_ext(beataml_cohort, "beataml_cohort")
primary <- beataml_cohort[in_primary == TRUE]
stopifnot(!anyDuplicated(primary$PatientID), all(primary$Sex == "Male"))
found$high_loss <- specimens[eligible == TRUE & has_drug == TRUE & f >= bp$high_loss_fraction, .N]
log_msg("eligible males with drug data: %d specimens, %d patients; primary cohort %d specimens (%d with f >= %.1f, %d with f > 0); stages: %s",
  nrow(candidates), uniqueN(candidates$PatientID), nrow(primary), sum(primary$f >= bp$high_loss_fraction), bp$high_loss_fraction, sum(primary$f > 0),
  paste(sprintf("%s %d", names(table(primary$stage)), table(primary$stage)), collapse = ", "))
print(primary[f > 0, .(SpecimenID, f = round(f, 2), karyotype_y_fraction = round(karyotype_y_fraction, 2), evidence_tier, stage, age, cbf_fusion, blasts, karyotype = substr(karyotype, 1, 40))][order(-f)])

# --- counts expected from the exploration ----------------------------------------------------------------------
expected <- bp$expected_counts
agree <- vapply(names(expected), function(key) within_tolerance(found[[key]], expected[[key]]), logical(1))
print(data.table(count = names(expected), expected = unlist(expected), found = unlist(found[names(expected)]), within_tolerance = agree))
stopifnot(all(agree))

# --- drug response: specimens in rows, inhibitors in columns (area under the fitted curve) ------------------------
fits <- fits[dbgap_rnaseq_sample %in% beataml_cohort[eligible == TRUE, SpecimenID] & is.finite(auc)]
stopifnot(!anyDuplicated(fits[, .(dbgap_rnaseq_sample, inhibitor)]))
wide <- dcast(fits, dbgap_rnaseq_sample ~ inhibitor, value.var = "auc")
auc <- as.matrix(wide[, -1])
rownames(auc) <- wide$dbgap_rnaseq_sample
saveRDS(auc, file.path(ext_out_dir, "beataml_auc.rds"))

# --- the tests to come, counted before any fit --------------------------------------------------------------------
primary_auc <- auc[primary$SpecimenID, , drop = FALSE]
exposure <- primary$f
test_counts <- data.table(
  inhibitor = colnames(primary_auc), specimens = colSums(!is.na(primary_auc)),
  high_loss = colSums(!is.na(primary_auc) & exposure >= bp$high_loss_fraction), any_loss = colSums(!is.na(primary_auc) & exposure > 0),
  exposure_sum_of_squares = apply(primary_auc, 2, function(v) sum((exposure[!is.na(v)] - mean(exposure[!is.na(v)]))^2))
)
test_counts[, tested := specimens >= bp$min_specimens & high_loss >= bp$min_high_loss]
write_ext(test_counts, "beataml_test_counts")
log_msg("inhibitors: %d with any data, %d pass the gate (at least %d specimens and %d with f >= %.1f); median specimens %d, median with high loss %d",
  nrow(test_counts), sum(test_counts$tested), bp$min_specimens, bp$min_high_loss, bp$high_loss_fraction, as.integer(median(test_counts[tested == TRUE, specimens])), as.integer(median(test_counts[tested == TRUE, high_loss])))

# --- drug families and the fixed list of inhibitors shared with the cell-line screens -------------------------------
families <- as.data.table(suppressMessages(read_excel(file.path(data_dir, "beataml_drug_families.xlsx"), sheet = "drug_family")))[inhibitor %in% colnames(auc), .(inhibitor, family)]
families <- unique(families)
write_ext(families, "beataml_families")
cell_line_drugs <- read_result("drug_pan_results")
cell_names <- unique(toupper(cell_line_drugs$compound_name[!is.na(cell_line_drugs$compound_name)]))
name_keys <- function(inhibitor) {
  base <- trimws(sub("\\s*\\(.*\\)\\s*$", "", inhibitor))
  inside <- trimws(unlist(regmatches(inhibitor, gregexpr("(?<=\\()[^)]*(?=\\))", inhibitor, perl = TRUE))))
  unique(toupper(c(base, inside)))
}
matches <- rbindlist(lapply(colnames(auc), function(inhibitor) {
  hit <- intersect(name_keys(inhibitor), cell_names)
  if (length(hit) == 0) return(NULL)
  data.table(inhibitor = inhibitor, compound_name = hit[1])
}))
write_ext(matches, "beataml_compound_match")
log_msg("inhibitors matched to a cell-line compound by name: %d of %d; with a drug family: %d", nrow(matches), ncol(auc), uniqueN(families$inhibitor))

# --- freeze ------------------------------------------------------------------------------------------------------------
jsonlite::write_json(list(
  specimens_with_rna = nrow(specimens), males = sum(is_male), eligible = sum(specimens$eligible), female_reference = sum(specimens$female_reference),
  eligible_with_drug = nrow(candidates), patients = uniqueN(candidates$PatientID), primary = nrow(primary),
  primary_high_loss = sum(primary$f >= bp$high_loss_fraction), primary_any_loss = sum(primary$f > 0),
  floor = primary_fraction$floor[1], xist_midpoint = xist_midpoint, counts_found = found, counts_expected = expected,
  karyotype_rank_correlation = cor(agreement$f_no_floor, agreement$karyotype_y_fraction, method = "spearman"), karyotype_compared = nrow(agreement),
  inhibitors = nrow(test_counts), inhibitors_tested = sum(test_counts$tested), matched_to_cell_lines = nrow(matches),
  parameters_sha256 = ext_section_hash(c("shared", "beataml"))
), file.path(ext_out_dir, "beataml_prepare_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
write_extension_freeze(
  "beataml",
  paths = c(
    cohort = file.path(ext_out_dir, "beataml_cohort.csv"), auc = file.path(ext_out_dir, "beataml_auc.rds"), test_counts = file.path(ext_out_dir, "beataml_test_counts.csv"),
    families = file.path(ext_out_dir, "beataml_families.csv"), compound_match = file.path(ext_out_dir, "beataml_compound_match.csv")
  ),
  sections = c("shared", "beataml"),
  counts = list(primary = nrow(primary), primary_high_loss = sum(primary$f >= bp$high_loss_fraction), inhibitors_tested = sum(test_counts$tested)),
  replace = "--replace" %in% commandArgs(trailingOnly = TRUE)
)
log_msg("BeatAML2 preparation done")
