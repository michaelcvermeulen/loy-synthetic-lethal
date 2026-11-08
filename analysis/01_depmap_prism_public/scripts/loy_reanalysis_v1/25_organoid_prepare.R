# Stage 25: Sanger organoids (Herranz-Ors et al. 2026), preparation.
# One screen per organoid, Y status (RNA-led with a DNA check), the cohort, the
# matrices, the gene annotation and the count of tests to come, then a freeze
# that the screen stage asserts. No dependency value is compared with Y status
# here. Outputs: extension/organoid_y_calls, organoid_cohort,
# organoid_gene_annotation, organoid_test_counts, organoid_input_manifest,
# organoid_effect.rds, organoid_depleted.rds, organoid_effect_v11.rds,
# organoid_prepare_summary.json, extension_freeze_organoid.json.
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

op <- ext_params$organoid
sp <- ext_params$shared
dna_par <- params$dna_call
core <- sp$y_genes
study_dir <- file.path(raw_dir, "studies", "sanger_organoid")
table_path <- function(i, extension) file.path(study_dir, "supplementary_tables", sprintf("supplementary_table_%d_revision.%s", i, extension))
figshare_path <- function(name) file.path(study_dir, "figshare", name)
read_sheet <- function(i, sheet) {
  out <- as.data.table(suppressMessages(read_excel(table_path(i, "xlsx"), sheet = sheet)))
  out[!is.na(sample_ID)]
}
within_tolerance <- function(found, expected) abs(found - expected) <= sp$count_tolerance

# --- inputs ----------------------------------------------------------------------------
inputs <- as.data.table(ext_params$inputs)[grepl("sanger_organoid", file)]
manifest <- rbindlist(lapply(seq_len(nrow(inputs)), function(i) verify_download(file.path(root_dir, inputs$file[i]), inputs$checksum[i], inputs$algorithm[i])))
manifest[, source := inputs$source]
write_ext(manifest, "organoid_input_manifest")

# --- one screen per organoid: the paper's rule, the better QC rank ---------------------------
qc <- read_sheet(4, "summary_QC_metrics")
setorder(qc, sample_ID, -order)
chosen <- qc[!duplicated(sample_ID)]
dual <- qc[, .N, by = sample_ID][N > 1, sample_ID]
# For the follow-up: the v1_1 screen wherever an organoid has both.
chosen_v11 <- qc[, .SD[if ("v1_1" %in% library) library == "v1_1" else seq_len(.N) == 1], by = sample_ID]
log_msg("screens: %d for %d organoids; %d organoids screened with both libraries; chosen %s", nrow(qc), nrow(chosen), length(dual), paste(sprintf("%s %d", names(table(chosen$library)), table(chosen$library)), collapse = ", "))

# --- sample annotation ------------------------------------------------------------------
clinical <- read_sheet(2, "clinical_data")
setnames(clinical, c("ICD-O-3 Code and Morphology", "organoid_derived_from_metastasis", "age_at_attendance"), c("morphology", "from_metastasis", "age"))
wgs <- read_sheet(3, "data_obtained_from_WGS")[model == "organoid"]
assert_unique(wgs$sample_ID, "organoid WGS rows")
assert_unique(clinical$sample_ID, "organoid clinical rows")
meta <- merge(
  clinical[, .(sample_ID, individual_ID, primary_tumour_type, sex, age = as.numeric(age), morphology, from_metastasis)],
  wgs[, .(sample_ID, genome_sex = gender, ploidy, wgd = wholeGenomeDuplication, ms_status = msStatus, PGA, purity)],
  by = "sample_ID", all.x = TRUE, sort = FALSE
)
meta[, male := sex %in% "Male" & genome_sex %in% "MALE"]
meta[, female_reference := sex %in% "Female" & genome_sex %in% "FEMALE"]
meta[, sex_discordant := !(male | female_reference)]
# ICD-O behaviour code /3 is malignant; adenomas (/0) are not cancer models.
meta[, malignant := is.na(morphology) | !grepl("/[012]( |$)", morphology)]

events <- read_sheet(3, "genomic_driver_events")[model == "organoid" & Genomic_alterations %in% op$mutation_types]
profiled <- wgs$sample_ID
altered <- function(symbol) fifelse(meta$sample_ID %in% events[gene == symbol, sample_ID], 1, fifelse(meta$sample_ID %in% profiled, 0, NA_real_))
for (symbol in op$genotype_genes) meta[, (paste0("alt_", symbol)) := altered(symbol)]

# --- RNA evidence and call -----------------------------------------------------------------
tpm <- as.data.table(suppressMessages(read_excel(table_path(3, "xlsx"), sheet = "RNAseq_TPM")))
stopifnot(all(core %in% tpm$gene), !anyDuplicated(tpm$gene[tpm$gene %in% core]))
y_tpm <- t(as.matrix(tpm[match(core, gene), -(1:2)]))
colnames(y_tpm) <- core
stopifnot(!anyNA(y_tpm), all(y_tpm >= 0))
log_expr <- log2(y_tpm + 1)
with_rna <- intersect(meta$sample_ID, rownames(log_expr))
thresholds_cell_line <- jsonlite::fromJSON(file.path(derived_dir, "call_thresholds.json"))
# The cell-line threshold of expressed_mean_min log2(TPM + 1), as a share of the typical male level in TPM.
expressed_fraction <- (2^params$rna_call$expressed_mean_min - 1) / (2^thresholds_cell_line$expressed_male_median - 1)
rna <- call_rna_y(
  log_expr[with_rna, ], female_ids = meta[female_reference & sample_ID %in% with_rna, sample_ID], male_ids = meta[male & sample_ID %in% with_rna, sample_ID],
  expressed_fraction = expressed_fraction, depleted_max_on = op$depleted_max_genes_on, expressed_min_on = op$expressed_min_genes_on
)
log_msg("RNA rule on the organoid scale: female reference n=%d, female mean max %.3f, male reference %.3f, expressed threshold %.3f (%.1f%% of the male level)",
  sum(meta$female_reference & meta$sample_ID %in% with_rna), rna$thresholds$depleted_mean_max, rna$thresholds$male_reference, rna$thresholds$expressed_mean_min, 100 * expressed_fraction)
lenient_call <- fifelse(rna$mean_y <= rna$thresholds$depleted_mean_max, "depleted", fifelse(rna$mean_y >= rna$thresholds$expressed_mean_min, "expressed", "intermediate"))

# --- DNA evidence and call: chrY copy number relative to ploidy over the non-PAR region ----------
segments <- fread(figshare_path("CNV-segments_organoids.tsv"), select = c("sample_ID", "chromosome", "start", "end", "copyNumber", "ploidy", "model"), showProgress = FALSE)
segments <- segments[model == "organoid" & chromosome == "chrY"]
segments[, ratio := copyNumber / ploidy]
dna <- dna_segment_summary(segments[, .(sample = sample_ID, start, end, ratio)], start = dna_par$nonpar_start, end = dna_par$nonpar_end, loss_max = dna_par$loss_ratio_max, retained_min = dna_par$retained_ratio_min)
dna[, dna_call := call_dna_segments(dna, floor_min = dna_par$wgs_floor_fraction_min)]
setnames(dna, c("sample", "covered", "floor", "retained", "ratio"), c("sample_ID", "dna_y_covered", "dna_y_floor", "dna_y_retained", "dna_y_ratio"))

# --- integration and alternative definitions ------------------------------------------------------
calls <- merge(meta, dna, by = "sample_ID", all.x = TRUE, sort = FALSE)
at <- match(calls$sample_ID, with_rna)
calls[, has_rna := !is.na(at)]
calls[, rna_y_mean := rna$mean_y[at]]
calls[, rna_n_genes_on := rna$n_on[at]]
calls[, rna_call := rna$call[at]]
calls[, rna_call_lenient := lenient_call[at]]
for (g in core) calls[, (paste0("rna_", g)) := log_expr[with_rna, g][at]]
primary <- integrate_y_calls(calls$rna_call, calls$dna_call, calls$male)
calls[, y_call := primary$call]
calls[, y_tier := primary$tier]
as_call <- function(rna_state) fifelse(rna_state == "depleted", "Y_negative", fifelse(rna_state == "expressed", "Y_positive", NA_character_))
calls[, y_call_rna := fifelse(male, as_call(rna_call), NA_character_)]
calls[, y_call_concordant := fifelse(y_tier %in% "concordant", y_call, NA_character_)]
calls[, y_call_lenient := integrate_y_calls(rna_call_lenient, dna_call, male)$call]
calls[, y_call_genome_male := integrate_y_calls(rna_call, dna_call, genome_sex %in% "MALE")$call]
inherited <- fread(file.path(root_dir, "derived", "sanger_organoid", "sex_call_rna_vs_wgs.csv"))
on_disk <- inherited[match(calls$sample_ID, inherited$sample_ID)]
calls[, y_call_on_disk := fifelse(on_disk$WGS_sex == "Male", fifelse(on_disk$loss_of_Y_expr, "Y_negative", "Y_positive"), NA_character_)]
# Continuous score: share of the Y-positive male level that is missing, median over the core genes.
expressed_males <- calls[male & rna_call %in% "expressed", sample_ID]
fraction <- y_fraction(log2(y_tpm[with_rna, ] + 0.1), male_ids = expressed_males, floor_sd = 0)
calls[, y_loss_fraction := fraction$f_raw[match(sample_ID, fraction$sample)]]
assert_unique(calls$sample_ID, "organoid_y_calls")
write_ext(calls, "organoid_y_calls")

# --- cohort of screened organoids -------------------------------------------------------------------
cohort <- merge(chosen[, .(sample_ID, library, TPR, AUC_ROC, COEFF, genes_depleted, SD_LFC, qc_order = order)], calls, by = "sample_ID", sort = FALSE)
stopifnot(nrow(cohort) == nrow(chosen), all(cohort$has_rna))
setnames(cohort, c("sample_ID", "primary_tumour_type", "sex"), c("ModelID", "OncotreeLineage", "Sex"))
cohort[, PatientID := sprintf("%.0f", individual_ID)]
cohort[, LibraryV1 := as.numeric(library == "v1_1")]
cohort[, MSI := as.numeric(ms_status == "MSI")]
cohort[, WGD := as.numeric(wgd)]
cohort[, tp53_altered := alt_TP53]
cohort[, metastasis := as.numeric(from_metastasis == "Yes")]
cohort[, eligible := male & malignant]
cohort[, called := eligible & !is.na(y_call)]
# One organoid per patient among the called ones: the better screen.
setorder(cohort, PatientID, -qc_order)
cohort[, first_of_patient := FALSE]
cohort[called == TRUE, first_of_patient := !duplicated(PatientID)]
cohort[, in_primary := called & first_of_patient]
cohort[, y_negative := fifelse(y_call == "Y_negative", 1, fifelse(y_call == "Y_positive", 0, NA_real_))]
setorder(cohort, ModelID)
write_ext(cohort, "organoid_cohort")

screened_males <- cohort[male == TRUE]
print(screened_males[, .N, by = .(rna_call, dna_call)][order(rna_call, dna_call)])
print(screened_males[, .N, by = .(y_tier, y_call)][order(y_tier)])
print(cohort[called == TRUE, .N, by = .(OncotreeLineage, y_call)][order(OncotreeLineage, y_call)])
print(cohort[in_primary == TRUE, .N, by = .(OncotreeLineage, y_call)][order(OncotreeLineage, y_call)])
log_msg("screened organoids %d; male by both sexes %d; malignant males %d; called %d (%d patients); primary cohort %d (%d Y-negative, %d Y-positive)",
  nrow(cohort), sum(cohort$male), sum(cohort$eligible), sum(cohort$called), uniqueN(cohort[called == TRUE, PatientID]),
  sum(cohort$in_primary), sum(cohort$in_primary & cohort$y_negative %in% 1), sum(cohort$in_primary & cohort$y_negative %in% 0))

# --- counts expected from the exploration; a different input needs a new decision -------------------------
expected <- op$expected_counts
contradictions <- screened_males[y_tier %in% "discordant", .N]
type_count <- function(type, call) cohort[called == TRUE & OncotreeLineage == type & y_call == call, .N]
found <- list(
  screened = nrow(cohort), males_concordant = nrow(screened_males),
  rna_depleted = screened_males[rna_call == "depleted", .N], rna_expressed = screened_males[rna_call == "expressed", .N], rna_intermediate = screened_males[rna_call == "intermediate", .N],
  dna_contradictions = contradictions,
  colorectal_negative = type_count("Colorectal", "Y_negative"), colorectal_positive = type_count("Colorectal", "Y_positive"),
  oesophageal_negative = type_count("Oesophageal", "Y_negative"), oesophageal_positive = type_count("Oesophageal", "Y_positive")
)
agree <- vapply(names(expected), function(key) within_tolerance(found[[key]], expected[[key]]), logical(1))
print(data.table(count = names(expected), expected = unlist(expected), found = unlist(found[names(expected)]), within_tolerance = agree))
stopifnot(found$screened == expected$screened, all(agree))
stopifnot(!any(cohort[in_primary == TRUE, Sex != "Male" | genome_sex != "MALE"]), !anyDuplicated(cohort[in_primary == TRUE, PatientID]))

# --- matrices: organoids in rows, genes in columns ------------------------------------------------------
to_matrix <- function(long, value, screens) {
  kept <- long[screens, on = .(sample_ID, library), nomatch = NULL]
  wide <- dcast(kept, sample_ID ~ gene, value.var = value)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$sample_ID
  mat
}
lfc <- fread(table_path(6, "csv"), sep = " ", showProgress = FALSE)
depletion <- fread(table_path(5, "csv"), sep = " ", showProgress = FALSE)
effect <- to_matrix(lfc, "LFC", chosen[, .(sample_ID, library)])
effect_v11 <- to_matrix(lfc, "LFC", chosen_v11[, .(sample_ID, library)])
depleted <- to_matrix(depletion, "is_depleted", chosen[, .(sample_ID, library)])
stopifnot(
  nrow(effect) == nrow(chosen), !anyNA(effect), !anyNA(depleted), !anyNA(effect_v11),
  identical(dimnames(effect), dimnames(depleted)), identical(dimnames(effect), dimnames(effect_v11)), all(depleted %in% c(0, 1)),
  setequal(rownames(effect), cohort$ModelID)
)
saveRDS(effect, file.path(ext_out_dir, "organoid_effect.rds"))
saveRDS(effect_v11, file.path(ext_out_dir, "organoid_effect_v11.rds"))
saveRDS(depleted, file.path(ext_out_dir, "organoid_depleted.rds"))
log_msg("matrices: %d organoids x %d genes; %d organoids differ in the v1_1 variant", nrow(effect), ncol(effect), sum(rowSums(effect != effect_v11) > 0))

# --- gene annotation ----------------------------------------------------------------------------------------
cell_line_genes <- read_result("gene_annotation")
position <- fread(figshare_path("position_genes.csv"))
gene_ann <- unique(lfc[, .(symbol = gene, ensembl = ensembl_id)])[match(colnames(effect), symbol)]
assert_unique(gene_ann$symbol, "organoid genes")
at <- match(gene_ann$symbol, cell_line_genes$symbol)
gene_ann[, feature := cell_line_genes$feature[at]]
gene_ann[, chromosome := sub("^chr", "", position$chromosome[match(symbol, position$gene)])]
gene_ann[is.na(chromosome), chromosome := cell_line_genes$chromosome[match(symbol, cell_line_genes$symbol)]]
gene_ann[, is_y_gene := chromosome %in% "Y" | cell_line_genes$is_y_gene[at] %in% TRUE]
gene_ann[, y_partner := params$xy_paralogs$y_gene[match(symbol, params$xy_paralogs$x_partner)]]
gene_ann[, common_essential := cell_line_genes$common_essential[at] %in% TRUE]
gene_ann[, cell_line_selective := cell_line_genes$selective[at] %in% TRUE]
gene_ann[, cell_line_in_scan := cell_line_genes$in_scan[at] %in% TRUE]
gene_ann[, in_scan := !is_y_gene]
primary_ids <- cohort[in_primary == TRUE, ModelID]
gene_ann[, n_depleted := colSums(depleted[primary_ids, ] == 1)]
gene_ann[, n_not_depleted := colSums(depleted[primary_ids, ] == 0)]
gene_ann[, selective := in_scan & n_depleted >= op$selective_min & n_not_depleted >= op$selective_min]
gene_ann[, seen_before := symbol %in% op$seen_before]
write_ext(gene_ann, "organoid_gene_annotation")
# Scale of the scores: Chronos puts the median common-essential gene at -1.
essential_median <- apply(effect[, gene_ann[common_essential == TRUE, symbol]], 1, median)
effect_scale <- -median(essential_median)
log_msg("genes: %d in the matrix, %d chrY removed, %d in scan, %d selective, %d also in the cell-line scan, %d cell-line selective; common-essential median %.3f (range %.3f to %.3f)",
  nrow(gene_ann), sum(gene_ann$is_y_gene), sum(gene_ann$in_scan), sum(gene_ann$selective), sum(gene_ann$cell_line_in_scan & gene_ann$in_scan), sum(gene_ann$cell_line_selective & gene_ann$in_scan),
  -effect_scale, min(essential_median), max(essential_median))

# --- the tests to come, counted before any fit -------------------------------------------------------------------
primary <- cohort[in_primary == TRUE]
two_types <- primary[OncotreeLineage %in% op$tumour_types]
cells_of <- function(frame, value) c(neg_in = sum(frame$y_negative == 1 & value %in% 1), pos_in = sum(frame$y_negative == 0 & value %in% 1), neg_out = sum(frame$y_negative == 1 & value %in% 0), pos_out = sum(frame$y_negative == 0 & value %in% 0))
panel <- c(setNames(paste0("alt_", op$genotype_genes), op$genotype_genes), WGD = "WGD", MSI = "MSI", `Metastasis-derived` = "metastasis")
context_rows <- rbindlist(c(
  lapply(op$tumour_types, function(type) cbind(data.table(class = "cancer_type", context = type), t(cells_of(two_types, as.numeric(two_types$OncotreeLineage == type))))),
  lapply(names(panel), function(label) cbind(data.table(class = "genotype", context = label), t(cells_of(primary, primary[[panel[[label]]]]))))
))
context_rows[, supported := pmin(neg_in, pos_in, neg_out, pos_out) >= sp$min_cell]
test_counts <- cbind(context_rows, data.table(
  primary_organoids = nrow(primary), pan_tests = sum(gene_ann$in_scan), selective_genes = sum(gene_ann$selective),
  combination_pan_genes = sum(gene_ann$cell_line_in_scan & gene_ann$in_scan), combination_type_genes = sum(gene_ann$cell_line_selective & gene_ann$in_scan)
))
write_ext(test_counts, "organoid_test_counts")
print(context_rows)

# --- freeze ------------------------------------------------------------------------------------------------------
summary_list <- list(
  screened = nrow(cohort), males = sum(cohort$male), eligible = sum(cohort$eligible), called = sum(cohort$called),
  primary = nrow(primary), primary_negative = sum(primary$y_negative == 1), primary_positive = sum(primary$y_negative == 0), patients = uniqueN(primary$PatientID),
  counts_found = found, counts_expected = expected,
  rna_thresholds = list(female_reference_n = sum(meta$female_reference & meta$sample_ID %in% with_rna), depleted_mean_max = rna$thresholds$depleted_mean_max,
    expressed_mean_min = rna$thresholds$expressed_mean_min, male_reference = rna$thresholds$male_reference, expressed_fraction = expressed_fraction,
    gene_max = as.list(rna$thresholds$gene_max)),
  tiers = as.list(table(screened_males$y_tier)),
  genes_in_scan = sum(gene_ann$in_scan), genes_selective = sum(gene_ann$selective), effect_scale = effect_scale,
  dual_library_organoids = length(dual), parameters_sha256 = ext_section_hash(c("shared", "organoid"))
)
jsonlite::write_json(summary_list, file.path(ext_out_dir, "organoid_prepare_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
write_extension_freeze(
  "organoid",
  paths = c(
    cohort = file.path(ext_out_dir, "organoid_cohort.csv"), y_calls = file.path(ext_out_dir, "organoid_y_calls.csv"),
    gene_annotation = file.path(ext_out_dir, "organoid_gene_annotation.csv"), effect = file.path(ext_out_dir, "organoid_effect.rds"),
    depleted = file.path(ext_out_dir, "organoid_depleted.rds"), effect_v11 = file.path(ext_out_dir, "organoid_effect_v11.rds")
  ),
  sections = c("shared", "organoid"),
  counts = list(primary_negative = sum(primary$y_negative == 1), primary_positive = sum(primary$y_negative == 0), pan_tests = sum(gene_ann$in_scan), selective_genes = sum(gene_ann$selective)),
  replace = "--replace" %in% commandArgs(trailingOnly = TRUE)
)
log_msg("organoid preparation done")
