# Stage 02: collect every piece of Y-chromosome evidence per model from the raw tables.
# No calls are made here and no outcome table (CRISPR, drugs) is read.
# Output: derived/loy_reanalysis_v1/y_evidence.csv, one row per model in Model.csv.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))

core <- params$y_genes_core
dna <- params$dna_call

model <- fread(input_path("model"), select = c(
  "ModelID", "PatientID", "StrippedCellLineName", "Sex", "OncotreeLineage",
  "OncotreePrimaryDisease", "OncotreeSubtype", "EngineeredModel", "ModelType"
))
assert_unique(model$ModelID, "Model.csv")

# Chromosome and pseudoautosomal flag from the HGNC cytoband ("Xp22.33 and Yp11.2" marks a PAR gene).
gene <- fread(input_path("gene"), select = c("symbol", "location"))
gene[, chromosome := fifelse(grepl("^([0-9]{1,2}|X|Y)[pq]", location), sub("^([0-9]{1,2}|X|Y)[pq].*$", "\\1", location), NA_character_)]
gene[, par := grepl("^X", location) & grepl(" and Y", location)]
x_genes <- gene[chromosome == "X" & par == FALSE, symbol]
log_msg("gene annotation: %d X non-PAR genes, %d PAR genes, %d Y genes", length(x_genes), sum(gene$par), sum(gene$chromosome == "Y", na.rm = TRUE))

symbol_matrix <- function(mat, symbols, prefix) {
  colnames(mat) <- gene_symbol(colnames(mat))
  stopifnot(all(symbols %in% colnames(mat)), !anyDuplicated(colnames(mat)[colnames(mat) %in% symbols]))
  out <- as.data.table(mat[, symbols, drop = FALSE], keep.rownames = "ModelID")
  setnames(out, symbols, paste0(prefix, symbols))
  out
}

# --- RNA ---------------------------------------------------------------------
rna <- read_default_profiles(input_path("expr_tpm_pc"), features = core)
audit_ids("expression default profiles", rownames(rna), model$ModelID)
rna_dt <- symbol_matrix(rna, core, "rna_")
rna_dt[, rna_y_mean := rowMeans(.SD), .SDcols = paste0("rna_", core)]

xist <- read_default_profiles(input_path("expr_tpm_all"), features = "XIST")
stopifnot(ncol(xist) == 1)
rna_dt <- merge(rna_dt, data.table(ModelID = rownames(xist), rna_XIST = xist[, 1]), by = "ModelID", all.x = TRUE)

counts <- read_default_profiles(input_path("expr_counts"))
rna_dt <- merge(
  rna_dt,
  data.table(ModelID = rownames(counts), rna_library_size = rowSums(counts, na.rm = TRUE)),
  by = "ModelID", all.x = TRUE
)
rm(counts)
invisible(gc())
log_msg("RNA evidence for %d models", nrow(rna_dt))

# --- DNA, WGS gene level -----------------------------------------------------
wgs <- read_default_profiles(input_path("cn_wgs"), features = c(core, x_genes))
audit_ids("WGS default profiles", rownames(wgs), model$ModelID)
wgs_dt <- symbol_matrix(wgs, core, "wgs_")
x_cols <- gene_symbol(colnames(wgs)) %in% x_genes
wgs_dt[, wgs_x_ratio := apply(wgs[, x_cols, drop = FALSE], 1, median, na.rm = TRUE)]
log_msg("WGS gene evidence for %d models (%d X genes used for the X ratio)", nrow(wgs_dt), sum(x_cols))

# --- DNA, WGS segments: fraction of evaluable non-PAR Y at the loss floor -----
seg <- fread(
  input_path("cn_segments"),
  select = c("ModelID", "IsDefaultEntryForModel", "Chrom", "Start", "End", "RelCopyRatio"),
  showProgress = FALSE
)
seg <- seg[Chrom == "chrY" & IsDefaultEntryForModel == "Yes"]
bin_size <- dna$bin_size
origin <- floor(dna$nonpar_start / bin_size) * bin_size
n_bins <- ceiling((dna$nonpar_end - origin) / bin_size)
seg_models <- unique(seg$ModelID)
bins <- matrix(NA_real_, length(seg_models), n_bins, dimnames = list(seg_models, NULL))
seg[, first_bin := pmax(1L, as.integer(floor((Start - origin) / bin_size)) + 1L)]
seg[, last_bin := pmin(n_bins, as.integer(ceiling((End - origin) / bin_size)))]
for (i in which(seg$first_bin <= seg$last_bin)) {
  bins[seg$ModelID[i], seg$first_bin[i]:seg$last_bin[i]] <- seg$RelCopyRatio[i]
}
# Mask windows where females, who have no Y, still show signal: mapping artefacts.
female_wgs <- intersect(model[Sex == "Female", ModelID], seg_models)
female_signal <- colMeans(bins[female_wgs, , drop = FALSE] > dna$loss_ratio_max, na.rm = TRUE)
mask <- !is.na(female_signal) & female_signal >= dna$female_mask_fraction
runs <- rle(mask)
run_end <- cumsum(runs$lengths)
run_start <- run_end - runs$lengths + 1
mask_windows <- data.table(
  start = origin + (run_start[runs$values] - 1) * bin_size,
  end = origin + run_end[runs$values] * bin_size
)
write_result(mask_windows, "y_female_artefact_windows")
log_msg("female-derived mask: %d windows, %.0f kb of %.0f kb", nrow(mask_windows), sum(mask) * bin_size / 1e3, n_bins * bin_size / 1e3)
evaluable <- bins[, !mask, drop = FALSE]
segment_dt <- data.table(
  ModelID = seg_models,
  wgs_y_covered_fraction = rowMeans(!is.na(evaluable)),
  wgs_y_floor_fraction = rowMeans(evaluable <= dna$loss_ratio_max, na.rm = TRUE),
  wgs_y_retained_fraction = rowMeans(evaluable >= dna$retained_ratio_min, na.rm = TRUE),
  wgs_y_segment_mean = rowMeans(evaluable, na.rm = TRUE)
)
wgs_dt <- merge(wgs_dt, segment_dt, by = "ModelID", all.x = TRUE)
rm(bins, evaluable)

# --- DNA, WES gene level (keyed by model condition) --------------------------
wes_header <- read_header(input_path("cn_wes"))
wes_keep <- wes_header[gene_symbol(wes_header) %in% c(core, x_genes)]
wes <- fread(input_path("cn_wes"), select = c("ModelConditionID", "IsDefaultEntryForMC", wes_keep), showProgress = FALSE)
# The default flag is "Yes" or "No_CDS-..." here, never "No".
wes <- wes[IsDefaultEntryForMC == "Yes"]
condition <- fread(input_path("model_condition"), select = c("ModelConditionID", "ModelID"))
audit_ids("WES model conditions", wes$ModelConditionID, condition$ModelConditionID)
wes <- merge(wes, condition, by = "ModelConditionID")
wes_core <- wes_keep[gene_symbol(wes_keep) %in% core]
stopifnot(length(wes_core) == length(core))
wes_x <- setdiff(wes_keep, wes_core)
wes[, wes_x_ratio := apply(as.matrix(.SD), 1, median, na.rm = TRUE), .SDcols = wes_x]
# A model with several default conditions takes the median across them.
wes_dt <- wes[, c(
  lapply(.SD, median, na.rm = TRUE),
  list(wes_n_conditions = .N)
), by = ModelID, .SDcols = c(wes_core, "wes_x_ratio")]
setnames(wes_dt, wes_core, paste0("wes_", gene_symbol(wes_core)))
log_msg("WES gene evidence for %d models (%d with more than one default condition)", nrow(wes_dt), sum(wes_dt$wes_n_conditions > 1))

# --- ploidy and genome doubling ----------------------------------------------
signatures <- fread(input_path("signatures"))
signatures <- signatures[IsDefaultEntryForModel == "Yes", .(ModelID, Ploidy, WGD)]
assert_unique(signatures$ModelID, "OmicsGlobalSignatures default profiles")

# --- assemble ----------------------------------------------------------------
evidence <- Reduce(
  function(left, right) merge(left, right, by = "ModelID", all.x = TRUE, sort = FALSE),
  list(model, rna_dt, wgs_dt, wes_dt, signatures)
)
evidence[, has_rna := !is.na(rna_y_mean)]
evidence[, has_wgs := !is.na(wgs_DDX3Y) | !is.na(wgs_ZFY)]
evidence[, has_wes := !is.na(wes_DDX3Y) | !is.na(wes_ZFY)]
stopifnot(nrow(evidence) == nrow(model))
assert_unique(evidence$ModelID, "y_evidence")
write_result(evidence, "y_evidence")
print(evidence[, .(models = .N, rna = sum(has_rna), wgs = sum(has_wgs), wes = sum(has_wes), any_dna = sum(has_wgs | has_wes)), by = Sex])
