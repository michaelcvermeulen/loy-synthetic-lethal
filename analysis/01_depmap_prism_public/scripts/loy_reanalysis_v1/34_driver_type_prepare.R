# Stage 34: driver status inside cancer types, preparation.
# One status per gene and model (altered with any of a driver-flagged mutation, an
# amplification or a deep deletion), the census of every gene-by-cancer-type pair at the
# three Oncotree levels, the list of splits that reach the gate with what each design can
# carry, the same for the Sanger organoids, a census for DepMap NextGen, the pairs to
# combine across sources and the gene families, then a freeze that stages 35 and 36
# assert. No dependency value is compared with Y status here. Outputs:
# extension/driver_type_events, driver_type_status, driver_type_organoid_status,
# driver_type_census, driver_type_census_summary, driver_type_splits,
# driver_type_split_lines, driver_type_organoid_splits, driver_type_organoid_split_lines,
# driver_type_pairs, driver_type_genes, driver_type_prepare_summary.json,
# extension_freeze_driver_type.json.
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
assert_extension_freeze("nextgen", dir = ext_dir)
dir.create(ext_out_dir, showWarnings = FALSE, recursive = TRUE)

dp <- ext_params$driver_type
sp <- ext_params$shared
op <- ext_params$organoid
ctx_par <- params$context
gate <- ctx_par$min_cell
stopifnot(gate == sp$min_cell)
within_tolerance <- function(found, expected) abs(found - expected) <= sp$count_tolerance
split_id_of <- function(level, type, gene) paste(level, gsub("[^A-Za-z0-9]+", "_", type), gene, sep = "__")

# --- cell lines: the frame of the context screen ----------------------------------------------
cohort <- read_result("cohort")
effect <- load_matrix("crispr_effect")
frame <- analysis_frame(cohort, effect)
rm(effect)
log_msg("cell lines: %d male CRISPR lines, %d Y-negative", nrow(frame), sum(frame$exposed == 1))

# --- cell lines: one status per gene, with the definitions of 05_build_cohort.R ---------------
as_flag <- function(x) toupper(as.character(x)) %in% c("TRUE", "Y", "YES", "1")
mut <- fread(input_path("mutations"), select = c(
  "ModelID", "IsDefaultEntryForModel", "HugoSymbol", "Hotspot", "OncogeneHighImpact", "TumorSuppressorHighImpact", "HessDriver"
), showProgress = FALSE)
mut <- mut[IsDefaultEntryForModel == "Yes"]
mut[, driver := as_flag(Hotspot) | as_flag(OncogeneHighImpact) | as_flag(TumorSuppressorHighImpact) | as_flag(HessDriver)]
mutation_pairs <- unique(mut[driver == TRUE & HugoSymbol != "" & ModelID %in% frame$ModelID, .(ModelID, gene = HugoSymbol)])
profiled <- fread(input_path("mutation_hotspot"), select = c("ModelID", "IsDefaultEntryForModel"))[IsDefaultEntryForModel == "Yes", unique(ModelID)]
rm(mut)

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
  from_wgs <- if (symbol %in% colnames(wgs)) wgs[match(frame$ModelID, rownames(wgs)), symbol] else rep(NA_real_, nrow(frame))
  from_wes <- if (symbol %in% names(wes)) wes[[symbol]][match(frame$ModelID, wes$ModelID)] else rep(NA_real_, nrow(frame))
  fifelse(!is.na(from_wgs), from_wgs, from_wes)
}
ratios <- sapply(cn_genes, copy_ratio)
rownames(ratios) <- frame$ModelID
pairs_where <- function(mask) {
  at <- which(mask, arr.ind = TRUE)
  data.table(ModelID = rownames(mask)[at[, 1]], gene = colnames(mask)[at[, 2]])
}
amp_ratios <- ratios[, ctx_par$amp_genes, drop = FALSE]
del_ratios <- ratios[, ctx_par$del_genes, drop = FALSE]
amplified <- pairs_where(!is.na(amp_ratios) & amp_ratios >= ctx_par$amp_ratio_min)
deleted <- pairs_where(!is.na(del_ratios) & del_ratios <= ctx_par$del_ratio_max)
copy_missing <- pairs_where(is.na(ratios))

# One row per altered model and gene, with the kinds of event that make it so.
event_table <- function(pairs_by_kind) {
  long <- rbindlist(lapply(names(pairs_by_kind), function(kind) cbind(pairs_by_kind[[kind]], kind = kind)))
  wide <- dcast(long, ModelID + gene ~ kind, fun.aggregate = length, value.var = "kind")
  for (kind in c("mutation", "amplification", "deletion")) {
    if (!kind %in% names(wide)) wide[, (kind) := 0L]
    wide[, (kind) := get(kind) > 0]
  }
  setcolorder(wide, c("ModelID", "gene", "mutation", "amplification", "deletion"))
  wide[]
}
events <- event_table(list(mutation = mutation_pairs, amplification = amplified, deletion = deleted))
genes <- sort(unique(events$gene))
status <- gene_status(frame$ModelID, genes, events, profiled, copy_missing)
log_msg("status: %d genes altered in at least one line; %d altered model-gene pairs (%d by mutation, %d by amplification, %d by deletion); %d missing",
  length(genes), nrow(events), sum(events$mutation), sum(events$amplification), sum(events$deletion), sum(is.na(status)))

# The rebuilt statuses must equal every alteration that the frozen context matrix kept.
context <- read_result("context_matrix")
definitions <- read_result("context_definitions")
stored <- definitions[kept == TRUE & class %in% c("driver", "amplification", "deletion")]
rebuilt <- function(class, symbol) {
  kind_events <- events[get(c(driver = "mutation", amplification = "amplification", deletion = "deletion")[[class]]) == TRUE]
  missing <- if (class == "driver") NULL else copy_missing
  drop(gene_status(frame$ModelID, symbol, kind_events, if (class == "driver") profiled else frame$ModelID, missing))
}
same_definition <- vapply(seq_len(nrow(stored)), function(i) {
  identical(unname(rebuilt(stored$class[i], stored$label[i])), as.numeric(context[[stored$context_id[i]]][match(frame$ModelID, context$ModelID)]))
}, logical(1))
if (!all(same_definition)) stop("rebuilt status differs from the stored context for: ", paste(stored$context_id[!same_definition], collapse = ", "))
log_msg("the rebuilt statuses equal all %d stored alteration contexts", nrow(stored))

write_ext(events, "driver_type_events")
write_ext(as.data.table(status, keep.rownames = "ModelID"), "driver_type_status")

# --- cell lines: census at the three levels -----------------------------------------------------
levels_used <- unlist(dp$levels)
census_cells <- rbindlist(lapply(names(levels_used), function(level) {
  out <- split_census(status, frame[[levels_used[[level]]]], frame$exposed, gate)
  out[, level := level]
}))
census_cells[, source := "Cell lines"]

# --- cell lines: the splits that reach the gate, and what each design can carry ----------------
covariates <- dp$covariates
splits <- census_cells[label %in% split_labels[1:2]]
splits[, split_id := split_id_of(level, type, gene)]
stopifnot(!anyDuplicated(splits$split_id))
members_of <- function(level, type) frame$ModelID[frame[[levels_used[[level]]]] %in% type]
members <- Map(members_of, splits$level, splits$type)
altered_of <- function(i) intersect(members[[i]], frame$ModelID[status[, splits$gene[i]] %in% 1])
altered <- lapply(seq_len(nrow(splits)), altered_of)
splits[, same_lines_as := same_lines_as(split_id, gene, members)]
splits[, fitted := is.na(same_lines_as)]
altered_key <- vapply(altered, function(lines) paste(sort(lines), collapse = ","), character(1))
splits[, altered_set := match(altered_key, unique(altered_key))]
set_first <- !duplicated(splits$altered_set)
set_group <- overlap_groups(altered[set_first], dp$shared_lines_fraction)
splits[, shared_group := set_group[match(altered_set, splits$altered_set[set_first])]]

design_facts <- function(i) {
  sframe <- split_frame(frame, frame$ModelID %in% members[[i]], status[, splits$gene[i]])
  designs <- split_designs(sframe, covariates, gate)
  alone <- build_design(sframe[ctx == 1], tested = "exposed", covariates = covariates)
  cells <- designs$cells
  donors <- function(mask) uniqueN(sframe$PatientID[mask])
  plain_inside <- function(design, tested) information_kept(design, tested, cells[["neg_in"]], cells[["pos_in"]])
  data.table(
    lines = nrow(sframe), residual_df = nrow(designs$inside$design) - ncol(designs$inside$design),
    donors_neg_alt = donors(sframe$exposed == 1 & sframe$ctx == 1), donors_pos_alt = donors(sframe$exposed == 0 & sframe$ctx == 1),
    donors_neg_wt = donors(sframe$exposed == 1 & sframe$ctx == 0), donors_pos_wt = donors(sframe$exposed == 0 & sframe$ctx == 0),
    kept_inside = plain_inside(designs$inside$design, "exposed_in"),
    kept_inside_alone = if (alone$estimable) plain_inside(alone$design, "exposed") else NA_real_,
    kept_specificity = if (designs$specificity && !is.null(designs$interaction)) sum(1 / cells) / solve(crossprod(designs$interaction$design))["interaction", "interaction"] else NA_real_,
    inside_unit = sqrt(solve(crossprod(designs$inside$design))["exposed_in", "exposed_in"]),
    specificity_fitted = designs$specificity && !is.null(designs$interaction)
  )
}
splits <- cbind(splits, rbindlist(lapply(seq_len(nrow(splits)), design_facts)))
# A split labelled as testable for specificity must have a design that can estimate it.
stopifnot(splits[label == split_labels[1], all(specificity_fitted)], splits[label == split_labels[2], !any(specificity_fitted)])
splits[, balanced_equivalent := 2 / inside_unit^2]
splits[, low_information := kept_inside < dp$information_flag]
splits[, exact_floor := exact_floor(neg_alt, pos_alt)]
splits[, robust_possible := exact_floor < sp$robust_permutation_p]
splits[, same_donor_twice := donors_neg_alt < neg_alt | donors_pos_alt < pos_alt]
# The same alteration in the other cancer types of that level: lines per Y group.
elsewhere <- function(i) {
  other <- !frame$ModelID %in% members[[i]] & status[, splits$gene[i]] %in% 1
  data.table(elsewhere_neg = sum(other & frame$exposed == 1), elsewhere_pos = sum(other & frame$exposed == 0))
}
splits <- cbind(splits, rbindlist(lapply(seq_len(nrow(splits)), elsewhere)))
setorder(splits, level, type, -neg_alt, -pos_alt, gene)
write_ext(splits, "driver_type_splits")

split_lines <- rbindlist(lapply(which(splits$fitted), function(i) {
  sframe <- split_frame(frame, frame[[levels_used[[splits$level[i]]]]] %in% splits$type[i], status[, splits$gene[i]])
  sframe[, .(split_id = splits$split_id[i], ModelID, exposed, ctx)]
}))
write_ext(split_lines, "driver_type_split_lines")
print(splits[fitted == TRUE, .(level, type, gene, neg_alt, pos_alt, neg_wt, pos_wt, label, kept_inside = round(kept_inside, 2), kept_alone = round(kept_inside_alone, 2), balanced = round(balanced_equivalent, 1), robust_possible, shared_group)], nrows = 100)

# --- Sanger organoids ----------------------------------------------------------------------------------
organoids <- read_ext("organoid_cohort", dir = ext_dir)[in_primary == TRUE]
organoids[, exposed := as.numeric(y_negative == 1)]
sheet <- as.data.table(suppressMessages(read_excel(
  file.path(raw_dir, "studies", "sanger_organoid", "supplementary_tables", "supplementary_table_3_revision.xlsx"), sheet = "genomic_driver_events"
)))
sheet <- sheet[!is.na(sample_ID) & model == "organoid"]
sheet <- sheet[sample_ID %in% organoids$ModelID]
listed <- sheet[Genomic_alterations %in% op$mutation_types]
is_deletion <- listed$Genomic_alterations %in% c("DEL", "HOM_DISRUPTION")
organoid_events <- event_table(list(
  mutation = listed[!is_deletion, .(ModelID = sample_ID, gene)], deletion = listed[is_deletion, .(ModelID = sample_ID, gene)],
  amplification = sheet[Genomic_alterations == dp$organoid_amplification_event & gene %in% ctx_par$amp_genes & CN / ploidy >= dp$organoid_amplification_ratio, .(ModelID = sample_ID, gene)]
))
# The catalogue covers the genes it lists for organoids; a gene it never lists has no status.
catalogue <- sort(unique(sheet$gene))
organoid_genes <- sort(unique(organoid_events$gene))
organoid_profiled <- organoids[!is.na(alt_TP53), ModelID]
organoid_status <- gene_status(organoids$ModelID, organoid_genes, organoid_events, organoid_profiled)
# The five genes of the stage-25 panel must come out as stage 25 stored them.
panel_same <- vapply(op$genotype_genes, function(symbol) {
  mutation_only <- drop(gene_status(organoids$ModelID, symbol, organoid_events[mutation == TRUE | deletion == TRUE], organoid_profiled))
  identical(unname(mutation_only), as.numeric(organoids[[paste0("alt_", symbol)]]))
}, logical(1))
if (!all(panel_same)) stop("organoid status differs from the stage-25 panel for: ", paste(op$genotype_genes[!panel_same], collapse = ", "))
write_ext(cbind(source = "Sanger organoids", organoid_events), "driver_type_organoid_events")
write_ext(as.data.table(organoid_status, keep.rownames = "ModelID"), "driver_type_organoid_status")

organoid_census <- split_census(organoid_status, organoids$OncotreeLineage, organoids$exposed, gate)
organoid_census[, `:=`(level = "tumour type", source = "Sanger organoids")]
organoid_splits <- organoid_census[label %in% split_labels[1:2]]
organoid_splits[, split_id := split_id_of("organoid", type, gene)]
organoid_facts <- function(i) {
  sframe <- split_frame(organoids, organoids$OncotreeLineage %in% organoid_splits$type[i], organoid_status[, organoid_splits$gene[i]])
  designs <- split_designs(sframe, dp$organoid_covariates, gate, lineage = FALSE, min_residual_df = sp$min_residual_df)
  cells <- designs$cells
  data.table(
    lines = nrow(sframe), residual_df = nrow(designs$inside$design) - ncol(designs$inside$design),
    kept_inside = information_kept(designs$inside$design, "exposed_in", cells[["neg_in"]], cells[["pos_in"]]),
    inside_unit = sqrt(solve(crossprod(designs$inside$design))["exposed_in", "exposed_in"]),
    specificity_fitted = designs$specificity && !is.null(designs$interaction)
  )
}
organoid_splits <- cbind(organoid_splits, rbindlist(lapply(seq_len(nrow(organoid_splits)), organoid_facts)))
organoid_splits[, balanced_equivalent := 2 / inside_unit^2]
organoid_splits[, low_information := kept_inside < dp$information_flag]
organoid_splits[, exact_floor := exact_floor(neg_alt, pos_alt)]
write_ext(organoid_splits, "driver_type_organoid_splits")
organoid_lines <- rbindlist(lapply(seq_len(nrow(organoid_splits)), function(i) {
  sframe <- split_frame(organoids, organoids$OncotreeLineage %in% organoid_splits$type[i], organoid_status[, organoid_splits$gene[i]])
  sframe[, .(split_id = organoid_splits$split_id[i], ModelID, exposed, ctx)]
}))
write_ext(organoid_lines, "driver_type_organoid_split_lines")
print(organoid_splits[, .(type, gene, neg_alt, pos_alt, neg_wt, pos_wt, label, kept_inside = round(kept_inside, 2), balanced = round(balanced_equivalent, 1))])

# --- DepMap NextGen: census only -------------------------------------------------------------------------
nextgen <- read_ext("nextgen_cohort", dir = ext_dir)[in_primary == TRUE]
nextgen_status <- sapply(dp$nextgen_statuses, function(column) as.numeric(nextgen[[column]]))
rownames(nextgen_status) <- nextgen$ModelID
nextgen_census <- split_census(nextgen_status, nextgen$OncotreeLineage, as.numeric(nextgen$y_negative == 1), gate)
nextgen_census[, `:=`(level = "lineage", source = "NextGen")]

census <- rbind(census_cells, organoid_census, nextgen_census)
setcolorder(census, c("source", "level", "type", "gene"))
write_ext(census, "driver_type_census")
census_summary <- census[, .(
  pairs_with_an_altered_line = .N, genes = uniqueN(gene), cancer_types = uniqueN(type),
  specificity_testable = sum(label == split_labels[1]), altered_only = sum(label == split_labels[2]),
  below_the_gate = sum(label == split_labels[3]), one_y_group_missing = sum(label == split_labels[4]),
  y_effect_testable = sum(label %in% split_labels[1:2])
), by = .(source, level)]
write_ext(census_summary, "driver_type_census_summary")
print(census_summary)

# --- pairs for the combination, with the cancer-type pairs of stage 28 --------------------------------
combine <- as.data.table(ext_params$combine$pairs)
pairs <- rbindlist(lapply(seq_len(nrow(combine)), function(i) {
  context_row <- definitions[context_id == combine$cell_line_context[i]]
  stopifnot(nrow(context_row) == 1)
  cell <- splits[level == context_row$class & type == context_row$label]
  organoid <- organoid_splits[type == combine$organoid_type[i]]
  both <- intersect(cell$gene, organoid$gene)
  if (length(both) == 0) return(NULL)
  cell <- cell[match(both, gene)]
  data.table(
    pair = combine$label[i], gene = both, cell_line_level = context_row$class, cell_line_type = context_row$label,
    cell_split_id = fifelse(is.na(cell$same_lines_as), cell$split_id, cell$same_lines_as), organoid_split_id = organoid[match(both, gene), split_id],
    cell_line_msi_covariate = combine$cell_line_msi_covariate[i]
  )
}))
write_ext(pairs, "driver_type_pairs")
print(pairs)

# --- gene families: those of stages 12 and 26 ---------------------------------------------------------
gene_ann <- read_result("gene_annotation")
pan_hits <- read_result("gene_pan_hits")
focus <- union(pan_hits$feature, gene_ann[in_scan == TRUE & !is.na(y_partner), feature])
cell_universe <- union(gene_ann[selective == TRUE, feature], focus)
organoid_ann <- read_ext("organoid_gene_annotation", dir = ext_dir)
organoid_selective <- organoid_ann[selective == TRUE, symbol]
organoid_universe <- union(organoid_selective, union(organoid_ann[cell_line_selective == TRUE & in_scan == TRUE, symbol], sp$x_partners))
organoid_universe <- intersect(organoid_universe, organoid_ann[in_scan == TRUE, symbol])
gene_families <- rbind(
  data.table(source = "Cell lines", feature = cell_universe, symbol = gene_symbol(cell_universe), selective = cell_universe %in% gene_ann[selective == TRUE, feature], focused = cell_universe %in% focus),
  data.table(source = "Sanger organoids", feature = organoid_universe, symbol = organoid_universe, selective = organoid_universe %in% organoid_selective, focused = organoid_universe %in% sp$x_partners)
)
write_ext(gene_families, "driver_type_genes")

# --- counts expected from the exploration; a different input needs a new decision -------------------------
expected <- dp$expected_counts
count_at <- function(level_name, labels) census_cells[level == level_name & label %in% labels, .N]
found <- list(
  lines = nrow(frame), y_negative = sum(frame$exposed == 1), genes_altered = length(genes),
  pairs_lineage = census_cells[level == "lineage", .N], pairs_disease = census_cells[level == "disease", .N], pairs_subtype = census_cells[level == "subtype", .N],
  inside_lineage = count_at("lineage", split_labels[1:2]), inside_disease = count_at("disease", split_labels[1:2]), inside_subtype = count_at("subtype", split_labels[1:2]),
  specificity_lineage = count_at("lineage", split_labels[1]), specificity_disease = count_at("disease", split_labels[1]), specificity_subtype = count_at("subtype", split_labels[1]),
  distinct_splits = sum(splits$fitted), distinct_specificity = splits[fitted == TRUE, sum(label == split_labels[1])], distinct_altered_sets = uniqueN(splits$altered_set),
  organoid_splits = nrow(organoid_splits), organoid_specificity = organoid_splits[, sum(label == split_labels[1])],
  nextgen_splits = nextgen_census[label %in% split_labels[1:2], .N], combined_pairs = nrow(pairs)
)
agree <- vapply(names(expected), function(key) within_tolerance(found[[key]], expected[[key]]), logical(1))
print(data.table(count = names(expected), expected = unlist(expected), found = unlist(found[names(expected)]), within_tolerance = agree))
stopifnot(found$lines == expected$lines, all(agree))

# --- freeze ---------------------------------------------------------------------------------------------
fitted_splits <- splits[fitted == TRUE]
summary_list <- list(
  lines = nrow(frame), y_negative = sum(frame$exposed == 1), genes_altered = length(genes),
  altered_pairs = nrow(events), by_kind = list(mutation = sum(events$mutation), amplification = sum(events$amplification), deletion = sum(events$deletion)),
  stored_contexts_reproduced = nrow(stored), counts_found = found, counts_expected = expected,
  settings = list(
    gate = gate, information_flag = dp$information_flag, shared_lines_fraction = dp$shared_lines_fraction, amplification_ratio = ctx_par$amp_ratio_min, deletion_ratio = ctx_par$del_ratio_max,
    amplification_genes = length(ctx_par$amp_genes), deletion_genes = length(ctx_par$del_genes), robust_permutation_p = sp$robust_permutation_p, robust_leave_one_out_p = sp$robust_leave_one_out_p,
    null_runs = dp$null_runs, followup_permutations = sp$followup_permutations
  ),
  splits = list(
    rows = nrow(splits), fitted = nrow(fitted_splits), specificity = fitted_splits[, sum(label == split_labels[1])], altered_sets = uniqueN(splits$altered_set),
    shared_groups = uniqueN(splits$shared_group), cancer_types = uniqueN(fitted_splits$type), genes = sort(unique(fitted_splits$gene)),
    smaller_altered_cell_below_10 = fitted_splits[, sum(pmin(neg_alt, pos_alt) < 10)],
    kept_inside_median = median(fitted_splits$kept_inside), kept_inside_min = min(fitted_splits$kept_inside), low_information = sum(fitted_splits$low_information),
    kept_alone_median = median(fitted_splits$kept_inside_alone, na.rm = TRUE), kept_alone_min = min(fitted_splits$kept_inside_alone, na.rm = TRUE), kept_alone_below_flag = sum(fitted_splits$kept_inside_alone < dp$information_flag, na.rm = TRUE),
    robust_not_possible = sum(!fitted_splits$robust_possible), same_donor_twice = sum(fitted_splits$same_donor_twice),
    below_five_donors = fitted_splits[, sum(pmin(donors_neg_alt, donors_pos_alt) < gate)]
  ),
  organoids = list(models = nrow(organoids), catalogue_genes = length(catalogue), genes_altered = length(organoid_genes), splits = nrow(organoid_splits), specificity = organoid_splits[, sum(label == split_labels[1])],
    curated_genes_outside_catalogue = setdiff(cn_genes, catalogue)),
  nextgen = list(models = nrow(nextgen), splits = found$nextgen_splits),
  pairs = nrow(pairs),
  gene_families = list(cell_lines = length(cell_universe), cell_lines_selective = sum(cell_universe %in% gene_ann[selective == TRUE, feature]), organoids = length(organoid_universe), organoids_selective = length(intersect(organoid_universe, organoid_selective))),
  parameters_sha256 = ext_section_hash(c("shared", "driver_type"))
)
jsonlite::write_json(summary_list, file.path(ext_out_dir, "driver_type_prepare_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
frozen <- function(name) file.path(ext_out_dir, name)
write_extension_freeze(
  "driver_type",
  paths = c(
    status = frozen("driver_type_status.csv"), organoid_status = frozen("driver_type_organoid_status.csv"), census = frozen("driver_type_census.csv"),
    splits = frozen("driver_type_splits.csv"), split_lines = frozen("driver_type_split_lines.csv"),
    organoid_splits = frozen("driver_type_organoid_splits.csv"), organoid_split_lines = frozen("driver_type_organoid_split_lines.csv"),
    pairs = frozen("driver_type_pairs.csv"), genes = frozen("driver_type_genes.csv")
  ),
  sections = c("shared", "driver_type"),
  counts = list(splits = nrow(fitted_splits), specificity = summary_list$splits$specificity, organoid_splits = nrow(organoid_splits), pairs = nrow(pairs), cell_line_genes = length(cell_universe)),
  replace = "--replace" %in% commandArgs(trailingOnly = TRUE)
)
log_msg("driver-by-type preparation done")
