# Helpers for the extension stages (25 to 36): organoids, the RNAi context
# screen, the combination across sources, BeatAML2, DepMap NextGen and driver status
# inside cancer types. Source after common.R,
# scan_helpers.R, context_helpers.R and power_helpers.R.
#
# Nothing here writes to the frozen analysis folder. Extension tables live in
# derived/<analysis>/extension/ and carry one of the prefixes below.
# LOY_EXTENSION_PARAMETERS and LOY_EXTENSION_OUTPUT point a trial run at another
# parameter file and output folder.

ext_dir <- file.path(derived_dir, "extension")
ext_parameter_file <- Sys.getenv("LOY_EXTENSION_PARAMETERS", unset = file.path(script_dir, "extension_parameters.json"))
ext_out_dir <- Sys.getenv("LOY_EXTENSION_OUTPUT", unset = ext_dir)
ext_params <- jsonlite::fromJSON(ext_parameter_file)
ext_prefixes <- c("organoid_", "rnai_", "combined_", "beataml_", "nextgen_", "driver_type_", "extension_")

robust_sd <- function(x) 1.4826 * median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)

# --- guarded output -----------------------------------------------------------------

write_ext <- function(dt, name, dir = ext_out_dir) {
  if (!any(startsWith(name, ext_prefixes))) {
    stop(sprintf("%s: an extension table must start with one of %s", name, paste(ext_prefixes, collapse = " ")))
  }
  if (identical(normalizePath(dir, mustWork = FALSE), normalizePath(derived_dir, mustWork = FALSE))) {
    stop("extension tables are never written to the frozen folder")
  }
  write_result(dt, name, dir = dir)
}

read_ext <- function(name, dir = ext_out_dir) read_result(name, dir = dir)

# Write a table unless the stored one already holds the same values. A refit with lm
# differs in its last digits between compute nodes (AVX-512 or not), so a table that
# only repeats an earlier result would otherwise change its bytes from run to run.
# Numbers are compared one by one, relative to their own size; anything else exactly.
# Returns TRUE when the table was written and FALSE when the stored one was kept.
write_result_stable <- function(dt, name, dir = derived_dir, tolerance = 1e-8) {
  paths <- file.path(dir, paste0(name, c(".csv", ".rds")))
  if (all(file.exists(paths))) {
    stored <- as.data.table(readRDS(paths[2]))
    same_column <- function(column) {
      a <- stored[[column]]
      b <- dt[[column]]
      if (!is.numeric(a) || !is.numeric(b)) return(identical(as.character(a), as.character(b)))
      identical(is.na(a), is.na(b)) && all(abs(a - b)[!is.na(a)] <= tolerance * pmax(abs(a), abs(b))[!is.na(a)])
    }
    if (identical(dim(stored), dim(dt)) && identical(names(stored), names(dt)) && all(vapply(names(dt), same_column, logical(1)))) {
      log_msg("kept %s as stored: the values agree to %g of their size", name, tolerance)
      return(invisible(FALSE))
    }
  }
  write_result(dt, name, dir = dir)
  invisible(TRUE)
}

# --- second freeze --------------------------------------------------------------------
# Each preparation stage freezes the parameter sections it depends on and the
# tables it prepared; the screen stages assert that freeze before any test.

ext_section_hash <- function(sections) {
  text <- as.character(jsonlite::toJSON(ext_params[sections], auto_unbox = TRUE, digits = NA, null = "null"))
  digest::digest(text, algo = "sha256", serialize = FALSE)
}

ext_freeze_path <- function(section, dir = ext_out_dir) file.path(dir, sprintf("extension_freeze_%s.json", section))

write_extension_freeze <- function(section, paths, sections, counts = list(), dir = ext_out_dir, replace = FALSE) {
  stopifnot(length(names(paths)) == length(paths), all(file.exists(paths)))
  current <- list(
    parameter_sections = sections, parameters_sha256 = ext_section_hash(sections),
    paths = as.list(paths), sha256 = as.list(vapply(paths, file_sha256, character(1)))
  )
  target <- ext_freeze_path(section, dir)
  if (file.exists(target)) {
    previous <- jsonlite::fromJSON(target)
    if (identical(previous$parameters_sha256, current$parameters_sha256) && identical(previous$sha256, current$sha256)) {
      log_msg("extension freeze '%s' already in place and unchanged", section)
      return(invisible(TRUE))
    }
    if (!replace) stop(sprintf("a different extension freeze '%s' already exists; rerun with --replace to archive it", section))
    archived <- sub("\\.json$", sprintf("_archived_%s.json", format(Sys.time(), "%Y%m%dT%H%M%S")), target)
    file.rename(target, archived)
    log_msg("previous extension freeze archived as %s", basename(archived))
  }
  jsonlite::write_json(c(
    list(created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), slurm_job_id = Sys.getenv("SLURM_JOB_ID"), section = section),
    current, list(counts = counts)
  ), target, auto_unbox = TRUE, pretty = TRUE)
  log_msg("extension freeze written: %s", target)
  invisible(TRUE)
}

assert_extension_freeze <- function(section, dir = ext_out_dir) {
  target <- ext_freeze_path(section, dir)
  if (!file.exists(target)) stop(sprintf("extension freeze '%s' is missing; run its preparation stage first", section))
  frozen <- jsonlite::fromJSON(target)
  if (!identical(ext_section_hash(frozen$parameter_sections), frozen$parameters_sha256)) {
    stop(sprintf("extension parameters (%s) changed after the '%s' freeze", paste(frozen$parameter_sections, collapse = ", "), section))
  }
  for (key in names(frozen$paths)) {
    if (!identical(file_sha256(frozen$paths[[key]]), frozen$sha256[[key]])) stop(sprintf("%s changed after the '%s' freeze", key, section))
  }
  log_msg("extension freeze verified: %s", section)
  invisible(TRUE)
}

# --- downloads --------------------------------------------------------------------------

verify_download <- function(path, checksum, algorithm) {
  if (!file.exists(path)) stop(sprintf("%s is missing", path))
  found <- digest::digest(file = path, algo = algorithm)
  if (!identical(tolower(found), tolower(checksum))) {
    stop(sprintf("%s: %s checksum %s does not match the recorded %s", path, algorithm, found, checksum))
  }
  data.table(file = path, bytes = file.size(path), algorithm = algorithm, checksum = found)
}

# --- Y status on another dataset's scale ---------------------------------------------------

# The RNA rule of 03_call_y_status.R for any expression table (samples in rows,
# the core Y genes in columns, log2(TPM + 1)). With `thresholds` it applies the
# given cut-offs. Without, it derives them from the dataset's own reference
# samples: a gene is "on" above its female maximum, "depleted" means a mean inside
# the female range, and "expressed" means a mean of at least `expressed_fraction`
# of the typical male level in linear TPM.
# `ceiling_fraction` is for a table whose females show no Y signal at all, where
# the female range would leave out males with a trace of signal. A gene is then
# "on" above its female maximum or that share of its typical male level,
# whichever is larger, and the mean ceiling is raised the same way. The typical
# male is taken over the males with enough genes on, so the ceilings and that
# set of males are settled together.
call_rna_y <- function(log_expr, female_ids = NULL, male_ids = NULL, expressed_fraction = NULL,
                       depleted_max_on = 1, expressed_min_on = 6, thresholds = NULL, ceiling_fraction = NULL) {
  mean_y <- rowMeans(log_expr)
  if (is.null(thresholds)) {
    stopifnot(all(female_ids %in% rownames(log_expr)), all(male_ids %in% rownames(log_expr)), !is.null(expressed_fraction))
    # One male among the reference females would lift every threshold.
    midpoint <- (median(mean_y[female_ids]) + median(mean_y[male_ids])) / 2
    odd <- female_ids[mean_y[female_ids] > midpoint]
    if (length(odd) > 0) stop("reference females with male-level Y expression: ", paste(odd, collapse = ", "))
    gene_max <- apply(log_expr[female_ids, , drop = FALSE], 2, max)
    on_reference <- rowSums(sweep(log_expr[male_ids, , drop = FALSE], 2, gene_max, ">"))
    male_reference <- median(mean_y[male_ids][on_reference >= expressed_min_on])
    thresholds <- list(
      gene_max = gene_max, depleted_mean_max = max(mean_y[female_ids]),
      expressed_mean_min = log2(expressed_fraction * (2^male_reference - 1) + 1), male_reference = male_reference
    )
    if (!is.null(ceiling_fraction)) {
      as_share <- function(level) log2(ceiling_fraction * (2^level - 1) + 1)
      for (pass in seq_len(20)) {
        reference <- on_reference >= expressed_min_on
        raised <- pmax(gene_max, as_share(apply(log_expr[male_ids[reference], , drop = FALSE], 2, median)))
        on_reference <- rowSums(sweep(log_expr[male_ids, , drop = FALSE], 2, raised, ">"))
        if (identical(on_reference >= expressed_min_on, reference)) break
        if (pass == 20) stop("the male reference did not settle")
      }
      male_reference <- median(mean_y[male_ids][reference])
      thresholds <- list(
        gene_max = raised, depleted_mean_max = max(thresholds$depleted_mean_max, as_share(male_reference)),
        expressed_mean_min = log2(expressed_fraction * (2^male_reference - 1) + 1), male_reference = male_reference,
        female_gene_max = gene_max, female_mean_max = thresholds$depleted_mean_max, ceiling_fraction = ceiling_fraction, passes = pass
      )
    }
  }
  n_on <- rowSums(sweep(log_expr, 2, thresholds$gene_max[colnames(log_expr)], ">"))
  depleted <- mean_y <= thresholds$depleted_mean_max & n_on <= depleted_max_on
  expressed <- mean_y >= thresholds$expressed_mean_min & n_on >= expressed_min_on
  call <- ifelse(depleted, "depleted", ifelse(expressed, "expressed", "intermediate"))
  list(
    call = setNames(call, rownames(log_expr)), n_on = setNames(as.integer(n_on), rownames(log_expr)),
    mean_y = mean_y, thresholds = thresholds
  )
}

# RNA leads; DNA fills in where RNA is missing and vetoes a clear contradiction.
# The same rule as integrate_calls() in 03_call_y_status.R.
integrate_y_calls <- function(rna, dna, male) {
  rna_as <- fifelse(rna == "depleted", "Y_negative", fifelse(rna == "expressed", "Y_positive", NA_character_))
  dna_as <- fifelse(dna == "loss", "Y_negative", fifelse(dna == "retained", "Y_positive", NA_character_))
  rna_definite <- !is.na(rna_as)
  dna_definite <- !is.na(dna_as)
  dna_ambiguous <- !is.na(dna) & dna %in% c("partial", "mosaic")
  tier <- rep(NA_character_, length(rna))
  tier[rna_definite & dna_definite & rna_as == dna_as] <- "concordant"
  tier[rna_definite & dna_definite & rna_as != dna_as] <- "discordant"
  tier[rna_definite & dna_ambiguous] <- "rna_led_dna_ambiguous"
  tier[rna_definite & is.na(dna)] <- "rna_only"
  tier[is.na(rna) & dna_definite] <- "dna_only"
  tier[!is.na(rna) & rna == "intermediate"] <- "rna_intermediate"
  tier[is.na(rna) & dna_ambiguous] <- "dna_ambiguous_no_rna"
  tier[is.na(rna) & is.na(dna)] <- "no_data"
  call <- fifelse(
    tier %in% c("concordant", "rna_led_dna_ambiguous", "rna_only"), rna_as,
    fifelse(tier == "dna_only", dna_as, NA_character_)
  )
  tier[!male] <- "not_male"
  call[!male] <- NA_character_
  list(call = call, tier = tier)
}

# chrY copy-number segments (columns sample, start, end, ratio; ratio is copy
# number relative to the sample's ploidy) summarised over a region: the share of
# the region covered, the share of the covered part at the loss floor and at
# retained level, and the length-weighted mean ratio.
dna_segment_summary <- function(segments, start, end, loss_max, retained_min) {
  region <- c(start, end)
  seg <- as.data.table(segments)[, .(sample, start, end, ratio)]
  seg[, width := pmax(0, pmin(end, region[2]) - pmax(start, region[1]))]
  seg <- seg[width > 0]
  seg[, .(
    covered = sum(width) / diff(region), floor = sum(width[ratio <= loss_max]) / sum(width),
    retained = sum(width[ratio >= retained_min]) / sum(width), ratio = sum(width * ratio) / sum(width)
  ), by = sample]
}

call_dna_segments <- function(summary, floor_min, covered_min = 0.5) {
  out <- ifelse(summary$floor >= floor_min, "loss",
    ifelse(summary$retained >= floor_min, "retained", ifelse(summary$floor > 0 & summary$retained > 0, "partial", "mosaic")))
  out[summary$covered < covered_min] <- NA_character_
  out
}

# The DNA call from a gene-level table (samples in rows, the core Y genes in
# columns, copy number relative to ploidy): loss or retained when at least
# `min_genes` genes say so, "partial" when some genes are lost and others kept,
# "mosaic" when the level sits in between, and no call with fewer than
# `min_genes` values. `ratio` is the mean over the genes with a value.
dna_gene_call <- function(values, loss_max, retained_min, min_genes) {
  values <- as.matrix(values)
  n_genes <- rowSums(!is.na(values))
  n_loss <- rowSums(values <= loss_max, na.rm = TRUE)
  n_retained <- rowSums(values >= retained_min, na.rm = TRUE)
  call <- ifelse(n_loss >= min_genes, "loss", ifelse(n_retained >= min_genes, "retained", ifelse(n_loss > 0 & n_retained > 0, "partial", "mosaic")))
  call[n_genes < min_genes] <- NA_character_
  ratio <- rowMeans(values, na.rm = TRUE)
  ratio[n_genes == 0] <- NA_real_
  data.table(
    sample = rownames(values), n_genes = as.integer(n_genes), n_loss = as.integer(n_loss), n_retained = as.integer(n_retained),
    ratio = unname(ratio), call = unname(call)
  )
}

# --- donors -------------------------------------------------------------------------------

# The donor part of an HCMI model name ("HCM-SANG-0299-C15-B" gives
# "HCM-SANG-0299"); missing for any other name.
hcmi_donor <- function(name) {
  pattern <- "^(HCM-[A-Z]+-[0-9]+)-.*$"
  out <- sub(pattern, "\\1", name)
  out[is.na(name) | !grepl(pattern, name)] <- NA_character_
  out
}

# Jaccard index between models (rows of `mutated`, a 0/1 or count matrix of
# calls) on the calls carried by at most `max_fraction` of the models. Two models
# of one donor share most of their rare calls; unrelated models share few. A
# model without a rare call is similar to nobody.
genotype_similarity <- function(mutated, max_fraction) {
  carried <- as.matrix(mutated) > 0
  stopifnot(!anyNA(carried), !is.null(rownames(carried)))
  frequency <- colMeans(carried)
  rare <- carried[, frequency > 0 & frequency <= max_fraction, drop = FALSE] * 1
  shared <- tcrossprod(rare)
  either <- outer(rowSums(rare), rowSums(rare), "+") - shared
  out <- ifelse(either > 0, shared / either, 0)
  diag(out) <- 1
  dimnames(out) <- list(rownames(carried), rownames(carried))
  out
}

# One donor label per model. Two models are joined when they share the value of
# any key (a missing or empty value never matches) or when their similarity
# reaches `cutoff`; joins chain. `similarity` may cover only some of the models.
# The label is the first ID of the group in C order, so it does not depend on the
# order of the rows.
donor_groups <- function(ids, keys = list(), similarity = NULL, cutoff = NULL) {
  stopifnot(!anyNA(ids), !anyDuplicated(ids))
  parent <- seq_along(ids)
  find <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  join <- function(i, j) {
    a <- find(i)
    b <- find(j)
    if (a != b) parent[max(a, b)] <<- min(a, b)
  }
  for (key in keys) {
    stopifnot(length(key) == length(ids))
    known <- !is.na(key) & key != ""
    for (members in split(seq_along(ids)[known], key[known])) for (j in members[-1]) join(members[1], j)
  }
  if (!is.null(similarity)) {
    stopifnot(!is.null(cutoff), identical(rownames(similarity), colnames(similarity)))
    at <- match(rownames(similarity), ids)
    pairs <- which(similarity >= cutoff & upper.tri(similarity), arr.ind = TRUE)
    for (k in seq_len(nrow(pairs))) if (!anyNA(at[pairs[k, ]])) join(at[pairs[k, 1]], at[pairs[k, 2]])
  }
  root <- vapply(seq_along(ids), find, integer(1))
  label <- vapply(split(ids, root), function(x) sort(x, method = "radix")[1], character(1))
  setNames(unname(label[as.character(root)]), ids)
}

# For each model, the candidate that stands for the same donor in another cohort: the
# candidate of the same name if there is one, otherwise the donor's preferred candidate,
# otherwise its first by name. A model without a donor ID has no partner, and a candidate
# without one is nobody's partner.
partner_model <- function(name, donor, candidates, candidate_donor, preferred = rep(FALSE, length(candidates))) {
  stopifnot(length(name) == length(donor), length(candidates) == length(candidate_donor), length(candidates) == length(preferred))
  vapply(seq_along(name), function(i) {
    if (is.na(donor[i])) return(NA_character_)
    own <- which(!is.na(candidate_donor) & candidate_donor == donor[i])
    if (length(own) == 0) return(NA_character_)
    if (name[i] %in% candidates[own]) return(name[i])
    candidates[own[order(!preferred[own], candidates[own], method = "radix")][1]]
  }, character(1))
}

# TRUE for the one model kept per group: the best `quality` (smaller is
# better, as for NNMD), ties broken by ID.
best_of_group <- function(ids, group, quality) {
  stopifnot(!anyNA(ids), !anyNA(group), !anyNA(quality), !anyDuplicated(ids))
  ordered <- order(group, quality, ids, method = "radix")
  keep <- logical(length(ids))
  keep[ordered[!duplicated(group[ordered])]] <- TRUE
  keep
}

# Fraction of a mixed sample that has lost Y. `level` holds log2 values per
# sample (rows) and Y gene (columns), either the gene itself or the gene minus
# its X partner. Each gene is taken relative to the median male sample, the
# sample's ratio is the median over genes, and the fraction is one minus that
# ratio, rescaled so that the female level (if given) is 1. Values inside
# `floor_sd` robust SDs of the male bulk are set to 0 in `f`.
y_fraction <- function(level, male_ids, female_ids = NULL, floor_sd = 2) {
  stopifnot(all(male_ids %in% rownames(level)), all(female_ids %in% rownames(level)))
  reference <- apply(level[male_ids, , drop = FALSE], 2, median)
  ratio <- apply(2^sweep(level, 2, reference, "-"), 1, median)
  zero <- if (is.null(female_ids)) 0 else median(ratio[female_ids])
  f_raw <- pmin(1, pmax(0, (1 - ratio) / (1 - zero)))
  floor_value <- 1 - 2^(-floor_sd * robust_sd(log2(ratio[male_ids])))
  data.table(
    sample = rownames(level), ratio = unname(ratio), f_raw = unname(f_raw),
    f = unname(ifelse(f_raw < floor_value, 0, f_raw)), floor = floor_value
  )
}

# Fraction of cells in clones without a Y, from an ISCN karyotype string.
# A clone is read as "count,sex chromosomes,aberrations[cells]". "idem" and "sl"
# inherit the first clone, "sdlN" the Nth sideline, and then apply their own +Y
# and -Y. A clone written with a Y in its sex chromosomes is Y-positive whatever
# else it lists (a tetraploid XXY,-Y still has one). Two or more X, no Y and no
# -Y is a female clone (donor cells), flagged and not counted as Y loss.
parse_karyotype_y <- function(text) {
  count_of <- function(pattern, tokens) sum(tokens == pattern)
  one <- function(raw) {
    out <- list(
      parsed = FALSE, cells = NA_integer_, cells_without_y = NA_integer_, fraction = NA_real_,
      mentions_y_loss = NA, female_clone = NA, chimerism = NA, counts_missing = NA
    )
    if (is.na(raw)) return(out)
    # Spaces (also no-break ones) and FISH results appended after the karyotype are dropped.
    string <- sub("\\.(nuc)?ish.*$", "", gsub("[[:space:]\u00a0]+", "", raw))
    if (!grepl("^[0-9]", string)) return(out)
    out$mentions_y_loss <- grepl("-Y", string, fixed = TRUE)
    out$chimerism <- grepl("//", string, fixed = TRUE)
    clones <- strsplit(gsub("//", "/", string, fixed = TRUE), "/", fixed = TRUE)[[1]]
    clones <- clones[nzchar(clones)]
    # The cell count is the last "[n]" or "[cpn]" of a clone; text after it is a note.
    has_count <- grepl("\\[(cp)?[0-9]+\\]", clones)
    cells <- rep(NA_integer_, length(clones))
    cells[has_count] <- as.integer(sub("^.*\\[(cp)?([0-9]+)\\][^\\[]*$", "\\2", clones[has_count]))
    y <- integer(length(clones))
    female <- logical(length(clones))
    for (i in seq_along(clones)) {
      tokens <- strsplit(sub("\\[.*$", "", clones[i]), ",", fixed = TRUE)[[1]]
      if (length(tokens) < 2 || !grepl("^[0-9]", tokens[1])) return(out)
      sex <- tokens[2]
      if (grepl("^[XY]+$", sex)) {
        y[i] <- nchar(gsub("X", "", sex, fixed = TRUE))
        female[i] <- y[i] == 0 && nchar(sex) >= 2 && count_of("-Y", tokens) == 0
      } else if (grepl("^(idem|sl|sdl[0-9]*)(x[0-9]+)?$", sex)) {
        # "x2" marks a doubled clone, which carries its parent's sex chromosomes twice.
        copies <- if (grepl("x[0-9]+$", sex)) as.integer(sub("^.*x", "", sex)) else 1L
        reference <- sub("x[0-9]+$", "", sex)
        parent <- if (grepl("^sdl", reference)) 1L + max(1L, as.integer(sub("sdl", "", reference)), na.rm = TRUE) else 1L
        if (parent >= i) return(out)
        y[i] <- y[parent] * copies + count_of("+Y", tokens) - count_of("-Y", tokens)
        female[i] <- female[parent]
      } else {
        return(out)
      }
    }
    lost <- y <= 0 & !female
    out$parsed <- TRUE
    out$female_clone <- any(female)
    out$counts_missing <- any(!has_count)
    if (!any(!has_count)) {
      out$cells <- sum(cells)
      out$cells_without_y <- sum(cells[lost])
      out$fraction <- out$cells_without_y / out$cells
    } else if (length(clones) == 1) {
      out$fraction <- as.numeric(lost)
    }
    out
  }
  rbindlist(lapply(text, function(raw) c(list(text = raw), one(raw))))
}

# --- context fits for genes with different line sets ----------------------------------------

# Genes (rows) grouped by the set of lines (columns) that have a value.
pattern_groups <- function(mat) {
  present <- !is.na(mat)
  key <- apply(present, 1, function(v) paste(which(!v), collapse = ","))
  lapply(unname(split(rownames(mat), key)), function(features) {
    list(features = features, lines = colnames(mat)[present[features[1], ]])
  })
}

# limma's moderated t (robust empirical Bayes, no trend) from per-gene
# estimates, unscaled standard deviations, residual SDs and degrees of freedom.
moderated_test <- function(estimate, unscaled, sigma, df_residual) {
  squeezed <- squeezeVar(sigma^2, df_residual, robust = TRUE)
  df_total <- pmin(df_residual + squeezed$df.prior, sum(df_residual, na.rm = TRUE))
  se <- unscaled * sqrt(squeezed$var.post)
  t <- estimate / se
  data.table(t = t, p = 2 * pt(-abs(t), df_total), df = df_total, se = se)
}

# The two fits of the context screen (interaction on every line, Y effect inside
# the context) for rows of `mat` that may each have their own set of lines.
# Genes sharing a line set are fitted together on those lines; a gene is tested
# only if every Y-by-context cell holds `min_cell` of its own lines; the variance
# prior is pooled once over all tested genes, separately for the two fits.
fit_context_patterns <- function(ctx_frame, mat, covariates, min_cell, min_residual_df = 3, weights = NULL, lineage = TRUE) {
  mat <- mat[, ctx_frame$ModelID, drop = FALSE]
  blank <- function(features) data.table(estimate = rep(NA_real_, length(features)), unit = NA_real_, sigma = NA_real_, df_residual = NA_real_)
  fit_one <- function(sub_frame, features, tested, fit_covariates) {
    design <- build_design(sub_frame, tested = tested, covariates = fit_covariates, lineage = lineage)
    if (!design$estimable || nrow(design$design) - ncol(design$design) < min_residual_df) return(blank(features))
    rows <- sub_frame$ModelID[design$rows]
    y <- mat[features, rows, drop = FALSE]
    weight_matrix <- if (is.null(weights)) NULL else matrix(weights[rows], nrow(y), ncol(y), byrow = TRUE)
    fit <- lmFit(y, design$design, weights = weight_matrix)
    data.table(estimate = fit$coefficients[, tested], unit = fit$stdev.unscaled[, tested], sigma = fit$sigma, df_residual = fit$df.residual)
  }
  parts <- rbindlist(lapply(pattern_groups(mat), function(group) {
    sub_frame <- ctx_frame[ModelID %in% group$lines]
    cells <- context_cells(sub_frame$exposed, sub_frame$ctx, "binary")
    out <- data.table(
      feature = group$features, n = nrow(sub_frame), n_neg_in = cells[["neg_in"]], n_pos_in = cells[["pos_in"]],
      n_neg_out = cells[["neg_out"]], n_pos_out = cells[["pos_out"]]
    )
    supported <- all(cells >= min_cell)
    interaction <- if (supported) fit_one(sub_frame, group$features, "interaction", c("exposed", "ctx", covariates)) else blank(group$features)
    inside <- if (supported) fit_one(sub_frame[ctx == 1], group$features, "exposed", covariates) else blank(group$features)
    setnames(interaction, paste0("interaction_", names(interaction)))
    setnames(inside, paste0("inside_", names(inside)))
    cbind(out, interaction, inside)
  }))
  for (test in c("interaction", "inside")) {
    column <- function(name) paste0(test, "_", name)
    ok <- is.finite(parts[[column("sigma")]])
    for (name in c("p", "se", "df")) parts[, (column(name)) := NA_real_]
    if (any(ok)) {
      moderated <- moderated_test(parts[[column("estimate")]][ok], parts[[column("unit")]][ok], parts[[column("sigma")]][ok], parts[[column("df_residual")]][ok])
      parts[ok, (column("p")) := moderated$p]
      parts[ok, (column("se")) := moderated$se]
      parts[ok, (column("df")) := moderated$df]
    }
  }
  parts[match(rownames(mat), feature)]
}

# --- weighted scan ---------------------------------------------------------------------------

# limma scan of one coefficient with a weight per sample. Returns the columns
# of limma_scan() plus the design unit (standard error per unit residual SD)
# and the standardised effect.
weighted_limma_scan <- function(mat_t, design, coef, weights = NULL, trend = TRUE) {
  stopifnot(ncol(mat_t) == nrow(design), coef %in% colnames(design))
  if (anyNA(mat_t)) stop("weighted_limma_scan needs a matrix without missing values")
  weight_matrix <- if (is.null(weights)) NULL else matrix(weights, nrow(mat_t), ncol(mat_t), byrow = TRUE)
  fit <- eBayes(lmFit(mat_t, design, weights = weight_matrix), trend = trend, robust = TRUE)
  tab <- topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = TRUE)
  unit <- unname(fit$stdev.unscaled[, coef])
  data.table(
    feature = rownames(mat_t), estimate = tab$logFC, ci_low = tab$CI.L, ci_high = tab$CI.R,
    t = tab$t, p = tab$P.Value, df = fit$df.total, n = ncol(mat_t),
    unit = unit, sigma_post = sqrt(fit$s2.post), smd = tab$t * unit
  )
}

# --- combination across sources -----------------------------------------------------------------

signed_z <- function(p, estimate) sign(estimate) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)

# Weighted combination of signed z-scores (rows = features, columns = sources).
# With weights 1 / design unit this is fixed-effect pooling of standardised
# effects. A missing source drops out of that feature's combination.
stouffer <- function(z, w) {
  z <- as.matrix(z)
  w <- if (is.matrix(w)) w else matrix(w, nrow(z), ncol(z), byrow = TRUE)
  w[is.na(z)] <- 0
  combined <- rowSums(w * replace(z, is.na(z), 0)) / sqrt(rowSums(w^2))
  list(z = combined, p = 2 * pnorm(-abs(combined)), n_sources = as.integer(rowSums(!is.na(z))))
}

# Difference between two sources' standardised effects, with its interval.
source_difference <- function(smd, unit) {
  smd <- matrix(smd, ncol = 2)
  unit <- if (is.matrix(unit)) unit else matrix(unit, nrow(smd), 2, byrow = TRUE)
  difference <- smd[, 1] - smd[, 2]
  se <- sqrt(rowSums(unit^2))
  data.table(
    difference = difference, se = se, ci_low = difference - qnorm(0.975) * se, ci_high = difference + qnorm(0.975) * se,
    p = 2 * pnorm(-abs(difference / se))
  )
}

# How a combined hit is carried: every source other than the strongest one must
# pass a one-sided test in the direction of the combined effect (BH across the
# hits) for "supported by both"; a source significantly on the other side makes
# it "sources disagree".
evidence_class <- function(z, w, alpha = 0.05) {
  z <- as.matrix(z)
  w_matrix <- if (is.matrix(w)) w else matrix(w, nrow(z), ncol(z), byrow = TRUE)
  direction <- sign(stouffer(z, w_matrix)$z)
  strongest <- max.col(abs(w_matrix * replace(z, is.na(z), 0)), ties.method = "first")
  other <- !is.na(z)
  other[cbind(seq_len(nrow(z)), strongest)] <- FALSE
  toward <- pnorm(direction * z, lower.tail = FALSE)
  toward[!other] <- NA
  adjusted <- matrix(bh(as.vector(toward)), nrow(z))
  against <- pnorm(-direction * z, lower.tail = FALSE)
  out <- rep("carried by one source", nrow(z))
  out[rowSums(other) > 0 & rowSums(other & adjusted >= alpha, na.rm = TRUE) == 0] <- "supported by both"
  out[rowSums(other & against < alpha, na.rm = TRUE) > 0] <- "sources disagree"
  out
}

# evidence_class for any number of sources: how many sources a feature has, how
# many support the combined direction (the strongest one, plus every other whose
# one-sided p-value in that direction passes BH over all the entries passed in),
# and a label. With two sources present the label means what evidence_class says.
evidence_summary <- function(z, w, alpha = 0.05) {
  z <- as.matrix(z)
  w_matrix <- if (is.matrix(w)) w else matrix(w, nrow(z), ncol(z), byrow = TRUE)
  direction <- sign(stouffer(z, w_matrix)$z)
  strongest <- max.col(abs(w_matrix * replace(z, is.na(z), 0)), ties.method = "first")
  other <- !is.na(z)
  other[cbind(seq_len(nrow(z)), strongest)] <- FALSE
  toward <- pnorm(direction * z, lower.tail = FALSE)
  toward[!other] <- NA
  adjusted <- matrix(bh(as.vector(toward)), nrow(z))
  against <- pnorm(-direction * z, lower.tail = FALSE)
  n_sources <- as.integer(rowSums(!is.na(z)))
  n_supporting <- as.integer(1 + rowSums(other & adjusted < alpha, na.rm = TRUE))
  label <- ifelse(n_sources == 1, "one source only",
    ifelse(n_supporting == n_sources, "supported by every source",
      ifelse(n_supporting == 1, "carried by one source", sprintf("supported by %d of %d sources", n_supporting, n_sources))))
  label[rowSums(other & against < alpha, na.rm = TRUE) > 0] <- "sources disagree"
  data.table(n_sources = n_sources, n_supporting = n_supporting, label = label)
}

# The bar a p-value must stay below for a family-wise error of `level`, from
# the smallest p-value of each null run: the largest value that fewer than that
# share of the null runs fall below.
family_wise_bar <- function(min_p, level = 0.05) {
  stopifnot(!anyNA(min_p), length(min_p) >= 1 / level)
  sort(min_p)[max(1, floor(level * length(min_p) + 1e-9))]
}

# What a second source says about an effect found in a first one. Replicated:
# the same sign and an adjusted p-value below alpha. Not replicated: not that,
# and significantly weaker than the reference (`difference_p` is the two-sided
# p-value of estimate against reference). Otherwise the source cannot tell.
replication_outcome <- function(estimate, adjusted_p, reference, difference_p, alpha = 0.05) {
  replicated <- sign(estimate) == sign(reference) & adjusted_p < alpha
  weaker <- sign(reference) * (reference - estimate) > 0
  ifelse(replicated, "replicated", ifelse(weaker & difference_p < alpha, "not replicated", "inconclusive"))
}

# --- pooled test across sources -----------------------------------------------------------------

# rank_normalise applied to each row of `mat` (features in rows, samples in
# columns) within each source's samples.
rank_within_source <- function(mat, source) {
  stopifnot(length(source) == ncol(mat), !anyNA(source))
  out <- mat
  for (one in unique(source)) {
    columns <- which(source == one)
    out[, columns] <- t(apply(mat[, columns, drop = FALSE], 1, rank_normalise))
  }
  out
}

# Cancer types that no source can test alone but the sources together can.
# `counts` has one row per type and source with n_negative and n_positive.
pooled_gate <- function(counts, min_group) {
  by_type <- as.data.table(counts)[, .(
    alone = any(n_negative >= min_group & n_positive >= min_group),
    together = sum(n_negative) >= min_group & sum(n_positive) >= min_group
  ), by = type]
  by_type[!alone & together, sort(type)]
}

# Difference between two groups within strata, with a permutation p-value.
# `mat` has features in rows and samples in columns, `exposed` is 0/1. The
# estimate is the coefficient of `exposed` in a regression with one term per
# stratum; a stratum holding one group adds nothing. The reference set is every
# arrangement of `exposed` within strata, enumerated when there are at most
# `max_arrangements` and sampled (`n_random` draws plus the observed one)
# otherwise. p is two-sided, twice the smaller tail, and cannot fall below
# `floor`. p_family is the share of arrangements whose largest standardised
# statistic over all features is at least the feature's own.
stratified_permutation_test <- function(mat, exposed, strata, max_arrangements = 2e6, n_random = 1e5, seed = 1, cells = 2e7) {
  mat <- as.matrix(mat)
  stopifnot(length(exposed) == ncol(mat), length(strata) == ncol(mat), all(exposed %in% c(0, 1)), !anyNA(strata))
  if (anyNA(mat)) stop("stratified_permutation_test needs a matrix without missing values")
  members <- split(seq_along(exposed), as.character(strata))
  members <- members[vapply(members, function(i) sum(exposed[i]) > 0 && sum(exposed[i]) < length(i), logical(1))]
  if (length(members) == 0) stop("no stratum holds both groups")
  size <- lengths(members)
  n_exposed <- vapply(members, function(i) sum(exposed[i]), numeric(1))
  per_stratum <- choose(size, n_exposed)
  n_arrangements <- prod(per_stratum)
  denominator <- sum(n_exposed * (size - n_exposed) / size)
  # One column per stratum, also when there is one feature or one stratum.
  over_strata <- function(f) rowSums(matrix(vapply(seq_along(members), f, numeric(nrow(mat))), nrow = nrow(mat)))
  # Sum over the exposed samples, minus its mean over arrangements.
  centre <- over_strata(function(s) n_exposed[s] * rowMeans(mat[, members[[s]], drop = FALSE]))
  observed <- over_strata(function(s) rowSums(mat[, members[[s]][exposed[members[[s]]] == 1], drop = FALSE])) - centre
  spread <- sqrt(over_strata(function(s) {
    y <- mat[, members[[s]], drop = FALSE]
    n_exposed[s] * (size[s] - n_exposed[s]) / (size[s] * (size[s] - 1)) * rowSums((y - rowMeans(y))^2)
  }))
  scale <- ifelse(spread > 0, spread, Inf)
  slack <- 1e-9 * spread
  statistic <- observed / scale
  exact <- n_arrangements <= max_arrangements
  if (exact) {
    sums <- lapply(seq_along(members), function(s) {
      chosen <- combn(size[s], n_exposed[s])
      indicator <- matrix(0, size[s], ncol(chosen))
      indicator[cbind(as.vector(chosen), rep(seq_len(ncol(chosen)), each = nrow(chosen)))] <- 1
      mat[, members[[s]], drop = FALSE] %*% indicator
    })
    lower <- upper <- numeric(nrow(mat))
    largest <- numeric(n_arrangements)
    chunk <- max(1, floor(cells / n_arrangements))
    for (start in seq(1, nrow(mat), by = chunk)) {
      rows <- start:min(nrow(mat), start + chunk - 1)
      values <- sums[[1]][rows, , drop = FALSE]
      for (s in seq_along(sums)[-1]) {
        values <- values[, rep(seq_len(ncol(values)), times = ncol(sums[[s]])), drop = FALSE] +
          sums[[s]][rows, rep(seq_len(ncol(sums[[s]])), each = ncol(values)), drop = FALSE]
      }
      values <- values - centre[rows]
      lower[rows] <- rowMeans(values <= observed[rows] + slack[rows])
      upper[rows] <- rowMeans(values >= observed[rows] - slack[rows])
      largest <- pmax(largest, apply(abs(values) / scale[rows], 2, max))
    }
    p_family <- vapply(abs(statistic), function(z) mean(largest >= z - 1e-9), numeric(1))
    floor_p <- 2 / n_arrangements
  } else {
    set.seed(seed)
    lower <- upper <- numeric(nrow(mat))
    largest <- numeric(n_random)
    chunk <- max(1, floor(cells / nrow(mat)))
    for (start in seq(1, n_random, by = chunk)) {
      draws <- start:min(n_random, start + chunk - 1)
      arrangement <- matrix(0, ncol(mat), length(draws))
      for (s in seq_along(members)) {
        # The first n_exposed of a random order of the stratum: a uniform draw of that many samples.
        picked <- apply(matrix(runif(size[s] * length(draws)), size[s]), 2, order)[seq_len(n_exposed[s]), , drop = FALSE]
        arrangement[cbind(members[[s]][as.vector(picked)], rep(seq_along(draws), each = n_exposed[s]))] <- 1
      }
      values <- mat %*% arrangement - centre
      lower <- lower + rowSums(values <= observed + slack)
      upper <- upper + rowSums(values >= observed - slack)
      largest[draws] <- apply(abs(values) / scale, 2, max)
    }
    lower <- (1 + lower) / (n_random + 1)
    upper <- (1 + upper) / (n_random + 1)
    p_family <- vapply(abs(statistic), function(z) (1 + sum(largest >= z - 1e-9)) / (n_random + 1), numeric(1))
    floor_p <- 2 / (n_random + 1)
  }
  list(
    table = data.table(
      feature = rownames(mat), estimate = unname(observed / denominator), statistic = unname(statistic),
      p = unname(pmin(1, 2 * pmin(lower, upper))), p_family = unname(p_family)
    ),
    exact = exact, arrangements = n_arrangements, floor = floor_p, n_random = if (exact) NA_real_ else n_random,
    strata = data.table(stratum = names(members), n = as.integer(size), n_exposed = as.integer(n_exposed), arrangements = per_stratum)
  )
}

# --- least squares and permutation for the drug screen ----------------------------------------------

# fit_ols for every column of y_mat (rows = samples); missing responses are
# dropped per column.
scan_ols <- function(y_mat, x, z) {
  out <- rbindlist(lapply(seq_len(ncol(y_mat)), function(j) as.data.table(fit_ols(y_mat[, j], x, z))))
  out[, feature := colnames(y_mat)]
  setcolorder(out, "feature")
  out[]
}

permutation_matrix <- function(n, n_perm, seed) {
  set.seed(seed)
  vapply(seq_len(n_perm), function(b) sample.int(n), integer(n))
}

# The permutations of all samples restricted to the rows in `keep`: each column
# becomes a permutation of seq_along(keep) that keeps the full one's order, so
# responses with different missing samples share the same shuffles.
restrict_permutations <- function(perms, keep) {
  matrixStats::colRanks(perms[keep, , drop = FALSE], ties.method = "first", preserveShape = TRUE)
}

# Freedman-Lane permutation test of one coefficient (x) adjusted for z: the
# residuals of the model without x are shuffled, added back to its fit, and the
# least-squares t of x is recomputed. `perms` holds one permutation per column.
freedman_lane <- function(y, x, z, perms) {
  stopifnot(length(y) == nrow(perms), !anyNA(y), !anyNA(x), !anyNA(z))
  decomposition <- qr(z)
  basis <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
  x_res <- qr.resid(decomposition, x)
  sxx <- sum(x_res^2)
  resid <- qr.resid(decomposition, y)
  df <- length(y) - decomposition$rank - 1
  statistic <- function(shuffled) {
    numerator <- drop(crossprod(x_res, shuffled))
    rss <- colSums(shuffled^2) - colSums(crossprod(basis, shuffled)^2) - numerator^2 / sxx
    (numerator / sxx) / sqrt(rss / df / sxx)
  }
  observed <- statistic(matrix(resid, ncol = 1))
  chunks <- split(seq_len(ncol(perms)), ceiling(seq_len(ncol(perms)) / 5000))
  t_perm <- unlist(lapply(chunks, function(cols) statistic(matrix(resid[perms[, cols]], nrow = length(y)))), use.names = FALSE)
  list(t = observed, p = (1 + sum(abs(t_perm) >= abs(observed))) / (length(t_perm) + 1), t_perm = t_perm)
}

# --- follow-up of a hit -------------------------------------------------------------------------------

# Least-squares estimate of x adjusted for z (optionally weighted), with the two
# checks every hit gets: the estimate and p-value with each sample left out in
# turn, and a permutation p-value from shuffling x within `strata`.
# A shuffle that reproduces the observed arrangement ties with it, and rounding
# then decides whether it counts, differently on different processors. With
# `tie_tolerance` a shuffle within that relative distance below the observed
# statistic counts as a tie. The default of 0 is the comparison as stages 26 and
# 30 ran it.
followup_checks <- function(y, x, z, strata, weights = NULL, n_perm = 20000, seed = 1, tie_tolerance = 0) {
  ok <- is.finite(y) & is.finite(x) & rowSums(!is.finite(z)) == 0
  root <- if (is.null(weights)) rep(1, sum(ok)) else sqrt(weights[ok])
  y <- y[ok] * root
  z <- z[ok, , drop = FALSE] * root
  x <- x[ok]
  strata <- strata[ok]
  full <- fit_ols(y, x * root, z)
  loo <- vapply(seq_along(y), function(i) {
    one <- fit_ols(y[-i], (x * root)[-i], z[-i, , drop = FALSE])
    c(one$estimate, one$p)
  }, numeric(2))
  set.seed(seed)
  decomposition <- qr(z)
  y_res <- qr.resid(decomposition, y)
  shuffled <- qr.resid(decomposition, vapply(seq_len(n_perm), function(b) permute_within(x, strata), numeric(length(y))) * root)
  sxx <- colSums(shuffled^2)
  numerator <- drop(crossprod(shuffled, y_res))
  t_perm <- (numerator / sxx) / sqrt((sum(y_res^2) - numerator^2 / sxx) / (length(y) - decomposition$rank - 1) / sxx)
  list(
    estimate = full$estimate, se = full$se, p = full$p, n = full$n,
    permutation_p = (1 + sum(abs(t_perm) >= abs(full$t) * (1 - tie_tolerance), na.rm = TRUE)) / (n_perm + 1),
    loo_sign_stable = all(sign(loo[1, ]) == sign(full$estimate)), loo_max_p = max(loo[2, ]),
    loo_min_abs = min(abs(loo[1, ])), loo_max_abs = max(abs(loo[1, ]))
  )
}

# The pieces that followup_checks summarises, kept so that a statistic can be combined
# across sources: the least-squares t of x adjusted for z (optionally weighted), the same t
# with each sample left out, and under shuffles of x within strata. For one seed the
# shuffles are those of followup_checks. They are drawn `chunk` at a time, so a million of
# them never sit in memory together, and the order of the draws does not depend on the
# chunk. `loo = FALSE` skips the leave-one-out refits.
followup_statistics <- function(y, x, z, strata, weights = NULL, n_perm = 20000, seed = 1, chunk = 20000, loo = TRUE) {
  ok <- is.finite(y) & is.finite(x) & rowSums(!is.finite(z)) == 0
  root <- if (is.null(weights)) rep(1, sum(ok)) else sqrt(weights[ok])
  y <- y[ok] * root
  z <- z[ok, , drop = FALSE] * root
  x <- x[ok]
  strata <- strata[ok]
  full <- fit_ols(y, x * root, z)
  left_out <- if (loo) {
    vapply(seq_along(y), function(i) {
      one <- fit_ols(y[-i], (x * root)[-i], z[-i, , drop = FALSE])
      c(one$t, one$df)
    }, numeric(2))
  } else {
    matrix(NA_real_, 2, length(y))
  }
  set.seed(seed)
  decomposition <- qr(z)
  y_res <- qr.resid(decomposition, y)
  df <- length(y) - decomposition$rank - 1
  perm_t <- numeric(n_perm)
  for (start in seq(1, n_perm, by = chunk)) {
    draws <- start:min(n_perm, start + chunk - 1)
    shuffled <- qr.resid(decomposition, vapply(draws, function(b) permute_within(x, strata), numeric(length(y))) * root)
    sxx <- colSums(shuffled^2)
    numerator <- drop(crossprod(shuffled, y_res))
    perm_t[draws] <- (numerator / sxx) / sqrt((sum(y_res^2) - numerator^2 / sxx) / df / sxx)
  }
  list(estimate = full$estimate, t = full$t, df = full$df, n = full$n, loo_t = left_out[1, ], loo_df = left_out[2, ], perm_t = perm_t, perm_df = df)
}

# Follow-up of a statistic combined across sources. `parts` is a named list of
# followup_statistics results (NULL for a source without the feature) and `w` holds the
# combination weights under the same names. Each t becomes the z of the same two-sided
# p-value and the sources are pooled as stouffer does. The permutation p-value shuffles
# every source at once, and counts a shuffle within `tie_tolerance` below the observed
# statistic as a tie; leaving one sample out changes its own source only.
combine_followup <- function(parts, w, tie_tolerance = 1e-9) {
  parts <- parts[!vapply(parts, is.null, logical(1))]
  w <- unlist(w[names(parts)])
  scale <- sqrt(sum(w^2))
  as_z <- function(t, df) signed_z(2 * pt(-abs(t), df), t)
  z <- vapply(parts, function(part) as_z(part$t, part$df), numeric(1))
  observed <- sum(w * z) / scale
  n_perm <- min(vapply(parts, function(part) length(part$perm_t), numeric(1)))
  shuffled <- Reduce(`+`, lapply(names(parts), function(name) {
    one <- as_z(parts[[name]]$perm_t[seq_len(n_perm)], parts[[name]]$perm_df)
    w[[name]] * replace(one, !is.finite(one), 0)
  })) / scale
  loo <- unlist(lapply(names(parts), function(name) (sum(w * z) - w[[name]] * z[[name]] + w[[name]] * as_z(parts[[name]]$loo_t, parts[[name]]$loo_df)) / scale), use.names = FALSE)
  list(
    z = observed, p = 2 * pnorm(-abs(observed)), n_sources = length(parts), source_z = z,
    permutation_p = (1 + sum(abs(shuffled) >= abs(observed) * (1 - tie_tolerance))) / (n_perm + 1),
    loo_sign_stable = all(sign(loo) == sign(observed)), loo_max_p = max(2 * pnorm(-abs(loo))), loo_min_abs = min(abs(loo)), n_loo = length(loo)
  )
}

# --- report text -------------------------------------------------------------------------------------

# Text of each top-level section of a rendered report, named by its heading.
report_sections <- function(html_path) {
  clean <- function(x) trimws(gsub("\\s+", " ", gsub("<[^>]+>", " ", x)))
  html <- paste(readLines(html_path, warn = FALSE), collapse = " ")
  parts <- strsplit(html, "<h1[^>]*>", perl = TRUE)[[1]][-1]
  setNames(vapply(sub("^.*?</h1>", "", parts, perl = TRUE), clean, character(1), USE.NAMES = FALSE), clean(sub("</h1>.*$", "", parts)))
}

# The text of every paragraph, list item, table row, caption and heading of a rendered
# report, in document order. A list item that holds a paragraph is returned once.
report_blocks <- function(html_path) {
  doc <- xml2::read_html(html_path)
  block <- "self::p or self::li or self::tr or self::caption or self::h1 or self::h2 or self::h3 or self::h4"
  nodes <- xml2::xml_find_all(doc, sprintf("//body//*[%s][not(.//p or .//li or .//tr)]", block))
  cell_text <- function(node) if (xml2::xml_name(node) == "tr") paste(xml2::xml_text(xml2::xml_children(node)), collapse = " ") else xml2::xml_text(node)
  text <- trimws(gsub("\\s+", " ", vapply(nodes, cell_text, character(1))))
  text[nzchar(text)]
}

# The lines of one top-level section of an R Markdown file. Lines inside code chunks are
# never read as headings, so an R comment cannot end or start a section.
rmd_section <- function(lines, heading) {
  fence <- cumsum(grepl("^```", lines)) %% 2 == 1
  top <- which(grepl("^# ", lines) & !fence & !grepl("^```", lines))
  start <- top[lines[top] == paste("#", heading)]
  if (length(start) != 1) stop("section not found exactly once: ", heading)
  end <- c(top[top > start], length(lines) + 1)[1] - 1
  if (end <= start) character() else lines[(start + 1):end]
}

# Numbers typed into the prose of R Markdown lines: code chunks and inline R are dropped,
# and so are digits that are part of a name (TP53, BeatAML2, Figure 11a, t(8;21), inv(16)).
prose_numbers <- function(lines) {
  fence <- cumsum(grepl("^```", lines)) %% 2 == 1
  prose <- gsub("`[^`]*`", " ", lines[!fence & !grepl("^```", lines)])
  prose <- gsub("\\b(t|inv|del|dup)\\([0-9XY;.pq]+\\)", " ", prose)
  found <- regmatches(prose, gregexpr("(?<![A-Za-z0-9.,])[0-9](?:[0-9.,]*[0-9])?(?![A-Za-z0-9])", prose, perl = TRUE))
  unlist(found)
}

# --- driver status inside cancer types (stages 34 to 36) ---------------------------------------

# Status of each gene in each model: 1 altered (any listed event), 0 not altered, NA where
# the data that could show an alteration is missing. `events` holds one row per altered
# model and gene. A model without mutation calls is missing for every gene; `copy_missing`
# lists the model-gene pairs of the curated copy-number genes that have no copy ratio.
gene_status <- function(models, genes, events, mutation_profiled, copy_missing = NULL) {
  status <- matrix(0, length(models), length(genes), dimnames = list(models, genes))
  status[!models %in% mutation_profiled, ] <- NA
  at <- function(pairs) cbind(match(pairs$ModelID, models), match(pairs$gene, genes))
  if (!is.null(copy_missing)) {
    gaps <- copy_missing[ModelID %in% models & gene %in% genes]
    if (nrow(gaps) > 0) status[at(gaps)] <- NA
  }
  altered <- unique(events[ModelID %in% models & gene %in% genes, .(ModelID, gene)])
  if (nrow(altered) > 0) status[at(altered)] <- 1
  status
}

split_labels <- c("specificity testable", "Y effect among altered only", "below the gate", "one Y group missing")
split_label <- function(neg_alt, pos_alt, neg_wt, pos_wt, min_cell) {
  fifelse(pmin(neg_alt, pos_alt, neg_wt, pos_wt) >= min_cell, split_labels[1],
    fifelse(pmin(neg_alt, pos_alt) >= min_cell, split_labels[2],
      fifelse(pmin(neg_alt, pos_alt) == 0, split_labels[4], split_labels[3])))
}

# Lines per cell (Y status by gene status) for every gene and every cancer type with at
# least one altered line. `status` has one row per line and one column per gene; a line
# with a missing status counts in no cell, and a line without a type in no type.
split_census <- function(status, type, exposed, min_cell) {
  stopifnot(nrow(status) == length(type), length(type) == length(exposed))
  in_type <- !is.na(type) & type != ""
  types <- sort(unique(type[in_type]))
  member <- outer(type, types, function(a, b) !is.na(a) & a == b) * 1
  count <- function(mask, value) {
    cells <- crossprod(member * mask, (status == value & !is.na(status)) * 1)
    dimnames(cells) <- list(types, colnames(status))
    cells
  }
  cells <- list(neg_alt = count(exposed == 1, 1), pos_alt = count(exposed == 0, 1), neg_wt = count(exposed == 1, 0), pos_wt = count(exposed == 0, 0))
  keep <- which(cells$neg_alt + cells$pos_alt > 0, arr.ind = TRUE)
  out <- data.table(type = types[keep[, 1]], gene = colnames(status)[keep[, 2]])
  for (name in names(cells)) out[, (name) := as.numeric(cells[[name]][keep])]
  out[, label := split_label(neg_alt, pos_alt, neg_wt, pos_wt, min_cell)]
  out[order(type, gene)]
}

# For each split, the id of the first split with the same gene and the same lines (NA for
# that first one): a disease that is its whole lineage is fitted once.
same_lines_as <- function(id, gene, members) {
  key <- paste(gene, vapply(members, function(lines) paste(sort(lines), collapse = ","), character(1)))
  first <- id[match(key, key)]
  fifelse(first == id, NA_character_, first)
}

# Groups of sets that share at least `cut` of their members (Jaccard), joined through chains.
overlap_groups <- function(sets, cut) {
  group <- seq_along(sets)
  for (i in seq_along(sets)) {
    for (j in seq_len(i - 1)) {
      shared <- length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]]))
      if (shared >= cut) group[group == group[i]] <- group[j]
    }
  }
  match(group, unique(group))
}

# The lines of one cancer type with a known status, and the Y columns of the split model:
# Y among altered lines, Y among unaltered lines, and the interaction of the screen.
split_frame <- function(frame, in_type, status) {
  out <- copy(frame)[, ctx := as.numeric(status)]
  out <- out[in_type & !is.na(ctx)]
  out[, exposed_in := exposed * ctx]
  out[, exposed_out := exposed * (1 - ctx)]
  out[, interaction := exposed * ctx]
  out[]
}

# Share of a plain two-group comparison (n_a against n_b lines) that a design keeps for
# its tested column once the other columns are in the model.
information_kept <- function(design, tested, n_a, n_b) {
  (1 / n_a + 1 / n_b) / solve(crossprod(design))[tested, tested]
}

# The designs of the split model on one cancer type's lines:
#   score ~ library + status + Y among altered + Y among unaltered.
# The same model is written three ways so that each tested column comes first: `inside`
# (Y among altered lines), `outside` (Y among unaltered lines) and `interaction` (their
# difference, the screen's interaction). A split without unaltered lines keeps `inside`
# only. `specificity` says whether all four cells reach the gate.
split_designs <- function(sframe, covariates, min_cell, lineage = TRUE, min_residual_df = 3) {
  cells <- context_cells(sframe$exposed, sframe$ctx, "binary")
  both_out <- cells[["neg_out"]] >= 1 && cells[["pos_out"]] >= 1
  status_term <- if (cells[["neg_out"]] + cells[["pos_out"]] > 0) "ctx" else character()
  one <- function(tested, others) {
    design <- build_design(sframe, tested = tested, covariates = c(others, covariates), lineage = lineage)
    if (!design$estimable || nrow(design$design) - ncol(design$design) < min_residual_df) return(NULL)
    design
  }
  inside <- one("exposed_in", c(if (both_out) "exposed_out", status_term))
  if (is.null(inside)) stop("split cannot be fitted: the Y effect among altered lines is not estimable")
  list(
    cells = cells, specificity = all(cells >= min_cell), inside = inside,
    outside = if (both_out) one("exposed_out", c("exposed_in", status_term)),
    interaction = if (both_out) one("interaction", c("exposed", status_term))
  )
}

# The split model fitted as the screen fits a context (limma's moderated t, no trend on
# rank-normalised rows). `mat` holds features in rows and the screen's lines in columns.
# Specificity (the interaction) is reported only when all four cells reach the gate.
fit_split <- function(sframe, mat, covariates, min_cell, lineage = TRUE, trend = FALSE) {
  designs <- split_designs(sframe, covariates, min_cell, lineage = lineage)
  # A fresh table each time: the columns are renamed in place below.
  blank <- function() data.table(estimate = rep(NA_real_, nrow(mat)), se = NA_real_, p = NA_real_, df = NA_real_)
  scan_one <- function(design, tested) {
    if (is.null(design)) return(blank())
    scan <- limma_scan(mat[, sframe$ModelID[design$rows], drop = FALSE], design$design, tested, trend = trend)
    scan[, .(estimate, se = abs(estimate / t), p, df)]
  }
  inside <- scan_one(designs$inside, "exposed_in")
  outside <- scan_one(designs$outside, "exposed_out")
  interaction <- scan_one(if (designs$specificity) designs$interaction, "interaction")
  setnames(inside, paste0("inside_", names(inside)))
  setnames(outside, paste0("outside_", names(outside)))
  setnames(interaction, paste0("interaction_", names(interaction)))
  cells <- designs$cells
  cbind(
    data.table(feature = rownames(mat), n = length(designs$inside$rows), n_neg_in = cells[["neg_in"]], n_pos_in = cells[["pos_in"]], n_neg_out = cells[["neg_out"]], n_pos_out = cells[["pos_out"]]),
    inside, outside, interaction
  )
}

# What stays the same in a split's fit when one gene is added: the real genes' residual
# variances and degrees of freedom, and the algebra that gives a new gene's estimate and
# residual variance (as context_fit_parts() does for a context of the screen).
split_fit_parts <- function(sframe, mat, covariates, min_cell, lineage = TRUE) {
  designs <- split_designs(sframe, covariates, min_cell, lineage = lineage)
  one_part <- function(design, tested) {
    x <- design$design
    lines <- sframe$ModelID[design$rows]
    y <- mat[, lines, drop = FALSE]
    complete <- rowSums(is.na(y)) == 0
    sigma <- rep(NA_real_, nrow(y))
    df_residual <- rep(NA_real_, nrow(y))
    for (part in list(complete, !complete)) {
      if (!any(part)) next
      fit <- lmFit(y[part, , drop = FALSE], x)
      sigma[part] <- fit$sigma
      df_residual[part] <- fit$df.residual
    }
    cross_inverse <- solve(crossprod(x))
    list(
      lines = lines, sigma2 = sigma^2, df_residual = df_residual,
      weights = drop(cross_inverse[tested, , drop = FALSE] %*% t(x)), unscaled = sqrt(cross_inverse[tested, tested]),
      basis = qr.Q(qr(x)), df_new = nrow(x) - ncol(x)
    )
  }
  list(
    inside = one_part(designs$inside, "exposed_in"),
    interaction = if (designs$specificity && !is.null(designs$interaction)) one_part(designs$interaction, "interaction"),
    lines_all = colnames(mat), cells = designs$cells
  )
}

# Simulated genes tested as if each alone were added to the split's fit, on enlarged
# copies of the design. `simulated_raw` holds raw-scale values for every line of the
# screen, one row per copy; the rows of one `group` are the copies of one gene. Each row
# is rank-normalised over all lines, as in the screen. The copies share the tested
# coefficient and keep their own nuisance terms, so the estimate is the mean of the
# copies' estimates, its standard error falls as the square root of their number, and
# the residual variance pools the copies and their disagreement. The variance prior is
# that of today's real genes plus the one added gene (limma's robust squeezeVar).
enlarged_split_tests <- function(parts, simulated_raw, group) {
  stopifnot(identical(colnames(simulated_raw), parts$lines_all), !anyNA(simulated_raw), length(group) == nrow(simulated_raw))
  groups <- unique(group)
  copies <- tabulate(match(group, groups))
  if (length(unique(copies)) != 1) stop("every group needs the same number of copies")
  k <- copies[1]
  index <- match(group, groups)
  simulated <- t(apply(simulated_raw, 1, rank_normalise))
  colnames(simulated) <- parts$lines_all
  one_test <- function(part) {
    blank <- rep(NA_real_, length(groups))
    if (is.null(part)) return(list(estimate = blank, sigma2 = blank, df = blank, unscaled = blank, p = blank))
    y <- simulated[, part$lines, drop = FALSE]
    estimate_copy <- drop(y %*% part$weights)
    rss_copy <- rowSums(y^2) - rowSums((y %*% part$basis)^2)
    estimate <- drop(rowsum(estimate_copy, index, reorder = TRUE)) / k
    between <- drop(rowsum((estimate_copy - estimate[index])^2, index, reorder = TRUE)) / part$unscaled^2
    df <- k * part$df_new + (k - 1)
    sigma2 <- (drop(rowsum(rss_copy, index, reorder = TRUE)) + between) / df
    unscaled <- part$unscaled / sqrt(k)
    df_pooled <- sum(part$df_residual, na.rm = TRUE) + df
    last <- length(part$sigma2) + 1
    p <- vapply(seq_along(groups), function(i) {
      squeezed <- squeezeVar(c(part$sigma2, sigma2[i]), c(part$df_residual, df), robust = TRUE)
      df_prior <- if (length(squeezed$df.prior) > 1) squeezed$df.prior[last] else squeezed$df.prior
      2 * pt(-abs(estimate[i] / unscaled / sqrt(squeezed$var.post[last])), min(df + df_prior, df_pooled))
    }, numeric(1))
    list(estimate = unname(estimate), sigma2 = unname(sigma2), df = rep(df, length(groups)), unscaled = rep(unscaled, length(groups)), p = p)
  }
  inside <- one_test(parts$inside)
  interaction <- one_test(parts$interaction)
  data.table(
    group = groups, copies = k,
    inside_estimate = inside$estimate, inside_sigma2 = inside$sigma2, inside_df = inside$df, inside_unscaled = inside$unscaled, inside_p = inside$p,
    interaction_estimate = interaction$estimate, interaction_p = interaction$p
  )
}

# One copy per gene: the p-values each simulated gene would get if it alone were added
# to the split's fit as it stands.
one_more_split_tests <- function(parts, simulated_raw) {
  out <- enlarged_split_tests(parts, simulated_raw, group = seq_len(nrow(simulated_raw)))
  out[, .(row = group, inside_p, interaction_p)]
}

# The multiple of today's lines at which an effect is found: the first multiple whose
# detection limit is at or below `effect`, interpolated between neighbouring multiples on
# the log scale. 1 when the effect is found today, Inf when no multiple reaches it. A
# limit of NA means that no shift was found often enough at that multiple.
lines_needed <- function(multiple, limit, effect) {
  stopifnot(length(multiple) == length(limit), !is.unsorted(multiple), multiple[1] == 1)
  vapply(effect, function(size) {
    reached <- which(!is.na(limit) & limit <= size)
    if (length(reached) == 0) return(Inf)
    first <- reached[1]
    if (first == 1) return(1)
    if (is.na(limit[first - 1])) return(multiple[first])
    slope <- (log(multiple[first]) - log(multiple[first - 1])) / (log(limit[first]) - log(limit[first - 1]))
    exp(log(multiple[first - 1]) + slope * (log(size) - log(limit[first - 1])))
  }, numeric(1))
}

# TRUE where lines_needed() could not interpolate, because no limit exists at the next
# smaller multiple: the multiple it returns is then an upper bound.
lines_needed_bound <- function(multiple, limit, effect) {
  vapply(effect, function(size) {
    reached <- which(!is.na(limit) & limit <= size)
    length(reached) > 0 && reached[1] > 1 && is.na(limit[reached[1] - 1])
  }, logical(1))
}

# The p-value the follow-up's shuffles give on average for the most extreme arrangement of
# Y status among the altered lines: the share of shuffles that land on that arrangement
# (or, when the two groups are equal in size, on its mirror image, which is as extreme in
# a design without covariates), plus the one added to the count. It is an average, not a
# bound: 20,000 shuffles land on one arrangement a random number of times, and with
# covariates the mirror image can count for less.
permutation_floor <- function(n_a, n_b, n_perm) (1 + (n_a == n_b)) / choose(n_a + n_b, n_a) + 1 / (n_perm + 1)

# Row order by p-value (within `first`, when given) in which tests that agree to `digits`
# significant digits keep the order given by `tie`. Two levels of one cancer type can give
# the same test, and which of two such p-values is the smaller in its last bit depends on
# the processor.
tied_order <- function(p, tie, first = NULL, digits = 10) {
  if (is.null(first)) order(signif(p, digits), tie, method = "radix") else order(first, signif(p, digits), tie, method = "radix")
}

# The split model by plain least squares, each gene with its own residual variance: the
# p-value of the Y effect among altered lines and, with four cells at the gate, of its
# difference from the unaltered lines. A gene with missing values is fitted on its own lines.
fit_split_ols <- function(sframe, mat, covariates, min_cell, lineage = TRUE) {
  designs <- split_designs(sframe, covariates, min_cell, lineage = lineage)
  one <- function(design, tested) {
    out <- rep(NA_real_, nrow(mat))
    if (is.null(design)) return(out)
    x <- design$design
    y <- mat[, sframe$ModelID[design$rows], drop = FALSE]
    nuisance <- x[, colnames(x) != tested, drop = FALSE]
    complete <- rowSums(is.na(y)) == 0
    if (any(complete)) out[complete] <- 2 * pt(-abs(fast_t(y[complete, , drop = FALSE], x[, tested], nuisance)), nrow(x) - ncol(x))
    for (i in which(!complete)) out[i] <- fit_ols(y[i, ], x[, tested], nuisance)$p
    out
  }
  data.table(feature = rownames(mat), inside_p_ols = one(designs$inside, "exposed_in"), interaction_p_ols = one(if (designs$specificity) designs$interaction, "interaction"))
}

# Plain least-squares p-values for rows of values on any scale, from a split's fit parts:
# no ranking and no borrowed variance. Used to ask what a test on the Chronos scale would
# find in a split.
ols_split_tests <- function(parts, values) {
  stopifnot(identical(colnames(values), parts$lines_all), !anyNA(values))
  one <- function(part) {
    if (is.null(part)) return(rep(NA_real_, nrow(values)))
    y <- values[, part$lines, drop = FALSE]
    estimate <- drop(y %*% part$weights)
    sigma2 <- (rowSums(y^2) - rowSums((y %*% part$basis)^2)) / part$df_new
    2 * pt(-abs(estimate / part$unscaled / sqrt(sigma2)), part$df_new)
  }
  data.table(row = seq_len(nrow(values)), inside_p = one(parts$inside), interaction_p = one(parts$interaction))
}
