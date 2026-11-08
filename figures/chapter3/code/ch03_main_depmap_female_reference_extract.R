# Extract female reference values for F01e from the same frozen 26Q1 matrix.
# No dependency models, effects or P values are fitted here.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
root <- dirname(dirname(dirname(dirname(script))))
export <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07")
derived <- file.path(export, "derived/loy_reanalysis_v1")
out <- file.path(root, "figures/chapter3/source_data/depmap_y_status_dependencies")
metadata <- read.csv(file.path(derived, "y_status_calls.csv"), stringsAsFactors = FALSE, check.names = FALSE)
effect <- readRDS(file.path(derived, "matrices/crispr_effect.rds"))
stopifnot(is.matrix(effect), !anyDuplicated(rownames(effect)), !anyDuplicated(metadata$ModelID))
symbols <- sub(" \\(.*$", "", colnames(effect))
genes <- c("DDX3X", "EIF1AX", "ZFX", "RPS4X")
stopifnot(all(vapply(genes, function(g) sum(symbols == g) == 1L, logical(1))))

# Match the cancer and engineering exclusions used for the male cohort.
# Female eligibility does not depend on RNA availability or a LOY call.
female <- metadata[!is.na(metadata$Sex) & metadata$Sex == "Female" &
                   !is.na(metadata$OncotreePrimaryDisease) & metadata$OncotreePrimaryDisease != "Non-Cancerous" &
                   (is.na(metadata$EngineeredModel) | metadata$EngineeredModel == ""), ]
female <- female[female$ModelID %in% rownames(effect), ]
female <- female[order(female$ModelID), ]
rows <- do.call(rbind, lapply(genes, function(gene) {
    values <- effect[female$ModelID, which(symbols == gene)]
    data.frame(symbol = gene, ModelID = female$ModelID, status = "Female", chronos = values,
               PatientID = female$PatientID, StrippedCellLineName = female$StrippedCellLineName,
               Sex = female$Sex, OncotreeLineage = female$OncotreeLineage,
               OncotreePrimaryDisease = female$OncotreePrimaryDisease,
               EngineeredModel = female$EngineeredModel, ModelType = female$ModelType,
               stringsAsFactors = FALSE)
}))
rows <- rows[is.finite(rows$chronos), ]
stopifnot(!anyDuplicated(rows[c("symbol", "ModelID")]))

# Confirm that the frozen matrix reproduces every displayed male measurement.
male <- read.csv(file.path(out, "fig2b_hit_boxplots.csv"), stringsAsFactors = FALSE)
matrix_values <- effect[cbind(match(male$ModelID, rownames(effect)), match(male$symbol, symbols))]
stopifnot(all(is.finite(matrix_values)), max(abs(matrix_values - male$chronos)) < 1e-12)
stopifnot(!any(rows$ModelID %in% male$ModelID))
write.csv(rows, file.path(out, "panel_e_female_reference.csv"), row.names = FALSE, na = "")
summary <- do.call(rbind, lapply(split(rows, rows$symbol), function(x) {
    data.frame(symbol = x$symbol[1], n = nrow(x), n_donors = length(unique(x$PatientID)),
               mean = mean(x$chronos), median = median(x$chronos), min = min(x$chronos), max = max(x$chronos))
}))
write.csv(summary, file.path(out, "panel_e_female_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
cat("Verified all", nrow(male), "male values against the same frozen matrix; no statistical models refitted.\n")
