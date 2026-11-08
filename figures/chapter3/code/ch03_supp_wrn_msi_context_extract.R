# Direct WRN-MSI follow-up: extract the frozen measurements and benchmark cohort.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
root <- dirname(dirname(dirname(dirname(script))))
derived <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1")
out <- file.path(root, "figures/chapter3/source_data/wrn_msi_context_followup")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cohort <- read.csv(file.path(derived, "cohort.csv"), check.names = FALSE)
context <- read.csv(file.path(derived, "context_matrix.csv"), check.names = FALSE)
effect <- readRDS(file.path(derived, "matrices/crispr_effect.rds"))
stopifnot(!anyDuplicated(cohort$ModelID), !anyDuplicated(context$ModelID))
cohort$MSI <- context$genomic__MSI[match(cohort$ModelID, context$ModelID)]
cohort$WRN <- effect[match(cohort$ModelID, rownames(effect)), "WRN (7486)"]
keep <- c("ModelID", "PatientID", "StrippedCellLineName", "Sex", "OncotreeLineage",
          "OncotreePrimaryDisease", "OncotreeSubtype", "y_call", "MSIScore", "MSI", "WRN",
          "LibraryKY", "LibraryHumagneCD", "in_primary")
points <- cohort[is.finite(cohort$WRN), keep]
points$benchmark_eligible <- complete.cases(points[, c("MSI", "WRN", "LibraryKY", "LibraryHumagneCD")])
points$y_called <- !is.na(points$y_call) & points$y_call != ""
stopifnot(all(points$Sex == "Male"), sum(points$benchmark_eligible & points$MSI == 1, na.rm = TRUE) == 27,
          sum(points$benchmark_eligible & points$MSI == 0, na.rm = TRUE) == 617)
write.csv(points, file.path(out, "wrn_msi_points.csv"), row.names = FALSE, na = "")
print(with(points, table(benchmark_eligible, y_called, MSI, useNA = "ifany")))
cat("Extracted frozen WRN measurements and MSI calls; no models fitted.\n")
