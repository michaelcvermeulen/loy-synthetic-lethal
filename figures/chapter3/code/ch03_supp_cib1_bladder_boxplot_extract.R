# Extract measured CIB1 scores for a descriptive preview of the retained within-cancer association.
# Use the saved discovery statistics; do not rerun the screen or change its testing family.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
derived <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1")
out <- file.path(root, "figures/chapter3/source_data/cib1_bladder_boxplot")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cohort <- read.csv(file.path(derived, "cohort.csv"), check.names = FALSE)
hits <- read.csv(file.path(derived, "power_alone_hits.csv"), check.names = FALSE)
hit <- hits[hits$symbol == "CIB1" & hits$context_id == "lineage__Bladder_Urinary_Tract", ]
stopifnot(nrow(hit) == 1)
mat <- readRDS(file.path(derived, "matrices/crispr_effect.rds"))
feature <- hit$feature
stopifnot(feature %in% colnames(mat))
frame <- cohort[which(cohort$in_primary & cohort$OncotreeLineage == "Bladder/Urinary Tract" &
                     cohort$ModelID %in% rownames(mat)), ]
frame$chronos <- mat[match(frame$ModelID, rownames(mat)), feature]
frame <- frame[complete.cases(frame[, c("chronos", "y_negative", "LibraryKY", "LibraryHumagneCD")]), ]
stopifnot(nrow(frame) == 18, sum(frame$y_negative == 1) == hit$n_negative,
          sum(frame$y_negative == 0) == hit$n_positive)
# All library covariates are constant here, so the raw-scale fitted coefficient is the mean difference.
stopifnot(length(unique(frame$LibraryKY)) == 1, length(unique(frame$LibraryHumagneCD)) == 1)
delta <- mean(frame$chronos[frame$y_negative == 1]) - mean(frame$chronos[frame$y_negative == 0])
stopifnot(abs(delta - hit$inside) < 1e-9)
points <- data.frame(ModelID=frame$ModelID, PatientID=frame$PatientID,
                     cell_line=frame$StrippedCellLineName, lineage=frame$OncotreeLineage,
                     y_status=ifelse(frame$y_negative == 1, "LOY", "Y-positive"),
                     chronos=frame$chronos, LibraryKY=frame$LibraryKY, LibraryHumagneCD=frame$LibraryHumagneCD)
write.csv(points, file.path(out, "plotted_points.csv"), row.names=FALSE)
write.csv(hit[, c("context_id", "label", "feature", "symbol", "n_negative", "n_positive", "inside", "inside_p", "own_fdr")],
          file.path(out, "saved_within_cancer_result.csv"), row.names=FALSE)
cat("Verified raw mean difference:", delta, "\n")
print(points[, c("cell_line", "y_status", "chronos")], row.names=FALSE)
