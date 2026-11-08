# Descriptive RCOR3 follow-up: copy measurements and existing test results, without refitting.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
derived <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1")
out <- file.path(root, "figures/chapter3/source_data/rcor3_exploratory")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cohort <- read.csv(file.path(derived, "cohort.csv"), check.names = FALSE)
points <- list()
for (assay in c("CRISPR", "RNAi")) {
  matrix_name <- if (assay == "CRISPR") "crispr_effect" else "rnai"
  feature <- if (assay == "CRISPR") "RCOR3 (55758)" else "RCOR3"
  cat("Reading", assay, "matrix\n"); flush.console()
  mat <- readRDS(file.path(derived, "matrices", paste0(matrix_name, ".rds")))
  stopifnot(feature %in% colnames(mat))
  frame <- cohort[which(cohort$in_primary & cohort$ModelID %in% rownames(mat)), ]
  frame$score <- mat[match(frame$ModelID, rownames(mat)), feature]
  frame <- frame[is.finite(frame$score) & !is.na(frame$y_negative), ]
  points[[assay]] <- data.frame(assay=assay, ModelID=frame$ModelID, PatientID=frame$PatientID,
    cell_line=frame$StrippedCellLineName, lineage=frame$OncotreeLineage,
    disease=frame$OncotreePrimaryDisease, subtype=frame$OncotreeSubtype,
    y_status=ifelse(frame$y_negative == 1, "LOY", "Y-positive"), score=frame$score,
    LibraryKY=frame$LibraryKY, LibraryHumagneCD=frame$LibraryHumagneCD,
    WGD=frame$WGD, tp53_altered=frame$tp53_altered, CIN=frame$CIN,
    Aneuploidy=frame$Aneuploidy, Ploidy=frame$Ploidy,
    ScreenNNMD=frame$ScreenNNMD, ScreenDoublingTime=frame$ScreenDoublingTime)
  rm(mat); gc()
}
dat <- do.call(rbind, points)
stopifnot(nrow(points$CRISPR) == 568, nrow(points$RNAi) == 306)
pan <- read.csv(file.path(derived, "gene_pan_results.csv"), check.names=FALSE)
pan <- pan[pan$symbol == "RCOR3", ]
stopifnot(nrow(pan) == 1)
stopifnot(sum(points$CRISPR$y_status == "LOY") == pan$n_negative,
          sum(points$CRISPR$y_status == "Y-positive") == pan$n_positive,
          abs(mean(points$CRISPR$score[points$CRISPR$y_status == "LOY"]) - pan$mean_negative) < 1e-10,
          abs(mean(points$CRISPR$score[points$CRISPR$y_status == "Y-positive"]) - pan$mean_positive) < 1e-10)
probability <- readRDS(file.path(derived, "matrices/crispr_dependency.rds"))
points$CRISPR$dependency_probability <- probability[match(points$CRISPR$ModelID, rownames(probability)), "RCOR3 (55758)"]
stopifnot(sum(points$CRISPR$dependency_probability[points$CRISPR$y_status == "LOY"] >= .5, na.rm=TRUE) == 0,
          sum(points$CRISPR$dependency_probability[points$CRISPR$y_status == "Y-positive"] >= .5, na.rm=TRUE) == 1)
write.csv(points$CRISPR, file.path(out, "crispr_points.csv"), row.names=FALSE)
write.csv(points$RNAi, file.path(out, "rnai_points.csv"), row.names=FALSE)
for (name in c("gene_pan_results", "gene_pan_continuous", "gene_pan_rnai", "forest_confounders")) {
  tab <- read.csv(file.path(derived, paste0(name, ".csv")), check.names=FALSE)
  write.csv(tab[tab$symbol == "RCOR3", ], file.path(out, paste0("saved_", name, ".csv")), row.names=FALSE)
}
cat("Most negative CRISPR scores:\n")
print(head(points$CRISPR[order(points$CRISPR$score), c("cell_line", "lineage", "y_status", "score", "dependency_probability")], 12), row.names=FALSE)
cat("Group sizes:\n")
print(table(dat$assay, dat$y_status))
