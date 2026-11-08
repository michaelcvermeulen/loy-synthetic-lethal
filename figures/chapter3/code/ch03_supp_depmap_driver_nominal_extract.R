# Extract the two non-core tests with P < 1e-5, using frozen models and values.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
root <- dirname(dirname(dirname(dirname(script))))
derived <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1")
out <- file.path(root, "figures/chapter3/source_data/depmap_driver_nominal_candidates")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
screen_path <- file.path(derived, "final/driver_type_screen.csv.gz")
columns <- names(read.csv(gzfile(screen_path), nrows = 1, check.names = FALSE))
keep <- c("split_id", "feature", "symbol", "type", "gene", "inside", "inside_se", "inside_p", "inside_fdr",
          "inside_fdr_split", "outside", "outside_p", "interaction", "interaction_p", "interaction_fdr",
          "n_neg_in", "n_pos_in", "n_neg_out", "n_pos_out", "altered_set")
classes <- setNames(rep("NULL", length(columns)), columns)
classes[keep] <- NA_character_
screen <- read.csv(gzfile(screen_path), colClasses = classes, check.names = FALSE)
stopifnot(nrow(screen) == 414052)
nominal <- screen[is.finite(screen$inside_p) & screen$inside_p < 1e-5, ]
distinct <- nominal[!duplicated(nominal[c("symbol", "altered_set")]), ]
core <- c("DDX3X", "EIF1AX", "RPS4X", "ZFX")
stopifnot(nrow(nominal) == 17, nrow(distinct) == 14, sum(distinct$symbol %in% core) == 12)
saved <- distinct[!distinct$symbol %in% core, ]
saved <- saved[match(c("INSM1", "ATRX"), saved$symbol), ]
stopifnot(nrow(saved) == 2, all(saved$inside > 0), all(saved$inside_fdr > .05))
write.csv(saved, file.path(out, "saved_tests.csv"), row.names = FALSE, na = "")
write.csv(distinct, file.path(out, "distinct_tests_p_below_1e5.csv"), row.names = FALSE, na = "")

lines <- read.csv(file.path(derived, "final/driver_type_split_lines.csv"), check.names = FALSE)
cohort <- read.csv(file.path(derived, "cohort.csv"), check.names = FALSE)
effect <- readRDS(file.path(derived, "matrices/crispr_effect.rds"))
stopifnot(!anyDuplicated(cohort$ModelID), !anyDuplicated(rownames(effect)))
points <- do.call(rbind, lapply(seq_len(nrow(saved)), function(i) {
    hit <- saved[i, ]
    selected <- lines[lines$split_id == hit$split_id, ]
    meta <- cohort[match(selected$ModelID, cohort$ModelID), ]
    stopifnot(!anyNA(meta$ModelID), all(selected$exposed == meta$y_negative), all(meta$Sex == "Male"))
    selected$chronos <- effect[selected$ModelID, hit$feature]
    selected <- selected[is.finite(selected$chronos), ]
    stopifnot(sum(selected$ctx == 1 & selected$exposed == 1) == hit$n_neg_in,
              sum(selected$ctx == 1 & selected$exposed == 0) == hit$n_pos_in,
              sum(selected$ctx == 0 & selected$exposed == 1) == hit$n_neg_out,
              sum(selected$ctx == 0 & selected$exposed == 0) == hit$n_pos_out)
    meta <- cohort[match(selected$ModelID, cohort$ModelID), ]
    data.frame(symbol = hit$symbol, driver = hit$gene, cancer = hit$type,
               selected, y_status = ifelse(selected$exposed == 1, "LOY", "Y-positive"),
               cell_line = meta$StrippedCellLineName, PatientID = meta$PatientID,
               lineage = meta$OncotreeLineage, disease = meta$OncotreePrimaryDisease,
               subtype = meta$OncotreeSubtype, LibraryKY = meta$LibraryKY,
               LibraryHumagneCD = meta$LibraryHumagneCD, stringsAsFactors = FALSE)
}))
stopifnot(!anyDuplicated(points[c("symbol", "ModelID")]))
write.csv(points, file.path(out, "all_context_points.csv"), row.names = FALSE, na = "")
plot_points <- points
stopifnot(nrow(plot_points) == 106,
          sum(plot_points$symbol == "INSM1" & plot_points$ctx == 0 & plot_points$exposed == 0) == 21,
          sum(plot_points$symbol == "INSM1" & plot_points$ctx == 0 & plot_points$exposed == 1) == 31,
          sum(plot_points$symbol == "ATRX" & plot_points$ctx == 0 & plot_points$exposed == 0) == 9,
          sum(plot_points$symbol == "ATRX" & plot_points$ctx == 0 & plot_points$exposed == 1) == 0)
write.csv(plot_points, file.path(out, "plotted_points.csv"), row.names = FALSE, na = "")
print(saved[, c("symbol", "gene", "type", "n_pos_in", "n_neg_in", "inside", "inside_p", "inside_fdr")], row.names = FALSE)
cat("Extracted 106 plotted points for all driver-by-Y groups; the CDKN2A-unaltered LOY group is empty. No models refitted.\n")
