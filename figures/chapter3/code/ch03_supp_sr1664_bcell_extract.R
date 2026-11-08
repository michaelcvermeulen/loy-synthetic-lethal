# Extract the observed SR-1664 within-subset comparison without refitting models.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
root <- dirname(dirname(dirname(dirname(script))))
derived <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07/derived/loy_reanalysis_v1")
out <- file.path(root, "figures/chapter3/source_data/sr1664_bcell_response")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(data.table))
hits <- fread(file.path(root, "tables/chapter3/depmap_drug_within_subsets/within_subset_hits.csv"), encoding = "UTF-8")
hit <- hits[compound_name == "SR-1664" & dataset == "prism_REP1M" & context_id == "disease__Mature_B_Cell_Neoplasms"]
stopifnot(nrow(hit) == 1, hit$n_neg_in == 8, hit$n_pos_in == 16)
cohort <- as.data.table(readRDS(file.path(derived, "cohort.rds")))
context <- as.data.table(readRDS(file.path(derived, "context_matrix.rds")))
mat <- readRDS(file.path(derived, "matrices/prism_REP1M.rds"))
stopifnot(!anyDuplicated(cohort$ModelID), !anyDuplicated(context$ModelID), !anyDuplicated(rownames(mat)))
frame <- cohort[ModelID %in% rownames(mat) & !is.na(y_call)]
stopifnot(all(frame$y_call %in% c("Y_negative", "Y_positive")))
frame[, response := mat[ModelID, hit$feature]]
frame[, in_subset := context[[hit$context_id]][match(ModelID, context$ModelID)]]
points <- frame[in_subset == 1 & is.finite(response)]
points[, y_status := fifelse(y_call == "Y_negative", "LOY", "Y-positive")]
stopifnot(nrow(points) == 24, all(points$Sex == "Male"),
          points[y_status == "LOY", .N] == hit$n_neg_in,
          points[y_status == "Y-positive", .N] == hit$n_pos_in,
          uniqueN(points$OncotreeLineage) == 1)
# This subset has one lineage and no screening-round term, so its raw-scale
# saved coefficient must equal the observed group-mean difference.
delta <- points[y_status == "LOY", mean(response)] - points[y_status == "Y-positive", mean(response)]
stopifnot(abs(delta - hit$inside) < 1e-10)
fwrite(points[, .(ModelID, cell_line = StrippedCellLineName, PatientID, y_status, response,
                 lineage = OncotreeLineage, disease = OncotreePrimaryDisease, subtype = OncotreeSubtype)],
       file.path(out, "plotted_points.csv"))
fwrite(hit, file.path(out, "saved_test.csv"))
summary <- points[, .(models = .N, mean = mean(response), median = median(response),
                     minimum = min(response), maximum = max(response)), by = y_status]
fwrite(summary, file.path(out, "group_summary.csv"))
print(summary)
cat("Extracted 24 measured models; raw mean difference matches the saved effect. No models refitted.\n")
