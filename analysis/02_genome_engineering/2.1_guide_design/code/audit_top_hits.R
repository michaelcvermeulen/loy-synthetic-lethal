# Audit saved outputs only; does not reconstruct the missing shortlist filter.
# Rscript --vanilla audit_top_hits.R <sgRNA_guide_development> <site_level_RDS>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply project folder and site-level RDS path")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, "filter_audit")
dir.create(out, showWarnings = FALSE)
candidates <- read.csv(file.path(root, "candidate_guides.csv"), check.names = FALSE)
hits <- read.csv(file.path(root, "top_hits.csv"), check.names = FALSE)
condensed <- read.csv(file.path(root, "top_hits_condensed.csv"), check.names = FALSE)
sites <- readRDS(args[[2]])
stopifnot(!anyDuplicated(candidates$guide), !anyDuplicated(hits$guide),
          identical(hits[, names(condensed)], condensed))
eligible <- candidates[candidates$sumOther == 0 & candidates$sumY > 9, ]
stopifnot(all(hits$guide %in% eligible$guide))
joined <- candidates[match(hits$guide, candidates$guide), ]
stopifnot(all(hits$chrY_exact_hits == joined$sumY),
          all(hits$chrY_exact_forward == joined$chrY),
          all(hits$chrY_exact_reverse == joined$rev_chrY))

summarise_sites <- function(z) {
  z <- as.data.frame(z)
  stopifnot(length(unique(z$main_guide)) == 1L)
  data.frame(guide = unique(z$main_guide), site_rows = nrow(z),
             mismatch_1 = sum(z$mismatches == 1),
             mismatch_2 = sum(z$mismatches == 2),
             mismatch_3 = sum(z$mismatches == 3),
             mismatch_4 = sum(z$mismatches == 4),
             mismatch_mean_CFD = mean(z$CFD),
             mismatch_max_CFD = max(z$CFD),
             mismatch_sum_CFD = sum(z$CFD),
             chrY_rows = sum(grepl("^chrY( |$)", z$chr)),
             chrM_rows = sum(grepl("^chrM( |$)", z$chr)))
}
site_summary <- do.call(rbind, lapply(sites, summarise_sites))
stopifnot(setequal(site_summary$guide, hits$guide))
site_summary <- site_summary[match(hits$guide, site_summary$guide), ]
write.csv(site_summary, file.path(out, "saved_site_level_summary.csv"), row.names = FALSE)
fields <- c("mismatch_1", "mismatch_2", "mismatch_3", "mismatch_4",
            "mismatch_mean_CFD", "mismatch_max_CFD", "mismatch_sum_CFD")
checks <- do.call(rbind, lapply(fields, function(f) {
  delta <- site_summary[[f]] - hits[[f]]
  data.frame(field = f, max_absolute_difference = max(abs(delta)),
             rows_different_beyond_1e_minus_9 = sum(abs(delta) > 1e-9))
}))
write.csv(checks, file.path(out, "site_to_top_hits_checks.csv"), row.names = FALSE)
differences <- do.call(rbind, lapply(fields, function(f) {
  select <- abs(site_summary[[f]] - hits[[f]]) > 1e-9
  if (!any(select)) return(NULL)
  data.frame(guide = hits$guide[select], field = f,
             saved_top_hits_value = hits[[f]][select],
             site_level_value = site_summary[[f]][select])
}))
write.csv(differences, file.path(out, "site_to_top_hits_differences.csv"), row.names = FALSE)
chr <- unlist(lapply(sites, function(z) sub(" .*", "", z$chr)))
chromosome_counts <- as.data.frame(table(chr))
names(chromosome_counts) <- c("chromosome", "site_rows")
write.csv(chromosome_counts, file.path(out, "site_chromosome_counts.csv"), row.names = FALSE)
observations <- data.frame(
  observation = c("candidate_guides", "sumOther_zero", "sumOther_zero_and_sumY_gt9",
                  "top_hits", "top_hits_without_1_or_2_mismatch_sites", "top_hits_with_exactly10_Y_hits",
                  "top_hits_with_at_most20_Y_hits", "site_rows", "site_chrY_rows",
                  "top_hits_sorted_by_decreasing_exact_Y_hits"),
  value = c(nrow(candidates), sum(candidates$sumOther == 0), nrow(eligible), nrow(hits),
            sum(hits$mismatch_1 == 0 & hits$mismatch_2 == 0), sum(hits$chrY_exact_hits == 10),
            sum(hits$chrY_exact_hits <= 20), sum(site_summary$site_rows), sum(site_summary$chrY_rows),
            all(diff(hits$chrY_exact_hits) <= 0)))
write.csv(observations, file.path(out, "filter_observations.csv"), row.names = FALSE)
print(observations, row.names = FALSE)
print(checks, row.names = FALSE)
