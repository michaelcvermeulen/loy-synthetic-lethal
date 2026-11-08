# Stage 18 (plan Stage 7c): does Y status keep its importance for a hit once a
# forest can also use lineage, library, WGD, TP53, aneuploidy, LOH, ploidy, MSI,
# CIN and screen quality?
# The real Y label competes with shadow labels shuffled within lineage x WGD x
# TP53 strata, which carry the same confounder structure. A hit "survives
# flexible adjustment" if the real label beats every shadow in at least 80% of
# donor subsamples under both seeds. The same procedure on random non-hit genes
# gives the false-survival rate.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "forest_helpers.R"))
assert_freeze()

cohort <- read_result("cohort")
gene_ann <- read_result("gene_annotation")
effect <- readRDS(file.path(matrix_dir, "crispr_effect.rds"))
pan_hits <- read_result("gene_pan_hits")
continuous <- read_result("gene_pan_continuous")
drug_hits <- read_result("drug_pan_hits")

frame <- cohort[in_primary == TRUE & ModelID %in% rownames(effect)]
numeric_covariates <- c(
  "LibraryKY", "LibraryHumagneCD", "WGD", "tp53_altered", "Aneuploidy", "LoHFraction", "Ploidy", "MSIScore", "CIN",
  "ScreenNNMD", "ScreenROCAUC", "CasActivity", "ScreenDoublingTime", "ScreenMedianEssentialDepletion"
)
covariates <- as.matrix(frame[, numeric_covariates, with = FALSE])
missing <- is.na(covariates)
indicators <- missing[, colSums(missing) > 0, drop = FALSE] * 1
colnames(indicators) <- paste0("missing_", colnames(indicators))
for (j in seq_len(ncol(covariates))) covariates[missing[, j], j] <- median(covariates[, j], na.rm = TRUE)
features <- cbind(y_status = frame$y_negative, lineage_columns(frame$OncotreeLineage), covariates, indicators)

# Strata for the shadow labels; a stratum too small to shuffle falls back to its lineage.
fine <- paste(frame$OncotreeLineage, frame$WGD, frame$tp53_altered)
size <- table(fine)[fine]
strata <- ifelse(size >= 4, fine, paste(frame$OncotreeLineage, "other"))
log_msg("confounder check: %d lines, %d features, %d shadow strata (%d lines in strata of one)", nrow(frame), ncol(features), length(unique(strata)), sum(table(strata)[strata] == 1))

donors <- unique(frame$PatientID)
check_outcome <- function(feature, seed_index) {
  d <- effect[frame$ModelID, feature]
  set.seed(forest_par$seeds[seed_index] + sum(utf8ToInt(feature)))
  wins <- vapply(seq_len(forest_par$confounder_subsamples), function(b) {
    keep <- frame$PatientID %in% sample(donors, round(forest_par$confounder_subsample_fraction * length(donors))) & is.finite(d)
    shadows <- vapply(seq_len(forest_par$confounder_shadows), function(s) permute_within(frame$y_negative[keep], strata[keep]), numeric(sum(keep)))
    colnames(shadows) <- paste0("shadow_", seq_len(ncol(shadows)))
    fit <- ranger(x = cbind(features[keep, , drop = FALSE], shadows), y = d[keep], num.trees = forest_par$trees, importance = "permutation", num.threads = 1, seed = b)
    importance <- fit$variable.importance
    c(importance[["y_status"]] > max(importance[colnames(shadows)]), importance[["y_status"]], rank(-importance)[["y_status"]])
  }, numeric(3))
  c(success = mean(wins[1, ]), importance = mean(wins[2, ]), rank = median(wins[3, ]))
}

hit_features <- union(pan_hits$feature, continuous[fdr < params$gene_scan$fdr, feature])
set.seed(params$seed + 7)
pool <- setdiff(gene_ann[in_scan == TRUE & is.na(y_partner), feature], hit_features)
pool <- pool[colSums(is.na(effect[frame$ModelID, pool])) == 0]
null_features <- sample(pool, forest_par$confounder_null_genes)
outcomes <- data.table(feature = c(hit_features, null_features), role = rep(c("hit", "random non-hit"), c(length(hit_features), length(null_features))))

results <- rbindlist(parallel::mclapply(seq_len(nrow(outcomes)), function(i) {
  seed1 <- check_outcome(outcomes$feature[i], 1)
  seed2 <- check_outcome(outcomes$feature[i], 2)
  data.table(
    feature = outcomes$feature[i], role = outcomes$role[i],
    success_seed1 = seed1[["success"]], success_seed2 = seed2[["success"]],
    y_importance = mean(c(seed1[["importance"]], seed2[["importance"]])), y_rank = mean(c(seed1[["rank"]], seed2[["rank"]]))
  )
}, mc.cores = n_cores))
results[, symbol := gene_symbol(feature)]
results[, survives := success_seed1 >= forest_par$confounder_survival_fraction & success_seed2 >= forest_par$confounder_survival_fraction]
write_result(results[order(role, -success_seed1)], "forest_confounders")
print(results[role == "hit"])
summary_tab <- results[, .(outcomes = .N, surviving = sum(survives), survival_rate = mean(survives), median_success = median((success_seed1 + success_seed2) / 2)), by = role]
summary_tab[, drug_hits_checked := 0L]
write_result(summary_tab, "forest_confounders_summary")
print(summary_tab)
log_msg("drug hits to check: %d (the pan-cancer drug scan produced none at FDR < %.2f)", nrow(drug_hits), params$drug_scan$fdr)
