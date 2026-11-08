# Compact summary of the expanded screen GSEA (read-only on saved results);
# writes the robust-set and theme tables beside the source data.
library(data.table)
d <- "A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter4/source_data/screen_expanded_gsea_20261021"
r <- fread(file.path(d, "expanded_fgsea_all.csv"))
options(width = 200)
s <- function(x) signif(x, 2)
short <- function(p) substr(p, 1, 58)

ht <- r[grepl("^HT1080", ranking)]
rob <- ht[, .(n_sig = sum(FDR_global <= 0.05), n = .N,
              same_sign = length(unique(sign(NES))) == 1), by = pathway]
rob <- rob[n == 4 & n_sig == 4 & same_sign, pathway]
robust <- r[ranking == "HT1080 drugZ" & pathway %in% rob,
            .(collection, pathway, NES, FDR_global, leadingEdge)][order(NES)]
robust <- merge(robust, r[ranking == "A549 drugZ", .(pathway, A549_NES = NES, A549_FDR = FDR_global)],
                by = "pathway", all.x = TRUE)[order(NES)]
fwrite(robust, file.path(d, "robust_in_all_four_ht1080_rankings.csv"))
cat("Robust in all four HT1080 rankings:", nrow(robust), " (sensitizer",
    sum(robust$NES > 0), ", resistor", sum(robust$NES < 0), ")\n")
print(robust[, .(collection, set = short(pathway), NES = s(NES), q = s(FDR_global),
                 A549 = s(A549_NES))], nrows = 200)

cat("\nA549 drugZ sets at FDR <= 0.05:\n")
print(r[ranking == "A549 drugZ" & FDR_global <= 0.05,
        .(collection, set = short(pathway), NES = s(NES), q = s(FDR_global))])

theme <- r[grepl("SAGA|STAGA|CHOLESTEROL_BIOSYN|STEROL_BIOSYN|SREBF|SREBP|MYC_TARGETS|MEDIATOR_COMPLEX|HISTONE_ACETYLTRANSFERASE_COMPLEX|^GOCC_ATAC|LIPID_HOMEOSTASIS", pathway)]
tw <- dcast(theme, pathway ~ ranking, value.var = c("NES", "FDR_global"))
keep <- c("pathway", "NES_HT1080 drugZ", "FDR_global_HT1080 drugZ",
          "NES_HT1080 IsogenicZ", "FDR_global_HT1080 IsogenicZ",
          "NES_A549 drugZ", "FDR_global_A549 drugZ")
tw <- tw[, ..keep]
fwrite(tw, file.path(d, "selected_theme_sets.csv"))
setnames(tw, c("set", "dz_NES", "dz_q", "iso_NES", "iso_q", "A549_NES", "A549_q"))
cat("\nSelected themes:\n")
print(tw[order(dz_q), lapply(.SD, function(x) if (is.numeric(x)) s(x) else short(x))], nrows = 200)
