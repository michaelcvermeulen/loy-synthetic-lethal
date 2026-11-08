# Re-tune the "use" flag on the existing T2T bin annotation
# without rebuilding GC/mappability/blacklist (which take time).
suppressPackageStartupMessages(library(QDNAseq))

bins <- readRDS("../t2t_bins_100kb.rds")
df   <- bins@data

# Slight tightening — keeps chrY MSY analyzable while excluding worst bins.
# bases >= 60%: mostly real sequence (vs 50% previously)
# mappability >= 55%: well-mappable 50-mers (vs 50% previously; gentle on chrY)
# blacklist < 20%: low difficulty-region overlap (vs 25% previously)
df$use <- df$bases >= 60 & df$mappability >= 55 & df$blacklist < 20
bins@data <- df

n_use <- sum(df$use)
cat("Bins passing use filter:", n_use, "/", nrow(df),
    "(", round(100 * n_use / nrow(df), 1), "%)\n")
cat("By chromosome:\n")
print(tapply(df$use, df$chromosome, function(x) paste0(sum(x), "/", length(x))))

saveRDS(bins, "../t2t_bins_100kb.rds")
cat("Updated ../t2t_bins_100kb.rds\n")
