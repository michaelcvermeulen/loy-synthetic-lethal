script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
tss <- read.delim("qc/CTRL_11_tss_profile.tsv")
fragments <- read.delim("qc/CTRL_11_fragment_lengths.tsv")
flank <- mean(c(head(tss$insertions, 100), tail(tss$insertions, 100)))
tss$enrichment <- tss$insertions/flank
tss$bin <- floor((tss$offset+2000)/20)
tss <- aggregate(cbind(offset, enrichment) ~ bin, tss, mean)
fragments$fraction <- fragments$fragments/sum(fragments$fragments)
png("qc/pilot_profiles.png", width = 1800, height = 750, res = 180)
par(mfrow = c(1,2), mar = c(4.5,4.5,3,1))
plot(tss$offset, tss$enrichment, type = "l", col = "#0072B2", lwd = 1.5,
     xlab = "Distance from TSS (bp)", ylab = "Insertion enrichment over flanks", main = "CTRL_11 TSS profile")
abline(h = 1, col = "grey60", lty = 2)
show <- fragments$length <= 1000
plot(fragments$length[show], fragments$fraction[show], type = "l", col = "#0072B2", lwd = 1.2,
     xlab = "Fragment length (bp)", ylab = "Fraction of eligible fragments", main = "CTRL_11 fragment distribution")
dev.off()
