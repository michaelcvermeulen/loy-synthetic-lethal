args <- commandArgs(trailingOnly=TRUE)
out_dir <- normalizePath(args[1],winslash="/",mustWork=TRUE)
x <- readRDS(file.path(out_dir,"batch_adjustment_visualization_only.rds"))
stopifnot(nrow(x$unadjusted)==25526,identical(colnames(x$unadjusted),x$metadata$sample_key))
result <- list();variance <- list()
for (i in 1:2) {
 batch <- c("AUG_2025","MAY_2026")[i]; name <- paste0("batch",i)
 keep <- x$metadata$batch==batch; m <- x$metadata[keep,,drop=FALSE]
 stopifnot(sum(keep)==c(9,6)[i])
 # Retain common all-15 TMM exposure offsets; fit PCA separately within period.
 pc <- prcomp(t(x$unadjusted[,keep,drop=FALSE]),center=TRUE,scale.=FALSE)
 write.csv(cbind(m,pc$x[,1:4]),file.path(out_dir,paste0("PCA_",name,".csv")),row.names=FALSE)
 variance[[name]] <- data.frame(view=name,component=paste0("PC",seq_along(pc$sdev)),variance_fraction=pc$sdev^2/sum(pc$sdev^2))
 result[[name]] <- pc
}
v <- do.call(rbind,variance)
write.csv(v,file.path(out_dir,"PCA_within_batch_variance.csv"),row.names=FALSE)
saveRDS(result,file.path(out_dir,"within_batch_PCA.rds"))
writeLines(c("Batch1: August2025,9 libraries (3 each parental,CTRL,LOY).","Batch2: May2026,6 libraries (3 each CTRL,LOY).","Input: original all-15 TMM logCPM (prior count2), fixed25526 genes; no period adjustment or within-period renormalisation.","PCA fitted separately per period; centred by gene, not variance-scaled.",capture.output(sessionInfo())),file.path(out_dir,"within_batch_validation.txt"))
print(v[v$component %in% c("PC1","PC2"),])
