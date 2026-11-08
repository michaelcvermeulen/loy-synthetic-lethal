# Read-only audit of original annotation versus final QDNAseq feature metadata.
# No input RDS, count table or figure is changed.
args <- commandArgs(trailingOnly=TRUE)
projects <- normalizePath(args[[1]],winslash="/",mustWork=TRUE)
out <- normalizePath(args[[2]],winslash="/",mustWork=TRUE)
b <- readRDS(file.path(projects,"WGS/t2t_bins_100kb.rds"))
d <- b@data
d$bin_id <- rownames(d)
write.csv(d,file.path(out,"original_bin_annotation.csv"),row.names=FALSE)
for (cell in c("HT1080","A549","LNCaP")) {
  cn <- readRDS(file.path(projects,"WGS",paste0(cell,"_QDNAseq_cn.rds")))
  f <- slot(slot(cn,"featureData"),"data")
  f$bin_id <- rownames(f)
  write.csv(f,file.path(out,paste0(cell,"_final_feature_metadata.csv")),row.names=FALSE)
  cat(cell,"original flags",sum(d$use),"final flags",sum(f$use),"total",nrow(f),"\n")
}
