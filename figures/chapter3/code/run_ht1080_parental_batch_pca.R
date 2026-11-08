args <- commandArgs(trailingOnly=TRUE)
source_dir <- normalizePath(args[1],winslash="/",mustWork=TRUE)
out_dir <- normalizePath(args[2],winslash="/",mustWork=TRUE)
suppressPackageStartupMessages({library(edgeR);library(limma)})
x <- readRDS(file.path(source_dir,"outputs/models/01_ingest.rds"))
e <- readRDS(file.path(source_dir,"outputs/models/02_expression.rds"))
cfg <- yaml::read_yaml(file.path(source_dir,"config/analysis.yml"))
source(file.path(source_dir,"R/02_effects.R"))
meta <- x$meta
meta$treatment <- factor(meta$treatment,levels=c("CTRL","LOY","PARENTAL"))
meta$batch <- factor(meta$batch,levels=c("AUG_2025","MAY_2026"))
stopifnot(nrow(meta)==15,identical(colnames(x$counts),meta$sample_key),sum(e$norm$keep)==25526)
context <- normalize_counts(x$counts,x$annotation$main_background,interaction(meta$batch,meta$treatment,drop=TRUE),cfg)
expr <- context$logCPM[e$norm$keep,]
preserve <- model.matrix(~treatment,meta)
full <- model.matrix(~treatment+batch,meta)
stopifnot(qr(full)$rank==ncol(full))
adjusted <- removeBatchEffect(expr,batch=meta$batch,design=preserve)
stopifnot(all(is.finite(adjusted)))
# Every sample in one period receives the same gene-specific shift.
shift <- expr-adjusted
within_period_error <- max(sapply(levels(meta$batch),function(b) {
  z <- shift[,meta$batch==b,drop=FALSE]
  max(abs(sweep(z,1,z[,1],"-")))
}))
stopifnot(within_period_error < 1e-10)
pre_fit <- lmFit(expr,full)$coefficients
post_fit <- lmFit(adjusted,full)$coefficients
condition_error <- max(abs(pre_fit[,2:3]-post_fit[,2:3]))
remaining_batch <- max(abs(post_fit[,4]))
stopifnot(condition_error < 1e-10, remaining_batch < 1e-10)
pcs <- list(unadjusted=prcomp(t(expr),center=TRUE,scale.=FALSE),adjusted=prcomp(t(adjusted),center=TRUE,scale.=FALSE))
old <- read.csv(file.path(source_dir,"outputs/tables/PCA_context_all15.csv"))
old <- old[match(rownames(pcs$unadjusted$x),old$sample_key),]
for(j in 1:4) if(cor(pcs$unadjusted$x[,j],old[[paste0("PC",j)]])<0) {
 pcs$unadjusted$x[,j] <- -pcs$unadjusted$x[,j];pcs$unadjusted$rotation[,j] <- -pcs$unadjusted$rotation[,j]
}
stopifnot(max(abs(pcs$unadjusted$x[,1:4]-as.matrix(old[,paste0("PC",1:4)])))<1e-6)
variance <- do.call(rbind,lapply(names(pcs),function(name) {
 pc <- pcs[[name]];write.csv(cbind(meta,pc$x[,1:4]),file.path(out_dir,paste0("PCA_",name,".csv")),row.names=FALSE)
 data.frame(view=name,component=paste0("PC",seq_along(pc$sdev)),variance_fraction=pc$sdev^2/sum(pc$sdev^2))
}))
write.csv(variance,file.path(out_dir,"PCA_variance.csv"),row.names=FALSE)
write.csv(cbind(sample_key=meta$sample_key,full),file.path(out_dir,"model_design.csv"),row.names=FALSE)
write.csv(as.data.frame(table(meta$treatment,meta$batch)),file.path(out_dir,"population_by_period.csv"),row.names=FALSE)
write.csv(data.frame(sample_key=meta$sample_key,context$dge$samples,effective_library_reads=context$effective),file.path(out_dir,"normalization_all15.csv"),row.names=FALSE)
saveRDS(list(unadjusted=expr,period_adjusted_visualization_only=adjusted,metadata=meta,PCA=pcs,preserved_design=preserve),file.path(out_dir,"batch_adjustment_visualization_only.rds"))
writeLines(c("Method: limma::removeBatchEffect; additive period term; population preserved.","Parentals August only; shared period shift informed by CTRL/LOY; no interaction model.",paste("Design rank:",qr(full)$rank,"of",ncol(full)),paste("Maximum change in preserved population coefficients:",condition_error),paste("Maximum residual period coefficient:",remaining_batch),paste("Within-period difference preservation error:",within_period_error),capture.output(args(limma::removeBatchEffect)),capture.output(sessionInfo())),file.path(out_dir,"run_validation.txt"))
print(variance[variance$component %in% c("PC1","PC2"),])
