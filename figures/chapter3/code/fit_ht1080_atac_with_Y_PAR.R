suppressPackageStartupMessages({library(edgeR);library(data.table)})
run <- "A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS"
out <- "A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/source_data/ht1080_atac_chrY"
s <- fread(file.path(run,"config/samples.tsv"),data.table=FALSE)
s$condition <- factor(s$condition,levels=c("CTRL","LOY"))
p <- fread(file.path(run,"peaks/consensus_membership.tsv"),data.table=FALSE)
c <- fread(file.path(run,"counts/peak_counts.tsv"),data.table=FALSE)
n <- fread(file.path(run,"differential/normalization_factors.csv"),data.table=FALSE)
n <- n[n$arm=="primary_TMM",];n <- n[match(s$sample,n$sample),]
stopifnot(identical(p$peak_id,c$peak_id))
m <- as.matrix(c[,s$sample]);rownames(m)<-c$peak_id
design <- model.matrix(~condition,s)
y <- DGEList(m,lib.size=n$library_size)
keep <- filterByExpr(y,design=design)
# Preserve the exact primary tested universe and filter the additional regions by the same rule.
primary <- fread(file.path(run,"differential/primary_TMM_all.csv"),data.table=FALSE)
extra <- p$chrom=="chrY" | p$is_PAR %in% c(TRUE,"True","TRUE")
keep <- p$peak_id %in% primary$peak_id | (extra & keep)
y <- y[keep,,keep.lib.sizes=TRUE];y$samples$norm.factors <- n$factor
y <- estimateDisp(y,design,robust=TRUE)
fit <- glmQLFit(y,design,robust=TRUE)
test <- glmQLFTest(fit,coef=2)
t <- topTags(test,n=Inf,sort.by="none")$table
res <- cbind(p[keep,],t)
res$region <- ifelse(res$is_PAR %in% c(TRUE,"True","TRUE"),paste0(res$chrom," PAR"),ifelse(res$chrom=="chrY","chrY non-PAR","Other"))
write.csv(res,file.path(out,"DAR_expanded_Y_PAR.csv"),row.names=FALSE)
write.csv(data.frame(sample=s$sample,effective_size=y$samples$lib.size*y$samples$norm.factors),file.path(out,"DAR_normalization.csv"),row.names=FALSE)
saveRDS(list(y=y,fit=fit,design=design),file.path(out,"DAR_expanded_model.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"DAR_session.txt"))
print(table(res$region,res$FDR<.05));print(range(res$logFC))
