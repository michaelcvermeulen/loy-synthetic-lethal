options(warn=1,stringsAsFactors=FALSE)
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
suppressPackageStartupMessages({library(limma);library(edgeR);library(data.table)})
source("R/00_common.R")
setDTthreads(2)
paths <- c("outputs/models/01_ingest.rds","outputs/models/02_expression.rds",
           "outputs/tables/batchDE_batch_condition__combined.csv")
before <- vapply(paths,hash_file,character(1))
x <- readRDS(paths[1]);e <- readRDS(paths[2]);m <- x$core
# Companion analysis: retain existing expression filter, add expressed PAR loci.
a <- e$effects[e$effects$primary_eligible | (e$effects$PAR & e$effects$expressed),]
counts <- x$counts[a$gene_id,m$sample_key,drop=FALSE]
effective <- setNames(e$norm$effective,colnames(e$logCPM))[m$sample_key]
m$treatment <- factor(m$treatment,levels=c("CTRL","LOY"))
m$batch <- factor(m$batch,levels=c("AUG_2025","MAY_2026"))
D <- model.matrix(~batch+treatment,m);rownames(D) <- m$sample_key
stopifnot(qr(D)$rank==3,nrow(D)==12,sum(a$primary_eligible)==25526,
          all(counts==round(counts)),all(is.finite(effective)))
v <- voom(counts,design=D,lib.size=effective,normalize.method="none",plot=FALSE,save.plot=TRUE)
fit <- eBayes(lmFit(v,D),trend=FALSE,robust=TRUE)
z <- topTable(fit,coef="treatmentLOY",number=Inf,sort.by="none",adjust.method="BH",confint=TRUE)
z <- z[match(a$gene_id,rownames(z)),]
names(z)[names(z)=="P.Value"] <- "PValue"
names(z)[names(z)=="adj.P.Val"] <- "FDR"
z <- cbind(a[,c("gene_id","gene_name","gene_type","chrom","start","end","PAR","primary_eligible","mean_CPM")],z)
z$passes_cutoffs <- z$FDR<.05 & abs(z$logFC)>=1
stopifnot(all(is.finite(z$FDR)),max(abs(z$FDR-p.adjust(z$PValue,"BH")))<1e-12)
tab(z,"batchDE_batch_condition_with_PAR__combined")
tab(data.frame(sample_key=m$sample_key,D,effective_library_reads=effective,check.names=FALSE),"batchDE_with_PAR_design")
saveRDS(list(result=z,voom=v,fit=fit,design=D),"outputs/models/batch_condition_with_PAR.rds")
summary <- data.frame(tested=nrow(z),PAR_tested=sum(z$PAR),
 higher=sum(z$passes_cutoffs & z$logFC>0),lower=sum(z$passes_cutoffs & z$logFC<0),
 PAR_pass=sum(z$passes_cutoffs & z$PAR))
tab(summary,"batchDE_with_PAR_summary")
stopifnot(identical(before,vapply(paths,hash_file,character(1))))
fwrite(data.frame(path=paths,sha256=before),"provenance/batchDE_with_PAR_input_hashes.csv")
writeLines(c("Companion model includes expressed PAR features; original model and inputs unchanged.",
 "Same 12 core libraries, original TMM exposures, batch + condition design. BH across all included features.",
 "PAR membership inherited from pinned GRCh38 interval-overlap annotation; no inference of X/Y allele origin.",
 capture.output(summary),capture.output(sessionInfo())),"provenance/batchDE_with_PAR_validation.txt")
print(summary)
print(z[z$PAR,c("gene_name","chrom","logFC","FDR","passes_cutoffs")])
