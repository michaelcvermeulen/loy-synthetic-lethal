options(warn=1)
.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
root <- 'A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS'
out <- file.path(root,'integration/RNA_interaction_CN')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkgs <- c('DESeq2','data.table','jsonlite','digest','QDNAseq','DNAcopy','Rsamtools','BSgenome.Hsapiens.UCSC.hg38','fgsea','limma')
for(p in pkgs) cat(p,if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else'NOT INSTALLED','\n')
suppressPackageStartupMessages(library(DESeq2))
writeLines(capture.output(sessionInfo()),file.path(out,'preflight_session.txt'))
writeLines(capture.output(args(DESeq2::DESeq)),file.path(out,'DESeq_API.txt'))
if(requireNamespace('QDNAseq',quietly=TRUE)) {
  writeLines(capture.output(QDNAseq::getBinAnnotations),file.path(out,'QDNAseq_getBinAnnotations.txt'))
  cat('QDNAseq annotation libraries:',grep('QDNAseq',rownames(installed.packages()),value=TRUE),'\n')
}
