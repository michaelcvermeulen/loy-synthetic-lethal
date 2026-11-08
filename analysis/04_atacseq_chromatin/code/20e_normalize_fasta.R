script <- sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
fa <- 'config/TF_extension/hg38_primary.fa'; dest <- paste0(fa,'.LF')
stopifnot(file.exists(fa), !file.exists(dest))
input <- file(fa,'rb'); output <- file(dest,'wb')
removed <- 0; total <- 0
repeat {
 chunk <- readBin(input,'raw',n=4*1024*1024)
 if(!length(chunk)) break
 keep <- chunk != as.raw(13); removed <- removed+sum(!keep); total <- total+length(chunk)
 writeBin(chunk[keep],output)
}
close(input);close(output)
stopifnot(file.info(dest)$size==total-removed)
stopifnot(file.rename(fa,'config/TF_extension/hg38_primary_CRLF.fa'))
if(file.exists(paste0(fa,'.fai'))) stopifnot(file.rename(paste0(fa,'.fai'),'config/TF_extension/hg38_primary_CRLF.fa.fai'))
stopifnot(file.rename(dest,fa))
writeLines(paste('Normalized FASTA: removed',removed,'CR bytes; sequence bases preserved'),paste0(fa,'.LF_verified'))
suppressPackageStartupMessages({library(data.table);library(Biostrings);library(BSgenome.Hsapiens.UCSC.hg38)})
p <- fread('differential/primary_TMM_all.csv')
s <- getSeq(BSgenome.Hsapiens.UCSC.hg38,GRanges(p$chrom,IRanges(p$start+1L,p$end)))
fwrite(data.table(peak_id=p$peak_id,sequence=as.character(s)),'config/TF_extension/tested_peak_sequences.tsv',sep='\t')
cat('FASTA converted and independent peak-sequence reference saved\n')
