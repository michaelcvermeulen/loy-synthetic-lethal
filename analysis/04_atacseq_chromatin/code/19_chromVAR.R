options(warn=1)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash='/', mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table); library(chromVAR); library(motifmatchr)
 library(TFBSTools); library(BSgenome.Hsapiens.UCSC.hg38); library(SummarizedExperiment)
 library(BiocParallel); library(limma); library(ggplot2); library(pheatmap); library(jsonlite)})
register(SerialParam()); set.seed(20260908)
out <- 'motifs/TF_extension/chromVAR'; dir.create(out, recursive=TRUE, showWarnings=FALSE)
cfg <- 'config/TF_extension'
samples <- fread('config/samples.tsv', data.table=FALSE)
samples$condition <- factor(samples$condition, levels=c('CTRL','LOY')); samples$passage <- factor(samples$passage)
res <- fread('differential/primary_TMM_all.csv')
tab <- fread('counts/peak_counts.tsv'); stopifnot(!anyDuplicated(tab$peak_id))
sample_names <- samples$sample
counts <- as.matrix(tab[match(res$peak_id, tab$peak_id), ..sample_names])
rownames(counts) <- res$peak_id
stopifnot(identical(colnames(counts), sample_names), !anyNA(counts), all(counts >= 0), all(counts == round(counts)))
genome <- BSgenome.Hsapiens.UCSC.hg38
gr <- GRanges(res$chrom, IRanges(res$start+1L, res$end)); names(gr) <- res$peak_id
seqs <- getSeq(genome, gr)
fwrite(data.table(peak_id=res$peak_id,sequence=as.character(seqs)),file.path(cfg,'tested_peak_sequences.tsv'),sep='\t')
clean <- as.numeric(letterFrequency(seqs, 'N')) == 0
fwrite(data.table(peak_id=res$peak_id, included=clean), file.path(out,'peak_accounting.csv'))
res <- res[clean]; counts <- counts[clean,,drop=FALSE]; gr <- gr[clean]; seqs <- seqs[clean]
se <- SummarizedExperiment(assays=list(counts=counts), rowRanges=gr, colData=DataFrame(samples))
se <- addGCBias(se, genome=genome)
pfm <- readJASPARMatrix(file.path(cfg,'JASPAR2024_CORE_vertebrates.jaspar'), matrixClass='PFM')
ids <- vapply(pfm, ID, character(1)); names(pfm) <- ids
motif_names <- vapply(pfm, name, character(1))
meta <- fread(file.path(cfg,'Vierstra_v2.1beta_metadata.tsv'))
mapping <- rbindlist(lapply(seq_along(ids), function(i) {
 exact <- unique(meta[motif_id == ids[i], cluster])
 base <- sub('\\..*$', '', ids[i])
 hits <- if(length(exact)) exact else unique(meta[sub('\\..*$', '', motif_id) == base, cluster])
 status <- if(length(exact)==1) 'exact_version' else if(length(hits)==1) 'same_accession_other_version' else if(length(hits)>1) 'ambiguous' else 'unmapped'
 data.table(motif_id=ids[i], motif_name=motif_names[i], archetype=if(length(hits)==1) hits else paste0('unmapped_',ids[i]), mapping_status=status)
}))
fwrite(mapping, file.path(out,'motif_archetype_mapping.csv'))
saveRDS(pfm, file.path(cfg,'JASPAR2024_CORE_vertebrates.rds'))
cat('Scanning',length(pfm),'motifs in',nrow(se),'peaks\n')
motifs <- matchMotifs(pfm, seqs, out='matches', p.cutoff=5e-5)
mm <- motifMatches(motifs); colnames(mm) <- ids
saveRDS(mm, file.path(out,'motif_matches.rds'))
# Collapsed archetypes are unions of member motif-positive peaks, not minimum member P-values.
groups <- split(seq_len(ncol(mm)), mapping$archetype)
family_mm <- do.call(cbind,lapply(groups, function(ix) Matrix::rowSums(mm[,ix,drop=FALSE]) > 0))
family_mm <- Matrix::Matrix(family_mm, sparse=TRUE); colnames(family_mm) <- names(groups)
fwrite(rbindlist(lapply(names(groups),function(k) data.table(archetype=k,
 member_motifs=paste(ids[groups[[k]]],collapse=';'), member_names=paste(unique(motif_names[groups[[k]]]),collapse=';')))),file.path(out,'archetype_members.csv'))
cat('Sampling GC/accessibility matched backgrounds\n')
bg <- getBackgroundPeaks(se, niterations=50)
saveRDS(bg,file.path(out,'background_peaks.rds'))
all_summaries <- list()
for (kind in c('motifs','archetypes')) {
 mat <- if(kind=='motifs') mm else family_mm
 dev <- computeDeviations(se, annotations=mat, background_peaks=bg)
 saveRDS(dev, file.path(out,paste0(kind,'_deviations.rds')))
 z <- deviationScores(dev); d <- deviations(dev)
 fwrite(data.table(id=rownames(z),z),file.path(out,paste0(kind,'_zscores.csv')))
 fwrite(data.table(id=rownames(d),d),file.path(out,paste0(kind,'_bias_corrected_deviations.csv')))
 valid <- rowSums(!is.finite(z))==0 & apply(z,1,sd)>0
 z <- z[valid,,drop=FALSE]
 for (arm in c('independent','passage_sensitivity')) {
  design <- if(arm=='independent') model.matrix(~condition,samples) else model.matrix(~passage+condition,samples)
  stopifnot(qr(design)$rank==ncol(design))
  fit <- eBayes(lmFit(z,design))
  stats <- as.data.table(topTable(fit,coef='conditionLOY',number=Inf,sort.by='P'),keep.rownames='id')
  setnames(stats,'logFC','LOY_minus_CTRL_zscore')
  if(kind=='motifs') stats <- merge(stats,mapping,by.x='id',by.y='motif_id',all.x=TRUE,sort=FALSE)
  setorder(stats,adj.P.Val,P.Value)
  fwrite(stats,file.path(out,paste0(kind,'_',arm,'_limma.csv')))
  all_summaries[[paste(kind,arm)]] <- data.table(kind,arm,tested=nrow(stats),FDR005=sum(stats$adj.P.Val<.05),FDR010=sum(stats$adj.P.Val<.1))
  top <- head(stats$id,25)
  labels <- if(kind=='motifs') paste(mapping$motif_name[match(top,mapping$motif_id)],top) else vapply(top,function(k) {
   members <- unique(motif_names[groups[[k]]]); paste(k,paste(head(members,3),collapse='/'),if(length(members)>3)'...'else'')
  },character(1))
  heat <- z[top,,drop=FALSE]; rownames(heat) <- labels
  ann <- samples[,c('condition','passage')]; rownames(ann) <- sample_names
  for(ext in c('png','pdf')) pheatmap(heat,scale='none',cluster_cols=FALSE,cluster_rows=FALSE,
   annotation_col=ann,main=paste('chromVAR:',kind,arm,'(top ranked; inspect FDR)'),fontsize_row=8,
   filename=file.path(out,paste0(kind,'_',arm,'_heatmap.',ext)),width=12,height=9)
 }
 # Requested simple group test is included as a sensitivity result with its own BH correction.
 welch <- rbindlist(lapply(seq_len(nrow(z)),function(i){
  x <- z[i,samples$condition=='LOY']; y <- z[i,samples$condition=='CTRL']
  p <- tryCatch(t.test(x,y)$p.value,error=function(e) NA_real_)
  data.table(id=rownames(z)[i],LOY_minus_CTRL_zscore=mean(x)-mean(y),PValue=p)
 }))
 welch[,FDR:=p.adjust(PValue,'BH')]; fwrite(welch,file.path(out,paste0(kind,'_Welch_test.csv')))
}
summary <- rbindlist(all_summaries); fwrite(summary,file.path(out,'analysis_summary.csv')); print(summary)
writeLines(c('chromVAR 1.32.0 on unnormalized integer fragment counts over tested non-Y/PAR peaks without N bases.',
 'JASPAR2024 CORE vertebrate nonredundant motifs; motifmatchr p.cutoff=5e-5. Seed 20260908; 50 background iterations.',
 'Scores are chromVAR background-standardized deviations, not z-scores standardized across the six samples.',
 'Archetypes use the union of member motif-positive peaks. Exact JASPAR ID mapping is preferred; unique accession mappings across versions are labeled explicitly. Ambiguous/unmapped motifs remain separate.',
 'Limma compares sample-level z-scores, with BH within each motif/archetype analysis arm. Effect is a z-score difference, not log2 fold change.',
 'Independent-replicate and passage-adjusted models are reported separately. Technical batch metadata remain unconfirmed.',
 'Welch t-tests are additional sensitivity results. Heatmaps rank by limma FDR; displayed rows are not necessarily significant.',
 'Motif accessibility and RNA expression directions need not agree; these statistics do not prove TF binding.'),file.path(out,'METHODS.txt'))
writeLines(capture.output(sessionInfo()),'logs/chromVAR_session.txt')
writeLines('PASS: sample and peak alignment, integer counts, full-rank designs, reproducible motif/family mappings',file.path(out,'verification.txt'))
# Export identical reference sequence for external motif scanning/footprinting without another genome download.
fa <- file.path(cfg,'hg38_primary.fa')
if(!file.exists(paste0(fa,'.complete'))) {
 if(file.exists(fa) || file.exists(paste0(fa,'.export'))) stop('Incomplete FASTA export exists; inspect before resuming')
 for(chr in c(paste0('chr',1:22),'chrX','chrY')) {
  s <- DNAStringSet(genome[[chr]]); names(s) <- chr
  writeXStringSet(s,paste0(fa,'.export'),append=file.exists(paste0(fa,'.export')),compress=FALSE)
 }
 # Windows text-mode export contains CRLF; normalize bytes for external genomic tools.
 input <- file(paste0(fa,'.export'),'rb'); output <- file(fa,'wb')
 repeat { chunk <- readBin(input,'raw',n=4*1024*1024); if(!length(chunk))break; writeBin(chunk[chunk!=as.raw(13)],output) }
 close(input);close(output)
 writeLines('FASTA line endings normalized to LF',paste0(fa,'.LF_verified'))
 writeLines('BSgenome.Hsapiens.UCSC.hg38 primary chromosomes',paste0(fa,'.complete'))
}
