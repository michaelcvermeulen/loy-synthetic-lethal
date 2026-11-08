# Build measured peak features without reading RNA effect values.
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(jsonlite)
})
set.seed(42)
args <- commandArgs(trailingOnly = FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)][1])
root <- dirname(normalizePath(script, winslash = '/'))
out <- file.path(root, 'outputs/tables/chromatin_atlas')
prov <- file.path(root, 'provenance/chromatin_atlas')
config <- fromJSON(file.path(root, 'config/chromatin_atlas.json'))
marks <- fread(file.path(out, 'mark_inputs.csv'))
genes <- fread(file.path(out, 'genes_with_primary_TSS.csv'))
tss <- fread(file.path(out, 'distinct_TSS.csv.gz'))
sizes <- fread(file.path(out, 'reference_chromosomes.csv'))
headers <- fromJSON(file.path(prov, 'bigwig_headers.json'), simplifyVector = FALSE)
size_map <- setNames(sizes$length, sizes$chrom)

as_gr <- function(d) GRanges(d$chrom, IRanges(d$start0 + 1L, d$end0))
as_bed <- function(gr) data.table(chrom = as.character(seqnames(gr)),
                                  start0 = start(gr) - 1L, end0 = end(gr))
bl <- fread(config$blacklist, header = FALSE, select = 1:3)
setnames(bl, c('chrom', 'start0', 'end0'))
blacklist <- reduce(as_gr(bl))

overlap_summary <- function(regions, peaks) {
  hit <- findOverlaps(regions, peaks, ignore.strand = TRUE)
  bp <- numeric(length(regions)); n <- integer(length(regions))
  if (length(hit)) {
    q <- queryHits(hit); s <- subjectHits(hit)
    w <- pmin(end(regions)[q], end(peaks)[s]) - pmax(start(regions)[q], start(peaks)[s]) + 1
    totals <- data.table(q = q, w = w)[, .(bp = sum(w), n = .N), by = q]
    bp[totals$q] <- totals$bp; n[totals$q] <- totals$n
  }
  list(bp = bp, n = n)
}

prepare_regions <- function(d) {
  d <- copy(d)
  d[, supported_chromosome := chrom %in% names(size_map)]
  # Non-reference contigs remain in the gene table; they are not assigned false-negative marks.
  d <- d[supported_chromosome == TRUE]
  d[, start0 := pmax(0L, start0)]
  d[, end0 := pmin(end0, size_map[chrom])]
  stopifnot(all(d$end0 > d$start0))
  gr <- as_gr(d)
  d[, total_bp := end0 - start0]
  d[, unmasked_bp := total_bp - overlap_summary(gr, blacklist)$bp]
  d[, unmasked_fraction := unmasked_bp / total_bp]
  d[, usable := unmasked_fraction >= config$minimum_unmasked_fraction & unmasked_bp > 0]
  d
}

all_tss <- copy(tss)
all_tss[, `:=`(start0 = tss0-config$primary_promoter_half_width,
               end0 = tss0+config$primary_promoter_half_width+1L)]
all_tss <- prepare_regions(all_tss)
all_tss_gr <- as_gr(all_tss)
windows <- c(config$primary_promoter_half_width, config$sensitivity_promoter_half_widths)
prom <- rbindlist(lapply(windows, function(w) {
  d <- genes[, .(gene_id, gene_name, chrom, strand, transcript_id, tss0, in_RNA_rank_universe)]
  d[, `:=`(window = w, start0 = tss0-w, end0 = tss0+w+1L)]
  prepare_regions(d)
}))
prom_gr <- as_gr(prom)
bodies <- prepare_regions(genes[, .(gene_id, gene_name, chrom, strand,
                                     in_RNA_rank_universe, start0 = start-1L, end0 = end)])
body_gr <- as_gr(bodies)
stopifnot(all(genes[in_RNA_rank_universe == TRUE, gene_id] %in% bodies$gene_id))
peak_rows <- list(); main_rows <- list(); body_rows <- list(); qc <- list()
all_tss_features <- copy(all_tss)
region_dir <- file.path(out, 'regions')
dir.create(region_dir, showWarnings = FALSE)

for (i in seq_len(nrow(marks))) {
  m <- marks[i]
  message('Annotating ', m$mark)
  raw <- fread(m$peak_path, header = FALSE)
  stopifnot(ncol(raw) %in% c(9L, 10L), all(is.finite(raw$V9)),
            all(raw$V2 >= 0), all(raw$V3 > raw$V2))
  names(raw)[1:3] <- c('chrom', 'start0', 'end0')
  called <- raw[chrom %in% names(size_map) & V9 >= m$min_minuslog10q]
  relaxed <- raw[chrom %in% names(size_map) & V9 >= m$sensitivity_min_minuslog10q]
  stopifnot(all(called$end0 <= size_map[called$chrom]))
  unmasked <- setdiff(reduce(as_gr(called)), blacklist, ignore.strand = TRUE)
  unmasked_relaxed <- setdiff(reduce(as_gr(relaxed)), blacklist, ignore.strand = TRUE)
  fwrite(as_bed(unmasked), file.path(region_dir, paste0(m$mark, '_filtered_unmasked.bed.gz')),
         sep = '\t', col.names = FALSE, compress = 'gzip')
  seq_available <- names(headers[[m$mark]]$chromosomes)
  summarize <- function(d, gr) {
    z <- overlap_summary(gr, unmasked)
    rz <- overlap_summary(gr, unmasked_relaxed)
    usable <- d$usable & d$chrom %in% seq_available
    data.table(peak_covered_bp = z$bp,
               peak_covered_fraction = ifelse(usable, z$bp/d$unmasked_bp, NA_real_),
               n_merged_peak_segments = z$n,
               detected = ifelse(usable, z$bp > 0, NA),
               detected_relaxed_q = ifelse(usable, rz$bp > 0, NA),
               evaluable = usable)
  }
  tp <- summarize(all_tss, all_tss_gr)
  for (col in c('detected','detected_relaxed_q','peak_covered_fraction','evaluable')) {
    set(all_tss_features, j = paste(m$mark,col,sep = '__'), value = tp[[col]])
  }
  mp <- cbind(copy(prom), summarize(prom, prom_gr)); mp[, mark := m$mark]
  main_rows[[i]] <- mp
  bp <- cbind(copy(bodies), summarize(bodies, body_gr)); bp[, mark := m$mark]
  body_rows[[i]] <- bp
  qc[[i]] <- data.table(mark = m$mark, raw_peaks = nrow(raw),
                       q_threshold = m$min_minuslog10q, retained_peaks = nrow(called),
                       relaxed_q_threshold = m$sensitivity_min_minuslog10q,
                       retained_relaxed_peaks = nrow(relaxed),
                       median_peak_width = median(called$end0-called$start0),
                       merged_masked_segments = length(unmasked),
                       unmasked_peak_bp = sum(width(unmasked)),
                       bp_removed_by_blacklist = sum(width(reduce(as_gr(called))))-sum(width(unmasked)))
}
main <- rbindlist(main_rows)
body <- rbindlist(body_rows)
fwrite(main, file.path(out,'promoter_peak_features_long.csv.gz'), compress = 'gzip')
fwrite(body, file.path(out,'gene_body_peak_features_long.csv.gz'), compress = 'gzip')
fwrite(all_tss_features, file.path(out,'all_TSS_peak_features.csv.gz'), compress = 'gzip')
fwrite(rbindlist(qc), file.path(out,'peak_input_qc.csv'))

# Export exact unmasked pieces for the separate Linux bigWig reader.
signal <- rbindlist(list(
  prom[window == config$primary_promoter_half_width,
       .(gene_id,chrom,start0,end0,unmasked_bp,region_type='promoter')],
  bodies[, .(gene_id,chrom,start0,end0,unmasked_bp,region_type='gene_body')]))
signal[, region_id := paste(gene_id,region_type,sep=':')]
stopifnot(!anyDuplicated(signal$region_id))
sg <- as_gr(signal)
pieces <- setdiff(reduce(sg),blacklist,ignore.strand=TRUE)
h <- findOverlaps(sg,pieces,ignore.strand=TRUE)
q <- queryHits(h); s <- subjectHits(h)
piece_table <- data.table(region_id=signal$region_id[q],chrom=signal$chrom[q],
                         start0=pmax(signal$start0[q],start(pieces)[s]-1L),
                         end0=pmin(signal$end0[q],end(pieces)[s]))
check <- piece_table[,.(observed_bp=sum(end0-start0)),by=region_id]
stopifnot(all(check$observed_bp==signal$unmasked_bp[match(check$region_id,signal$region_id)]))
fwrite(signal,file.path(out,'signal_regions.csv'))
fwrite(piece_table,file.path(out,'signal_region_unmasked_pieces.csv.gz'),compress='gzip')

# Candidate distal regulatory regions: overlapping K4me1/K27ac peaks outside all annotated TSS windows.
k4 <- fread(file.path(region_dir,'H3K4me1_filtered_unmasked.bed.gz'),header=FALSE)
kac <- fread(file.path(region_dir,'H3K27ac_filtered_unmasked.bed.gz'),header=FALSE)
setnames(k4,c('chrom','start0','end0')); setnames(kac,c('chrom','start0','end0'))
candidate <- intersect(as_gr(k4),as_gr(kac),ignore.strand=TRUE)
distal <- setdiff(candidate,reduce(all_tss_gr),ignore.strand=TRUE)
distal <- distal[width(distal)>=100L]
fwrite(as_bed(distal),file.path(region_dir,'candidate_distal_H3K4me1_H3K27ac.bed.gz'),
       sep='\t',col.names=FALSE,compress='gzip')
write_json(list(primary_window=config$primary_promoter_half_width, promoter_rows=nrow(main),
                gene_body_rows=nrow(body), annotated_TSS=nrow(all_tss_features),
                candidate_distal_intervals=length(distal),
                distal_gene_links='Not assigned; no evidence-based links supplied'),
           file.path(prov,'peak_feature_summary.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(prov,'R_sessionInfo.txt'))
message('Peak features complete.')
