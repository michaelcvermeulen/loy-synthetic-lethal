read_gtf_genes <- function(path) {
  con <- gzfile(path, open = 'rt'); on.exit(close(con))
  gene_lines <- list(); k <- 0L
  repeat {
    lines <- readLines(con, n = 100000, warn = FALSE)
    if (!length(lines)) break
    hit <- grepl('^[^\t]+\t[^\t]+\tgene\t', lines, perl = TRUE)
    if (any(hit)) { k <- k + 1L; gene_lines[[k]] <- lines[hit] }
  }
  fields <- strsplit(unlist(gene_lines, use.names = FALSE), '\t', fixed = TRUE)
  assert(all(lengths(fields) == 9), 'Malformed GTF gene record')
  mat <- matrix(unlist(fields, use.names = FALSE), ncol = 9, byrow = TRUE)
  attribute <- function(key) {
    pattern <- paste0('(?:^|;[[:space:]]*)', key, ' "([^"]*)"')
    hits <- regmatches(mat[, 9], regexec(pattern, mat[, 9], perl = TRUE))
    vapply(hits, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  }
  a <- data.frame(gene_id = attribute('gene_id'), gene_name = attribute('gene_name'),
    gene_type = attribute('gene_type'), chrom = mat[, 1], start = as.integer(mat[, 4]),
    end = as.integer(mat[, 5]), strand = mat[, 7])
  assert(!anyNA(a$gene_id) && !anyNA(a$gene_type), 'Missing required gene attributes')
  a
}
ingest <- function(cfg, repo) {
  meta <- as.data.frame(data.table::fread(cfg$sample_manifest))
  raw_path <- file.path(repo, cfg$raw_counts)
  raw <- data.table::fread(raw_path)
  ids <- raw[[1]]
  assert(!anyDuplicated(ids), 'Duplicate exact counting gene IDs')
  counts <- as.matrix(raw[, -1, with = FALSE])
  rownames(counts) <- ids
  counts <- align_samples(counts, meta)
  assert(all(is.finite(counts) & counts >= 0 & counts == floor(counts)),
         'Counts must be nonnegative finite integers')
  storage.mode(counts) <- 'integer'
  audit <- list(); summaries <- list(); feature_metadata <- NULL
  for (i in seq_len(nrow(meta))) {
    path <- meta$source_count_file[i]
    f <- data.table::fread(path, skip = 1)
    assert(ncol(f) == 7, paste('Expected one counted library:', path))
    assert(setequal(f$Geneid, ids) && !anyDuplicated(f$Geneid), paste('Feature mismatch:', path))
    f <- f[match(ids, f$Geneid)]
    assert(all(f[[7]] == counts[, i]), paste('Source counts differ:', path))
    fm <- as.data.frame(f[, 1:6])
    if (is.null(feature_metadata)) feature_metadata <- fm else
      assert(identical(fm, feature_metadata), paste('Counting feature metadata differ:', path))
    h <- readLines(path, n = 1, warn = FALSE)
    assert(grepl('v2.1.1', h, fixed = TRUE) && grepl('gencode.v49', h, fixed = TRUE),
           paste('Unexpected counting version/annotation:', path))
    flags <- gsub('"', '', h, fixed = TRUE)
    assert(grepl(' -p ', flags, fixed = TRUE) && !grepl('--countReadPairs', flags, fixed = TRUE) &&
             !grepl(' -s ', flags, fixed = TRUE), paste('Counting options changed:', path))
    audit[[i]] <- data.frame(sample_key = meta$sample_key[i], source = path,
                            sha256 = hash_file(path), source_matches_matrix = TRUE,
                            count_unit = 'reads', strandedness = 'unstranded_default', command = h)
    s <- data.table::fread(paste0(path, '.summary'))
    names(s) <- c('status', 'reads'); s$sample_key <- meta$sample_key[i]
    assert(s$reads[s$status == 'Assigned'] == sum(counts[, i]), 'Assigned summary/count mismatch')
    summaries[[i]] <- as.data.frame(s)
  }
  message('All 6 source tables match; importing exact GENCODE v49 gene records')
  a <- read_gtf_genes(file.path(repo, cfg$gtf))
  assert(!anyDuplicated(a$gene_id), 'Duplicate GTF gene entries')
  assert(all(ids %in% a$gene_id), 'Unmapped exact counting IDs in GTF')
  a <- a[match(ids, a$gene_id), ]; rownames(a) <- NULL
  a$stable_gene_id <- stable_id(a$gene_id)
  a$feature_length_bp <- feature_metadata$Length
  assert(all(a$feature_length_bp > 0), 'Nonpositive counting feature length')
  fc_chr <- sub(';.*', '', feature_metadata$Chr)
  assert(all(a$chrom == fc_chr), 'GTF/count chromosome mismatch')
  par <- data.table::fread('config/hg38_PAR.bed', header = FALSE)
  a$PAR <- FALSE
  for (i in seq_len(nrow(par))) {
    a$PAR <- a$PAR | (a$chrom == par[[1]][i] & a$start <= par[[3]][i] & a$end > par[[2]][i])
  }
  a$main_background <- a$chrom %in% c(paste0('chr', 1:22), 'chrX') & !a$PAR
  a$ribosomal_protein <- grepl('^RP[SL][0-9]', a$gene_name)
  a$region <- ifelse(a$PAR, paste0(a$chrom, '_PAR'), a$chrom)
  cyto <- data.table::fread('resources/hg38_cytoBand.txt')
  a$arm <- NA_character_
  for (ch in unique(a$chrom)) {
    cb <- cyto[chrom == ch]
    if (!nrow(cb) || !any(startsWith(cb$band, 'p')) || !any(startsWith(cb$band, 'q'))) next
    p_end <- max(cb$end[startsWith(cb$band, 'p')], na.rm = TRUE)
    q_start <- min(cb$start[startsWith(cb$band, 'q')], na.rm = TRUE)
    ix <- which(a$chrom == ch)
    a$arm[ix] <- ifelse(a$end[ix] <= p_end, paste0(ch, 'p'),
                       ifelse(a$start[ix] > q_start, paste0(ch, 'q'), 'centromere_spanning'))
  }
  iy <- a$chrom == 'chrY' & !a$PAR
  a$region[iy] <- paste0(a$arm[iy], '_nonPAR')
  # DU145 libraries are passage-matched (P12/P13/P14) across arms; reuse the HT1080 pair helper.
  meta$passage_day_block <- meta$passage_block; meta$passage_number <- meta$passage
  core <- meta[meta$treatment %in% c('CTRL', 'LOY'), ]
  assert(nrow(core) == 6 && all(table(core$treatment) == 3), 'Expected 3 LOY and 3 CTRL libraries')
  pairs <- pair_indices(core)
  assert(nrow(pairs) == 3, 'Expected three passage-matched pairs')
  tab(meta, 'sample_manifest'); tab(pairs, 'passage_pairs')
  tab(a, 'gene_annotation_gencode49'); tab(do.call(rbind, audit), 'counting_audit')
  sums <- do.call(rbind, summaries); tab(sums, 'featurecounts_assignment_status')
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(a, row.names = a$gene_id),
    colData = S4Vectors::DataFrame(meta, row.names = meta$sample_key))
  saveRDS(se, 'outputs/models/raw_counts_SE.rds')
  list(counts = counts, annotation = a, meta = meta, core = core, pairs = pairs, assignment = sums)
}
