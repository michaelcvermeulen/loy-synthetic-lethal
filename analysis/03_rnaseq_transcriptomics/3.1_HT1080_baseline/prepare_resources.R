# Run explicitly when first preparing or intentionally updating frozen resources.
options(timeout = 300)
root <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/'))
setwd(root)
repo <- normalizePath('../../..', winslash = '/')
dir.create('resources', showWarnings = FALSE)
cfg <- yaml::read_yaml('config/analysis.yml')
if (!file.exists('resources/msigdb_all_collections.rds')) {
  message('Freezing all available human MSigDB collections and subcollections')
  catalogue <- msigdbr::msigdbr_collections(db_species = 'HS')
  data.table::fwrite(catalogue, 'resources/msigdb_all_collection_catalogue.csv')
  all_sets <- list(); metadata <- list(); membership <- list()
  cols <- c('gs_name', 'gs_collection', 'gs_subcollection', 'gs_id', 'gs_url', 'db_version')
  for (collection in unique(catalogue$gs_collection)) {
    message('MSigDB collection ', collection)
    z <- msigdbr::msigdbr(species = 'Homo sapiens', db_species = 'HS', collection = collection)
    all_sets[[collection]] <- lapply(split(z$gene_symbol, z$gs_name), unique)
    metadata[[collection]] <- unique(as.data.frame(z[, cols]))
    membership[[collection]] <- unique(data.frame(collection = collection,
      pathway = z$gs_name, symbol = z$gene_symbol))
  }
  metadata <- do.call(rbind, metadata)
  stopifnot(!anyDuplicated(metadata$gs_name), !anyDuplicated(metadata$gs_id),
    length(unique(metadata$db_version)) == 1L,
    nrow(metadata) == sum(catalogue$num_genesets))
  saveRDS(all_sets, 'resources/msigdb_all_collections.rds')
  data.table::fwrite(metadata, 'resources/msigdb_all_set_metadata.csv')
  data.table::fwrite(data.table::rbindlist(membership), 'resources/msigdb_all_membership.csv.gz',
    compress = 'gzip')
  rm(z, all_sets, metadata, membership); gc()
}
if (!file.exists('resources/msigdb_hallmark_reactome.rds')) {
  message('Freezing MSigDB human Hallmark and Reactome with database-version metadata')
  h <- msigdbr::msigdbr(species = 'Homo sapiens', collection = 'H')
  r <- msigdbr::msigdbr(species = 'Homo sapiens', collection = 'C2', subcollection = 'CP:REACTOME')
  sets <- list(hallmark = lapply(split(h$gene_symbol, h$gs_name), unique),
               reactome = lapply(split(r$gene_symbol, r$gs_name), unique))
  saveRDS(sets, 'resources/msigdb_hallmark_reactome.rds')
  cols <- c('gs_name', 'gs_collection', 'gs_subcollection', 'gs_id', 'gs_url', 'db_version')
  data.table::fwrite(unique(rbind(as.data.frame(h[, cols]), as.data.frame(r[, cols]))),
                      'resources/msigdb_set_metadata.csv')
  data.table::fwrite(rbind(data.frame(collection = 'hallmark', pathway = h$gs_name, symbol = h$gene_symbol),
                           data.frame(collection = 'reactome', pathway = r$gs_name, symbol = r$gene_symbol)),
                     'resources/msigdb_membership.csv')
}
if (!file.exists('resources/collectri_network.csv'))
  stopifnot(file.copy(file.path(repo, cfg$collectri_source), 'resources/collectri_network.csv'))
if (!file.exists('resources/progeny_model_top100.rds'))
  saveRDS(progeny::getModel(organism = 'Human', top = 100), 'resources/progeny_model_top100.rds')
if (!file.exists('resources/hg38_cytoBand.txt')) {
  dest <- 'resources/hg38_cytoBand.txt.gz'
  if (!file.exists(dest) || file.info(dest)$size == 0)
    download.file(cfg$cytoband_url, dest, mode = 'wb', quiet = TRUE)
  z <- read.table(gzfile(dest), sep = '\t', header = FALSE)
  names(z) <- c('chrom', 'start', 'end', 'band', 'stain')
  data.table::fwrite(z, 'resources/hg38_cytoBand.txt', sep = '\t')
}
files <- setdiff(list.files('resources', full.names = TRUE, recursive = TRUE), 'resources/manifest.csv')
data.table::fwrite(data.frame(path = files,
  sha256 = vapply(files, function(f) digest::digest(file = f, algo = 'sha256'), character(1)),
  frozen_at = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z')), 'resources/manifest.csv')
message('Frozen resources ready')
