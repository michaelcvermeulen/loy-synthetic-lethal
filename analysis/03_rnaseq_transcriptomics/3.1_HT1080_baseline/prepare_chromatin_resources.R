# Obtain descriptions from the same installed/cached release as the frozen membership.
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/')))
cfg <- yaml::read_yaml('config/chromatin_followup.yml')
catalogue <- as.data.frame(data.table::fread('resources/msigdb_all_set_metadata.csv'))
selected <- catalogue[grepl(cfg$screen_pattern, catalogue$gs_name, ignore.case = TRUE), ]
destination <- 'resources/msigdb_chromatin_descriptions.csv'
if (!file.exists(destination)) {
  columns <- c('gs_name', 'gs_id', 'gs_collection', 'gs_subcollection', 'gs_description',
    'gs_source_species', 'gs_pmid', 'gs_geoid', 'gs_exact_source', 'gs_url', 'db_version')
  out <- list()
  for (collection in unique(selected$gs_collection)) {
    message('Chromatin metadata: ', collection)
    z <- msigdbr::msigdbr(species = 'Homo sapiens', db_species = 'HS', collection = collection)
    out[[collection]] <- unique(as.data.frame(z[z$gs_name %in% selected$gs_name, columns]))
  }
  descriptions <- do.call(rbind, out)
  stopifnot(setequal(descriptions$gs_name, selected$gs_name),
    !anyDuplicated(descriptions$gs_name), identical(unique(descriptions$db_version), unique(catalogue$db_version)))
  data.table::fwrite(descriptions, destination)
}
# Verify existing frozen resources before extending their manifest.
manifest <- as.data.frame(data.table::fread('resources/manifest.csv'))
stopifnot(all(vapply(manifest$path, function(f) digest::digest(file = f, algo = 'sha256'), character(1)) == manifest$sha256))
if (!destination %in% manifest$path) {
  manifest <- rbind(manifest, data.frame(path = destination,
    sha256 = digest::digest(file = destination, algo = 'sha256'),
    frozen_at = format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z')))
  data.table::fwrite(manifest, 'resources/manifest.csv')
}
message('Chromatin descriptions frozen')
