# Extract same-release metadata from the installed MSigDB cache; original resource manifest is unchanged.
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/')))
source('R/00_common.R')
cfg <- yaml::read_yaml('config/epigenetic_deep_dive.yml')
catalogue <- as.data.frame(data.table::fread('resources/msigdb_all_set_metadata.csv'))
selected <- catalogue[grepl(cfg$screen_pattern, catalogue$gs_name, ignore.case = TRUE), ]
dir.create('resources/epigenetic_deep_dive', recursive = TRUE, showWarnings = FALSE)
destination <- 'resources/epigenetic_deep_dive/descriptions.csv'
if (!file.exists(destination)) {
  cols <- c('gs_name','gs_id','gs_collection','gs_subcollection','gs_description','gs_source_species',
    'gs_pmid','gs_geoid','gs_exact_source','gs_url','db_version')
  out <- lapply(unique(selected$gs_collection), function(collection) {
    message('Epigenetic source metadata: ', collection)
    # Require the cached release rather than silently obtaining different current membership.
    cached <- file.path(tools::R_user_dir('msigdbr', 'cache'), paste0('msigdb.2026.1.Hs.', collection, '.rds'))
    assert(file.exists(cached), paste('Missing frozen-release cache:', cached))
    z <- msigdbr::msigdbr(species = 'Homo sapiens', db_species = 'HS', collection = collection)
    unique(as.data.frame(z[z$gs_name %in% selected$gs_name, cols]))
  })
  descriptions <- do.call(rbind, out)
  assert(setequal(descriptions$gs_name, selected$gs_name) && !anyDuplicated(descriptions$gs_name), 'Metadata coverage mismatch')
  assert(identical(unique(descriptions$db_version), unique(catalogue$db_version)), 'MSigDB release mismatch')
  data.table::fwrite(descriptions, destination)
}
descriptions <- as.data.frame(data.table::fread(destination))
assert(setequal(descriptions$gs_name, selected$gs_name), 'Cached catalogue differs from frozen selection')
paths <- setdiff(list.files('resources/epigenetic_deep_dive', full.names = TRUE),
  'resources/epigenetic_deep_dive/manifest.csv')
data.table::fwrite(data.frame(path = paths, sha256 = unname(vapply(paths, hash_file, character(1)))),
  'resources/epigenetic_deep_dive/manifest.csv')
message('Epigenetic descriptions frozen: ', nrow(descriptions), ' sets')
