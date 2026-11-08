# Independent checks on universes, carried-over statistics, the shared rule and external polarity.
assert(length(x$genes) == 16165, paste('Common gene universe changed:', length(x$genes)))
assert(length(x$sets) == 22961, paste('Common gene-set universe changed:', length(x$sets)))
assert(!anyDuplicated(x$fg$HT1080$pathway) && all(x$fg$HT1080$contrast == 'combined'), 'HT1080 gene sets are not restricted to the combined contrast')
# Carried-over statistics equal the per-line source tables exactly.
for (L in LINES) {
  src <- x$fg[[L]]; d <- gs$sets
  i <- match(d$pathway, src$pathway)
  assert(isTRUE(all.equal(d[[paste0('NES_', L)]], src$NES[i])) && isTRUE(all.equal(d[[paste0('BH_', L)]], src$BH_across_collections[i])),
         paste('NES/BH not carried over unchanged for', L))
}
# The shared rule reproduces from the carried-over values.
nes <- as.matrix(gs$sets[, paste0('NES_', LINES)]); bh <- as.matrix(gs$sets[, paste0('BH_', LINES)])
rule <- rowSums(bh < cfg$alpha, na.rm = TRUE) >= 2 & apply(sign(nes), 1, function(v) length(unique(v)) == 1)
assert(identical(unname(rule), unname(gs$sets$class == 'shared')), 'Shared classification does not reproduce')
n_shared <- sum(gs$sets$class == 'shared'); n_div <- sum(gs$sets$class == 'divergent')
assert(n_shared == 46 && n_div == 23, sprintf('Shared/divergent counts differ from the survey: %d / %d', n_shared, n_div))
# Gene-level effects equal the per-line tables.
for (L in LINES) {
  e <- x$eff[[L]]
  assert(isTRUE(all.equal(gl$genes[[paste0('md_', L)]], e$mean_difference[match(gl$genes$gene_id, e$gene_id)])), paste('Gene effect mismatch for', L))
}
# External polarity: a Y gene must be lower with Y absent wherever it is tested.
pc <- ex$polarity[!is.na(ex$polarity$polarity_ok), ]
assert(all(pc$polarity_ok), paste('Polarity check failed for:', paste(pc$dataset[!pc$polarity_ok], collapse = ', ')))
assert(all(LINES %in% pc$dataset), 'In-house polarity not verified')
assert(all(c('PageLab_Fib_X0vXY', 'PageLab_LCL_X0vXY', 'XYY_integrated_flipped', 'GSE313572_HT1080_ctrl_YnegvYpos',
                  'Tahoe_A549_DMSO_YnegvYpos', 'Tahoe_J82_DMSO_YnegvYpos', 'Tahoe_SKMEL2_DMSO_YnegvYpos', 'Tahoe_HCT15_DMSO_YnegvYpos') %in% pc$dataset), 'Expected Y-gene polarity checks missing')
reg <- data.table::fread(cfg$external_registry)
assert(reg$polarity[reg$key == 'XYY_integrated_flipped'] == -1, 'XYY dataset is not sign-flipped')
hashes <- data.table::fread('outputs/tables/input_hashes.csv')
assert(all(vapply(hashes$path, hash_file, character(1)) == hashes$sha256), 'Inputs changed during the run')
jsonlite::write_json(list(status = 'passed', checked_at = format(Sys.time(), '%Y-%m-%d %H:%M:%S %z'),
  common_genes = length(x$genes), common_gene_sets = length(x$sets), shared_sets = n_shared, divergent_sets = n_div,
  polarity_checked_datasets = nrow(pc), shared_rule = cfg$shared_rule), 'provenance/validation.json', auto_unbox = TRUE, pretty = TRUE)
message('Scientific validation passed')
