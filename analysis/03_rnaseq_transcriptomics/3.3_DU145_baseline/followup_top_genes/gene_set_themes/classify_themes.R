# Group the significant DU145 gene sets into themes by leading-edge composition.
#
# A set's theme comes from what its leading-edge genes are, not from its name. Reference gene
# groups are fixed in advance: Hallmark IFN-gamma/IFN-alpha/inflammatory, Hallmark G2M/mitotic
# spindle/E2F, the DepMap GRHL2 top-200 co-expression neighbourhood, cytosolic and mitochondrial
# ribosomal proteins, and replication-dependent histones. Reads saved outputs only.
options(stringsAsFactors = FALSE)
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', args, value = TRUE))))
root <- normalizePath(file.path(here, '..', '..'))
fg <- data.table::fread(file.path(root, 'outputs/tables/fgsea_all_collections.csv'))
H <- readRDS(file.path(root, 'resources/msigdb_all_collections.rds'))$H
nbhd <- c(data.table::fread(file.path(here, '..', 'cross_dataset/task3b_GRHL2_top200_neighbourhood.csv'))$gene, 'GRHL2')

groups <- list(
  interferon_inflammatory = unique(unlist(H[c('HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                                              'HALLMARK_INFLAMMATORY_RESPONSE')])),
  mitotic = unique(unlist(H[c('HALLMARK_G2M_CHECKPOINT', 'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_E2F_TARGETS')])),
  GRHL2_epithelial = nbhd)
is_ribo <- function(g) grepl('^(RPL[0-9]|RPS[0-9]|RPLP[0-9]|RPSA|MRPL[0-9]|MRPS[0-9]|FAU|UBA52)', g)
is_hist <- function(g) grepl('^H(1|2A|2B|3|4)C[0-9]|^H1-[0-9]|^H2A[CZ]|^H3C|^H4C', g)

sig <- fg[BH_across_collections < 0.05]
le <- strsplit(sig$leadingEdge, ';', fixed = TRUE)
frac <- function(test) vapply(le, function(g) mean(test(g)), numeric(1))
sig$frac_ribosomal <- frac(is_ribo)
sig$frac_histone <- frac(is_hist)
for (nm in names(groups)) sig[[paste0('frac_', nm)]] <- frac(function(g) g %in% groups[[nm]])

# Composition artefacts first, then biological themes by the largest reference overlap (>= 20% of the leading edge).
bio <- as.matrix(sig[, paste0('frac_', names(groups)), with = FALSE])
best <- names(groups)[max.col(bio, ties.method = 'first')]
best_frac <- apply(bio, 1, max)
sig$theme <- ifelse(sig$frac_ribosomal >= 0.4, 'ribosomal/translation',
             ifelse(sig$frac_histone >= 0.3, 'histone cluster (6p22)',
             ifelse(best_frac >= 0.2, best,
             ifelse(sig$collection == 'C1', 'positional', 'other'))))
sig$status <- ifelse(sig$interpretable, 'survives both screens',
              ifelse(sig$survives_CN_screen, 'positionally concentrated', 'copy-number sensitive'))
keep <- c('pathway', 'collection', 'subcollection', 'NES', 'BH_across_collections', 'NES_CN_stable_only', 'status', 'theme',
          'frac_ribosomal', 'frac_histone', paste0('frac_', names(groups)), 'leading_edge_top_arm', 'leading_edge_top_arm_fraction', 'leadingEdge')
out <- sig[order(theme, BH_across_collections), ..keep]
data.table::fwrite(out, file.path(here, 'significant_sets_by_theme.csv'))
summ <- out[, .(sets = .N, down = sum(NES < 0), up = sum(NES > 0),
                survive_both = sum(status == 'survives both screens')), by = .(theme)][order(-survive_both, -sets)]
data.table::fwrite(summ, file.path(here, 'theme_summary.csv'))
print(summ)
cat('\nOther (survives both screens):\n')
print(out[theme == 'other' & status == 'survives both screens', .(pathway, collection, NES = round(NES, 2), BH = signif(BH_across_collections, 2))])
