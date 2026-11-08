# Runtime check only: these Monte Carlo results are not used for inference.
library(fgsea)
library(data.table)
r <- fread('outputs/tables/pathway_TF_effect_ranking.csv')
m <- fread('outputs/tables/chromatin_atlas/gene_set_membership_RNA_universe.csv')
g <- m[set_name == 'H3K4me3_PROMOTER', gene_name]
v <- setNames(r$combined, r$gene_name)
v <- v[order(-v, names(v))]
set.seed(42)
print(system.time(invisible(fgseaMultilevel(list(runtime_check = g), v,
  minSize = 15, maxSize = length(v) - 1L, nPermSimple = 1000, eps = 1e-10,
  BPPARAM = BiocParallel::SerialParam()))))
