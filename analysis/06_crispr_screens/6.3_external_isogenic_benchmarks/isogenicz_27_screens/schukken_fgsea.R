# Preranked fgsea of the Schukken et al. (2026) aneuploidy rankings, mirroring 4.3.3 (compare_with_ht1080.R):
# same 5,880 Hallmark/Reactome/GO:BP/GO:CC sets, ranking = -normZ (positive NES = sensitizer, i.e. more
# dependent in the aneuploid line), minSize 15, maxSize 500, eps 0, nPermSimple 10000, global BH per ranking,
# HGNC symbol updates with collisions dropped. Primary = matched endpoint drugZ; IsogenicZ as sensitivity.
#
#   Rscript schukken_fgsea.R
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})
set.seed(20261021)
proj <- "A:/Dropbox/LOY_paralogs/projects"
thesis <- file.path(proj, "SEPT_11_2026_THESIS")
out_dir <- file.path(thesis, "figures/chapter4/source_data/schukken_2026_rescored")
ext_dir <- file.path(thesis, "figures/chapter4/source_data/external_isogenic_screens_20261021")

current <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
update_symbols <- function(sym) {
  up <- ifelse(sym %in% current, sym, suppressWarnings(limma::alias2SymbolTable(sym, species = "Hs")))
  dup <- unique(up[!is.na(up) & duplicated(up)])
  up[up %in% dup] <- NA
  up
}

msig <- readRDS(file.path(proj, "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/resources/msigdb_all_collections.rds"))
meta <- fread(file.path(proj, "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/resources/msigdb_all_set_metadata.csv"))
pathways <- c(msig$H, msig$C2[meta[gs_subcollection == "CP:REACTOME", gs_name]],
              msig$C5[meta[gs_subcollection %in% c("GO:BP", "GO:CC"), gs_name]])
message(length(pathways), " candidate sets; fgsea keeps those with 15-500 ranked genes (5,880 in 4.3.3)")

scores <- fread(file.path(out_dir, "schukken_scores_long.csv.gz"))
summ <- fread(file.path(out_dir, "contrast_summary.csv"))
res <- rbindlist(lapply(c("matched_endpoint_drugZ", "IsogenicZ"), function(m) {
  rbindlist(lapply(summ$id, function(i) {
    x <- scores[method == m & id == i]
    x[, symbol := update_symbols(as.character(Gene))]
    x <- x[!is.na(symbol)]
    r <- sort(setNames(-x$normZ, x$symbol), decreasing = TRUE)
    f <- as.data.table(fgsea(pathways, r, minSize = 15, maxSize = 500, gseaParam = 1,
                             eps = 0, nPermSimple = 10000, nproc = 4))
    message(m, " ", i, " done (", length(r), " genes)")
    f[, .(method = m, id = i, library = summ[id == i, library], pathway, NES, pval, FDR_global = padj, size)]
  }))
}))
fwrite(res, file.path(out_dir, "schukken_fgsea_all.csv.gz"))

theme_sets <- c("GOCC_SAGA_COMPLEX", "GOCC_SWI_SNF_COMPLEX", "REACTOME_ATP_DEPENDENT_CHROMATIN_REMODELERS",
                "REACTOME_MITOCHONDRIAL_TRANSLATION", "GOCC_MITOCHONDRIAL_LARGE_RIBOSOMAL_SUBUNIT",
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_COMPLEX_I_BIOGENESIS",
                "REACTOME_TRANSLATION", "GOCC_RIBOSOME", "HALLMARK_MYC_TARGETS_V1",
                "HALLMARK_E2F_TARGETS", "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP")
ref <- fread(file.path(ext_dir, "theme_sets_all_screens.csv"))[id %in% c("HT1080 LOY", "A549 LOY")]
themes <- rbind(res[pathway %in% theme_sets], ref[, .(method = "primary drugZ", id, library = "Brunello genome-wide",
                                                       pathway, NES, pval, FDR_global, size)])
fwrite(themes, file.path(out_dir, "schukken_theme_sets.csv"))
message("done")
