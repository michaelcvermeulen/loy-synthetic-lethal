# Run from the analysis directory, with R helpers sourced. Tests target scientific failure modes.
testthat::test_that('source experiment definitions override misleading chromatin set names', {
  z <- readRDS('outputs/models/04_chromatin_followup.rds')
  d <- z$definitions
  testthat::expect_equal(d$source_perturbation[match(c('PRC2_EZH2_UP.V1_UP','PRC1_BMI_UP.V1_UP'),d$gs_name)],
    c('knockdown','knockdown'))
  testthat::expect_equal(d$source_perturbation[d$gs_name=='WANG_SMARCE1_TARGETS_UP'],'ectopic_expression')
  testthat::expect_equal(d$evidence_kind[d$gs_name=='MIKKELSEN_MCV6_HCP_WITH_H3K27ME3'],
    'occupancy_or_histone_mark_in_source_cells')
  testthat::expect_equal(d$gs_source_species[d$gs_name=='MIKKELSEN_MCV6_HCP_WITH_H3K27ME3'],'MM')
})

testthat::test_that('paired source-response arms exclude conflicting genes and retain direction', {
  effects <- matrix(c(3,-1,100,4,-2,200),nrow=3,
    dimnames=list(c('up','down','conflict'),c('August','May')))
  sets <- list(U=c('up','conflict'),D=c('down','conflict'))
  z <- chromatin_pair_means(effects,sets,'U','D')
  testthat::expect_equal(z$source_response_contrast,c(4,6))
  testthat::expect_equal(z$ambiguous_genes_excluded,c(1L,1L))
  testthat::expect_true(all(z$two_arm_pattern=='source_like_two_arm_pattern'))
  testthat::expect_equal(chromatin_pair_means(effects,sets,'D','U')$source_response_contrast,c(-4,-6))
})

testthat::test_that('chromatin sensitivities preserve original results and remove declared genes', {
  z <- readRDS('outputs/models/04_chromatin_followup.rds')
  original <- z$results[z$results$contrast=='combined',]
  baseline <- z$sensitivity[z$sensitivity$scenario=='original',]
  ix <- match(baseline$pathway,original$pathway)
  testthat::expect_equal(baseline$NES,original$NES[ix])
  testthat::expect_equal(baseline$pval,original$pval[ix])
  testthat::expect_equal(baseline$camera_P,original$camera_P[ix])
  testthat::expect_false(any(c('padj','BH_across_collections','FDR') %in% names(z$sensitivity)))
  universe <- data.table::fread('outputs/tables/chromatin_sensitivity_gene_universes.csv')
  removed <- unique(unlist(z$programs))
  testthat::expect_false(any(universe[scenario=='without_IFN_ECM_cycle' & retained]$gene_name %in% removed))
  testthat::expect_equal(z$sensitivity$retained_genes+z$sensitivity$removed_genes,
    baseline$retained_genes[match(z$sensitivity$pathway,baseline$pathway)])
  testthat::expect_true(all(z$sensitivity$status[z$sensitivity$retained_genes<15]=='below_minimum_represented_genes'))
  testthat::expect_equal(length(unique(z$regulator_matches$gene_id)),34L)
  testthat::expect_equal(nrow(z$regulator_matches),34L*6L)
  testthat::expect_equal(length(unique(z$pairs$contrast[z$pairs$summary_type=='individual_passage_match'])),6L)
})

testthat::test_that('CAMERA contrasts recover matched period effects without aliasing period and block', {
  x <- readRDS('outputs/models/01_ingest.rds'); meta <- x$core
  model <- camera_design(meta)
  testthat::expect_equal(model$residual_df, 4L)
  block_effect <- as.numeric(factor(meta$passage_day_block)) * 10
  y <- block_effect + (meta$treatment == 'LOY') * ifelse(meta$batch == 'AUG_2025', 2, 6)
  coef <- qr.coef(qr(model$design), y)
  testthat::expect_equal(as.numeric(crossprod(model$contrasts, coef)), c(4, 2, 6), tolerance = 1e-10)
  perm <- rev(seq_len(nrow(meta))); permuted <- camera_design(meta[perm, ])
  testthat::expect_equal(as.numeric(crossprod(permuted$contrasts, qr.coef(qr(permuted$design), y[perm]))),
    c(4, 2, 6), tolerance = 1e-10)
})

testthat::test_that('expanded tests cover the complete catalogue with a shared eligible universe', {
  p <- readRDS('outputs/models/03_pathways_regulators.rds')
  catalogue <- data.table::fread('resources/msigdb_all_collection_catalogue.csv')
  testthat::expect_setequal(unique(p$set_coverage$collection), unique(catalogue$gs_collection))
  testthat::expect_equal(nrow(p$set_coverage), sum(catalogue$num_genesets))
  testthat::expect_false(anyDuplicated(p$set_coverage$pathway) > 0)
  eligible <- p$set_coverage$pathway[p$set_coverage$eligible]
  for (contrast in c('combined', 'AUG_2025', 'MAY_2026')) {
    fg <- p$fgsea[p$fgsea$contrast == contrast, ]; cam <- p$camera[p$camera$contrast == contrast, ]
    testthat::expect_setequal(fg$pathway, eligible)
    testthat::expect_setequal(cam$pathway, eligible)
    testthat::expect_false(anyDuplicated(cam$pathway) > 0)
    testthat::expect_equal(fg$BH_across_collections, p.adjust(fg$pval, 'BH', n = length(eligible)))
    testthat::expect_equal(cam$BH_across_collections, p.adjust(cam$PValue, 'BH', n = length(eligible)))
    testthat::expect_equal(cam$NGenes[match(fg$pathway, cam$pathway)], fg$size)
  }
  model <- readRDS('outputs/models/camera_voom_model.rds'); e <- readRDS('outputs/models/02_expression.rds')
  testthat::expect_identical(colnames(model$voom$E), rownames(model$design))
  testthat::expect_equal(unname(model$effective_library_reads), unname(e$norm$effective))
  testthat::expect_identical(rownames(model$voom$E), e$mapping$gene_name)
  testthat::expect_true(all(is.finite(model$voom$weights) & model$voom$weights > 0))
})

testthat::test_that('missing enrichment P values do not silently shrink the correction family', {
  z <- data.frame(contrast = 'combined', collection = c('H', 'H', 'C2'), subcollection = '', p = c(0.01, NA, 0.5))
  ans <- enrichment_adjust(z, 'p')
  testthat::expect_equal(ans$BH_across_collections, c(0.03, NA, 0.75))
  testthat::expect_equal(ans$BH_within_collection, c(0.02, NA, 0.5))
})

testthat::test_that('ZMYND8 omission diagnostics preserve the original full-regulon score and raw counts', {
  z <- readRDS('outputs/models/04_ZMYND8_followup.rds'); x <- readRDS('outputs/models/01_ingest.rds')
  full <- z$omission[z$omission$scenario == 'all_targets', ]
  testthat::expect_equal(full$score, z$original$score[match(full$contrast, z$original$contrast)])
  testthat::expect_equal(z$omission$remaining_targets[z$omission$scenario == 'omit_MMP1_MMP3'], rep(6L, 3))
  testthat::expect_equal(z$sample$raw_reads,
    as.numeric(x$counts[cbind(match(z$sample$gene_id, rownames(x$counts)), match(z$sample$sample_key, colnames(x$counts)))]))
  testthat::expect_equal(z$targets$signed_combined_effect, z$targets$mor * z$targets$combined)
})

testthat::test_that('sample order cannot change passage contrasts', {
  meta <- data.frame(sample_key = c('c1','l1','c2','l2'), treatment = rep(c('CTRL','LOY'),2),
    batch = 'AUG_2025', passage_number = rep(1:2,each=2), passage_day_block = rep(c('b1','b2'),each=2))
  mat <- matrix(c(10,13,100,107),nrow=1,dimnames=list('g',meta$sample_key))
  pairs <- pair_indices(meta)
  testthat::expect_equal(as.numeric(matched_differences(mat,pairs)),c(3,7))
  perm <- c(4,1,3,2)
  testthat::expect_equal(matched_differences(mat[,perm,drop=FALSE],pairs),matched_differences(mat,pairs))
  testthat::expect_equal(align_samples(mat[,perm,drop=FALSE],meta),mat)
  bad <- rbind(meta,meta[1,]); testthat::expect_error(pair_indices(bad),'duplicated passage match')
})
testthat::test_that('period weights remain equal with unequal match counts', {
  z <- matrix(c(0,0,6,6,6),nrow=1,dimnames=list('g',letters[1:5]))
  pairs <- data.frame(batch=c('August','August','May','May','May'))
  testthat::expect_equal(as.numeric(rowMeans(period_average(z,pairs))),3)
  testthat::expect_false(isTRUE(all.equal(mean(z),3)))
})
testthat::test_that('symbol representatives retain locus identity and ignore differential statistics', {
  a <- data.frame(gene_id=c('ENSG2.1','ENSG1.1','ENSG3.1'),gene_name=c('DUP','DUP','ONE'))
  z <- choose_symbols(a,c(4,4,1),rep(TRUE,3))
  testthat::expect_equal(z$gene_id[z$representative],c('ENSG1.1','ENSG3.1'))
  testthat::expect_equal(nrow(z),3L)
  testthat::expect_equal(stable_id(c('ENSG1.9_PAR_Y','ENSG1.9')),c('ENSG1_PAR_Y','ENSG1'))
})
testthat::test_that('paralog totals use abundance denominators and flag absent baselines', {
  testthat::expect_equal(relative_pair_total(10,0,9,1),1)
  testthat::expect_equal(relative_pair_total(2,0,1,9),0.2)
  testthat::expect_true(is.na(relative_pair_total(1,1,0,0)))
})
testthat::test_that('real run preserves the counting annotation and confirmed replication structure', {
  x <- readRDS('outputs/models/01_ingest.rds'); e <- readRDS('outputs/models/02_expression.rds')
  testthat::expect_identical(rownames(x$counts),x$annotation$gene_id)
  testthat::expect_equal(nrow(x$core),12L)
  testthat::expect_equal(length(unique(x$core$starting_population_id)),2L)
  testthat::expect_equal(length(unique(x$core$generation_experiment_id)),1L)
  testthat::expect_equal(nrow(pair_indices(x$core)),6L)
  testthat::expect_true(all(x$annotation$arm[x$annotation$gene_name %in% c('BRCA2','RB1')] == 'chr13q'))
  testthat::expect_false(any(c('pvalue','p.value','p','fdr','FDR','padj') %in% names(e$effects)))
  testthat::expect_equal(length(unique(e$mapping$gene_name)),nrow(e$mapping))
  testthat::expect_true(all(is.finite(e$logCPM)))
})
testthat::test_that('August signature provenance excludes May selection and scores samples independently', {
  f <- readRDS('outputs/models/August_frozen_signature.rds')
  testthat::expect_true(all(grepl('AUG_2025',f$discovery_samples)))
  testthat::expect_false(any(grepl('MAY_2026',f$discovery_samples)))
  x <- readRDS('outputs/models/01_ingest.rds'); d <- readRDS('outputs/models/04_dosage_repeatability.rds')
  s <- x$core$sample_key[12]; v <- x$counts[f$background_ids,s]
  ranks <- (rank(v,ties.method='average')-1)/(length(v)-1)
  z <- mean(ranks[f$selected$gene_id[f$selected$direction>0]]) - mean(ranks[f$selected$gene_id[f$selected$direction<0]])
  testthat::expect_equal(z,d$signature_scores$signed_percentile_score[d$signature_scores$sample_key==s])
})
