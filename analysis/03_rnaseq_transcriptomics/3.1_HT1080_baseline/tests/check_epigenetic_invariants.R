# Validate the extension against frozen inputs and independent arithmetic identities.
local({
  z <- readRDS('outputs/models/04_epigenetic_deep_dive.rds')
  x <- readRDS('outputs/models/01_ingest.rds'); e <- readRDS('outputs/models/02_expression.rds')
  p <- readRDS('outputs/models/03_pathways_regulators.rds')
  testthat::test_that('post hoc selection retains original coverage and global multiplicity', {
    testthat::expect_equal(nrow(z$catalogue), 339L)
    testthat::expect_equal(sum(z$catalogue$eligible), 235L)
    testthat::expect_equal(sum(z$catalogue$focus), 94L)
    index <- match(paste(z$results$pathway,z$results$contrast), paste(p$fgsea$pathway,p$fgsea$contrast))
    testthat::expect_equal(z$results$BH_across_collections, p$fgsea$BH_across_collections[index])
    testthat::expect_equal(z$results$NES, p$fgsea$NES[index])
    index <- match(paste(z$results$pathway,z$results$contrast), paste(p$camera$pathway,p$camera$contrast))
    testthat::expect_equal(z$results$camera_global_BH, p$camera$BH_across_collections[index])
    index <- match(paste(z$results$pathway,z$results$contrast), paste(p$camera_sensitivity$pathway,p$camera_sensitivity$contrast))
    testthat::expect_equal(z$results$camera_estimated_global_BH, p$camera_sensitivity$BH_across_collections[index])
    testthat::expect_false(any(z$results$camera_estimated_global_BH[z$results$contrast=='combined'] < 0.05))
  })
  testthat::test_that('mark source semantics do not conflate histone gains with RNA gains', {
    a <- z$catalogue[grepl('^ACEVEDO_',z$catalogue$pathway), ]
    testthat::expect_true(all(a$evidence_kind=='source_differential_histone_mark'))
    testthat::expect_true(all(is.na(a$response_arm)))
    a <- z$atlas[!duplicated(z$atlas$pathway), ]
    testthat::expect_false(anyDuplicated(paste(a$row_label,a$mark)) > 0)
    testthat::expect_true(all(is.na(a$NES[!a$eligible])))
    testthat::expect_equal(a$mark[a$pathway=='MEISSNER_NPC_HCP_WITH_H3K4ME2'],'K4me2')
    testthat::expect_equal(a$mark[a$pathway=='MIKKELSEN_IPS_ICP_WITH_H3K4ME3_AND_H327ME3'],'K4me3 + K27me3')
    testthat::expect_equal(a$study[a$pathway=='MIKKELSEN_MEF_HCP_WITH_H3K27ME3'],'Mikkelsen 2007')
    testthat::expect_equal(a$study[a$pathway=='MIKKELSEN_MEF_ICP_WITH_H3K27ME3'],'Mikkelsen 2008')
  })
  testthat::test_that('partitions are exhaustive and means recover complete program effects', {
    a <- z$partitions; mark <- unique(unlist(z$sets[z$cfg$mark_comparison], use.names=FALSE))
    testthat::expect_setequal(c(a$EED_only,a$shared_EED_mark),z$sets[[z$cfg$primary_program]])
    testthat::expect_setequal(c(a$mark_only,a$shared_EED_mark),mark)
    testthat::expect_equal(anyDuplicated(unlist(a,use.names=FALSE)),0L)
    testthat::expect_equal(vapply(a,length,integer(1)),c(EED_only=148L,shared_EED_mark=6L,mark_only=497L))
    d <- z$part_effects[z$part_effects$contrast=='combined' & z$part_effects$partition!='mark_only', ]
    testthat::expect_equal(weighted.mean(d$mean_difference,d$genes),
      mean(z$values[z$sets[[z$cfg$primary_program]],'combined']))
    ev <- z$evidence
    testthat::expect_true(all(ev$gene_id==e$mapping$gene_id[match(ev$gene_name,e$mapping$gene_name)]))
    testthat::expect_true(all(vapply(names(z$sets),function(name)
      setequal(z$sets[[name]],ev$gene_name[ev$pathway==name]),logical(1))))
  })
  testthat::test_that('scores use full arms and omissions keep equal period weighting', {
    raw <- x$counts[e$mapping$gene_id,x$core$sample_key,drop=FALSE]
    ranked <- apply(raw,2,function(v)100*(rank(v,ties.method='average')-1)/(length(v)-1))
    rownames(ranked) <- e$mapping$gene_name
    up <- z$sets[['PRC2_EED_DN.V1_UP']]; dn <- z$sets[['PRC2_EED_DN.V1_DN']]
    score <- colMeans(ranked[up,,drop=FALSE])-colMeans(ranked[dn,,drop=FALSE])
    actual <- z$sample_scores[z$sample_scores$program=='PRC2_EED_DN.V1__signed', ]
    testthat::expect_equal(actual$percentile_score,unname(score[actual$sample_key]))
    delta <- score[x$pairs$LOY]-score[x$pairs$CTRL]
    expected <- vapply(seq_len(6),function(j) mean(tapply(delta[-j],x$pairs$batch[-j],mean)),numeric(1))
    testthat::expect_equal(unname(z$score_values['PRC2_EED_DN.V1__signed',paste0('omit_',x$pairs$block)]),expected)
    testthat::expect_equal(nrow(z$sample_scores),12L*12L*3L)
    testthat::expect_equal(length(unique(z$sample_scores$sample_key)),12L)
    testthat::expect_equal(ncol(z$values),15L)
  })
  testthat::test_that('influence exclusions retain fixed membership and correct denominators', {
    name <- z$cfg$primary_program; genes <- z$sets[[name]]
    influence <- z$influence[z$influence$pathway==name, ]
    testthat::expect_equal(sum(influence$contribution_to_mean),mean(z$values[genes,'combined']))
    testthat::expect_equal(sum(influence$absolute_effect_fraction),1)
    keep <- setdiff(genes,head(influence$gene_name,10))
    a <- z$sensitivity[z$sensitivity$pathway==name & z$sensitivity$scenario=='omit_top_10', ]
    testthat::expect_equal(a$mean_difference,unname(colMeans(z$values[keep,a$contrast,drop=FALSE])))
    testthat::expect_true(all(a$genes==length(keep)))
    testthat::expect_true(all(a$removed_absolute_effect_fraction==sum(head(influence$absolute_effect_fraction,10))))
    empty <- hdac_summary(z$values,character(),x$pairs)
    testthat::expect_true(all(is.na(empty$mean_difference)) && all(empty$genes==0))
  })
  testthat::test_that('post hoc readouts include opposing responses without entering score definitions', {
    a <- z$readouts
    testthat::expect_true(all(a$readout_eligible & a$selected_post_hoc))
    testthat::expect_true(all(a$mean_CPM>=2 & abs(a$combined)>=0.5))
    testthat::expect_true(all(a$pass_period & a$pass_loo & a$pass_normalization))
    testthat::expect_true(all(table(a$partition,a$direction)<=3))
    testthat::expect_equal(sum(a$partition=='EED_only' & a$direction=='higher'),1L)
    testthat::expect_setequal(unique(a$direction),c('higher','lower'))
    testthat::expect_equal(z$regulators$status[z$regulators$gene_name=='PALI1'],'not_mapped_as_symbol')
    testthat::expect_equal(sum(z$regulators$gene_name=='LCOR'),1L)
    testthat::expect_equal(nrow(z$regulators),length(z$cfg$regulators))
  })
})
