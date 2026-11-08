# Scientific checks: ambiguous membership, fixed ranks, contribution denominators and temporal weights.
testthat::test_that('HDAC partitions are mutually exclusive and conflicts take precedence', {
  src <- list(HDAC1 = 'h1', HDAC2 = 'h2', HDAC1_2 = 'both', HDAC3 = 'h3')
  m <- data.frame(pathway = c('h1_UP', 'h2_UP', 'h2_DN', 'both_UP', 'h1_UP', 'both_DN', 'h3_UP'),
    symbol = c('shared', 'shared', 'two', 'double', 'conflict', 'conflict', 'three'))
  z <- hdac_partition(m, src); rownames(z) <- z$gene_name
  testthat::expect_equal(z[c('shared', 'two', 'double', 'conflict', 'three'), 'partition'],
    c('concordant_shared', 'HDAC2_only_vs_HDAC1', 'double_knockdown_only', 'conflicting_direction', 'HDAC3_comparison_only'))
  testthat::expect_equal(nrow(z), 5L)
  testthat::expect_equal(z['conflict', 'partition_arm'], 'MIXED')
})
testthat::test_that('HDAC percentiles handle ties and never depend on other samples', {
  mat <- matrix(c(0, 0, 10, 20, 20, 10, 0, 0), 4,
    dimnames = list(c('a', 'b', 'c', 'd'), c('s1', 's2')))
  r <- hdac_percentiles(mat)
  testthat::expect_equal(unname(r[, 1]), c(100/6, 100/6, 200/3, 100))
  testthat::expect_equal(r[, 1, drop = FALSE], hdac_percentiles(mat[, 1, drop = FALSE]))
  score <- hdac_signed_score(r, c('c', 'd'), c('a', 'c'))
  testthat::expect_equal(unname(score['signed', 1]), 100 - 100/6)
  testthat::expect_equal(unname(score['signed', 2]), -100 + 100/6)
})
testthat::test_that('real HDAC summaries preserve source sets and the previous full-arm effects', {
  z <- readRDS('outputs/models/04_HDAC12_gene_drivers.rds')
  old <- readRDS('outputs/models/04_chromatin_followup.rds')
  e <- readRDS('outputs/models/02_expression.rds')
  testthat::expect_equal(z$background$gene_id, e$mapping$gene_id)
  testthat::expect_equal(sum(z$frozen$within_source_ambiguous), 14L)
  testthat::expect_true(all(z$frozen$source[z$frozen$within_source_ambiguous] == 'HDAC3'))
  testthat::expect_equal(length(unique(z$samples$sample_key)), 12L)
  for (source in names(z$cfg$sources)) {
    a <- z$influence[z$influence$source == source & z$influence$summary_arm == 'signed' & z$influence$scenario == 'original', ]
    b <- old$pairs[old$pairs$source_pair == z$cfg$sources[[source]], ]
    testthat::expect_equal(a$mean_difference, b$source_response_contrast[match(a$contrast, b$contrast)])
  }
  testthat::expect_true(all(is.na(z$evidence$leading_edge_combined[z$evidence$pathway == 'SENESE_HDAC3_TARGETS_DN'])))
  testthat::expect_equal(sum(z$coverage$genes), nrow(z$partitions))
})
testthat::test_that('HDAC contribution fractions and omissions use fixed genes without net-effect denominators', {
  z <- readRDS('outputs/models/04_HDAC12_gene_drivers.rds')
  c <- z$contributions[z$contributions$source == 'HDAC2' & z$contributions$summary_arm == 'signed', ]
  testthat::expect_equal(sum(c$absolute_effect_fraction), 1)
  testthat::expect_equal(c$contribution, c$weight * c$combined_difference)
  a <- z$influence[z$influence$source == 'HDAC2' & z$influence$summary_arm == 'signed', ]
  original <- a[a$scenario == 'original' & a$contrast == 'combined', ]
  testthat::expect_equal(sum(c$contribution), original$mean_difference)
  omitted <- a[a$scenario == 'omit_top_5', ]
  testthat::expect_equal(length(unique(omitted$removed_genes)), 1L)
  testthat::expect_equal(unique(omitted$removed_absolute_effect_fraction), sum(head(c$absolute_effect_fraction, 5)))
  genes <- z$evidence[z$evidence$source == 'HDAC2' & !z$evidence$gene_name %in% head(c$gene_name, 5), ]
  expected <- mean(genes$mean_difference[genes$source_arm == 'UP']) - mean(genes$mean_difference[genes$source_arm == 'DN'])
  testthat::expect_equal(omitted$mean_difference[omitted$contrast == 'combined'], expected)
  testthat::expect_equal(length(unique(z$shortlist$gene_id)), nrow(z$shortlist))
  testthat::expect_setequal(unique(z$shortlist$source_arm), c('UP', 'DN'))
  testthat::expect_false(any(c('pvalue', 'PValue', 'FDR', 'padj') %in% names(z$evidence)))
})
testthat::test_that('HDAC leave-one-match-out scores preserve equal period weighting', {
  z <- readRDS('outputs/models/04_HDAC12_gene_drivers.rds')
  a <- z$score_effects[z$score_effects$source == 'HDAC2' & z$score_effects$score_arm == 'signed', ]
  x <- setNames(a$percentile_difference, a$contrast)
  expected <- (mean(x[z$pairs$block[2:3]]) + mean(x[z$pairs$block[4:6]])) / 2
  testthat::expect_equal(unname(x[paste0('omit_', z$pairs$block[1])]), unname(expected))
  sample <- z$samples[z$samples$source == 'HDAC2' & z$samples$score_arm == 'signed', ]
  v <- setNames(sample$percentile_score, sample$sample_key)
  testthat::expect_equal(unname(x[z$pairs$block]), unname(v[z$pairs$LOY] - v[z$pairs$CTRL]))
})
