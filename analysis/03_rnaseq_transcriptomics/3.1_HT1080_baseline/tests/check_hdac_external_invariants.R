# These checks validate scientific contracts against the completed external and original models.
local({
  source <- readRDS('outputs/models/04_HDACEXT_source_models.rds')
  z <- readRDS('outputs/models/04_HDACEXT_comparison.rds')
  original <- readRDS('outputs/models/02_expression.rds')
  x <- readRDS('outputs/models/01_ingest.rds')
  cfg <- source$cfg
  testthat::test_that('external source counts, exclusions and within-stratum designs retain provenance', {
    raw <- as.data.frame(data.table::fread(cfg$counts))
    testthat::expect_equal(unname(source$counts), unname(as.matrix(raw[, source$metadata$sample_key])))
    testthat::expect_equal(dim(source$counts), c(60683L, 69L))
    testthat::expect_equal(sum(source$metadata$used_in_model), 63L)
    testthat::expect_setequal(unique(source$metadata$source_group[!source$metadata$used_in_model]), c('HAP_MS1', 'HAP_MS3'))
    testthat::expect_equal(length(unique(source$metadata$GSM)), 69L)
    testthat::expect_equal(source$model_audit$residual_df, c(22L, 20L))
    cc <- source$contrasts
    valid <- vapply(seq_len(nrow(cc)), function(i) {
      fit <- source$fits[[cc$sample_ID_stratum[i]]]
      d <- fit$design
      # Each contrast is a difference of actual source groups, without replicate-number pairing.
      con <- setNames(numeric(ncol(d)), colnames(d))
      con[cc$numerator[i]] <- 1; con[cc$denominator[i]] <- -1
      n <- fit$metadata$model_group == cc$numerator[i]
      k <- fit$metadata$model_group == cc$denominator[i]
      all(d %in% c(0, 1)) && all(rowSums(d) == 1) && qr(d)$rank == ncol(d) &&
        sum(n) == 3 && sum(k) == 3 && all((d %*% con)[n] == 1) && all((d %*% con)[k] == -1)
    }, logical(1))
    testthat::expect_true(all(valid))
    testthat::expect_equal(nrow(cc), 27L)
  })
  testthat::test_that('external gene sets use only frozen source statistics and keep both arms', {
    valid <- vapply(source$coverage$pathway, function(path) {
      a <- source$coverage[source$coverage$pathway == path, ]
      r <- source$results[source$results$contrast == a$external_contrast, ]
      expected <- r[is.finite(r$external_FDR_nonzero) & r$external_FDR_nonzero < cfg$signature_FDR_nonzero &
        abs(r$logFC) >= cfg$signature_min_abs_logFC & (r$logFC > 0) == (a$source_arm == 'UP'), ]
      expected <- expected[order(-expected$F, expected$gene_id), ]
      actual <- source$signatures[source$signatures$pathway == path, ]
      identical(actual$gene_id, expected$gene_id) &&
        identical(actual$selected_top250, seq_len(nrow(expected)) <= cfg$signature_cap_per_arm)
    }, logical(1))
    testthat::expect_true(all(valid))
    testthat::expect_true(all(vapply(split(source$results, source$results$contrast), function(r)
      isTRUE(all.equal(r$external_FDR_nonzero, p.adjust(r$PValue, 'BH'))), logical(1))))
    testthat::expect_equal(length(z$sets), 54L)
    testthat::expect_true(all(source$coverage$eligible))
    testthat::expect_true(all(vapply(z$contrasts$contrast, function(n)
      length(intersect(z$sets[[paste0(n, '__UP')]], z$sets[[paste0(n, '__DN')]])) == 0L, logical(1))))
  })
  testthat::test_that('cross-study mapping uses a unique shared locus background and unchanged HT1080 effects', {
    m <- z$common
    testthat::expect_equal(nrow(m), 14881L)
    testthat::expect_false(anyDuplicated(m$gene_id) > 0L || anyDuplicated(m$HT1080_symbol) > 0L)
    testthat::expect_true(all(m$expressed_S47 & m$expressed_S52 & m$chromosome_agrees))
    testthat::expect_equal(m$gene_id, original$mapping$stable_gene_id[match(m$HT1080_gene_id, original$mapping$gene_id)])
    testthat::expect_equal(unname(z$HT1080_effects[, x$pairs$block]), unname(original$differences[m$HT1080_gene_id, ]))
    combined <- rowMeans(cbind(rowMeans(z$HT1080_effects[, x$pairs$block[1:3]]),
      rowMeans(z$HT1080_effects[, x$pairs$block[4:6]])))
    loo1 <- (rowMeans(z$HT1080_effects[, x$pairs$block[2:3]]) + rowMeans(z$HT1080_effects[, x$pairs$block[4:6]])) / 2
    testthat::expect_equal(unname(z$HT1080_effects[, 'combined']), unname(combined))
    testthat::expect_equal(unname(z$HT1080_effects[, paste0('omit_', x$pairs$block[1])]), unname(loo1))
    hashes <- as.data.frame(data.table::fread('provenance/HDACEXT_frozen_input_hashes.csv'))
    testthat::expect_equal(unname(vapply(hashes$path, hash_file, character(1))), hashes$sha256)
  })
  testthat::test_that('both response arms and all masks use the declared genes and exposures', {
    valid <- vapply(names(z$masks), function(scenario) {
      keep <- setdiff(z$common$gene_id, z$masks[[scenario]])
      a <- z$arm_means[z$arm_means$scenario == scenario, ]
      all(vapply(seq_len(nrow(a)), function(i) {
        n <- a$external_contrast[i]; c <- a$contrast[i]
        u <- intersect(z$sets[[paste0(n, '__UP')]], keep)
        d <- intersect(z$sets[[paste0(n, '__DN')]], keep)
        isTRUE(all.equal(a$UP_mean[i], mean(z$HT1080_effects[u, c]))) &&
          isTRUE(all.equal(a$DN_mean[i], mean(z$HT1080_effects[d, c]))) &&
          a$UP_genes[i] == length(u) && a$DN_genes[i] == length(d)
      }, logical(1)))
    }, logical(1))
    testthat::expect_true(all(valid))
    testthat::expect_equal(z$arm_means$signed_mean_difference, z$arm_means$UP_mean - z$arm_means$DN_mean)
    testthat::expect_equal(z$arm_means$source_like_two_arm_pattern, z$arm_means$UP_mean > 0 & z$arm_means$DN_mean < 0)
    testthat::expect_equal(z$enrichment$size, z$enrichment$camera_fixed_genes)
    testthat::expect_equal(z$enrichment$size, z$enrichment$camera_estimated_genes)
    testthat::expect_true(all(vapply(seq_len(nrow(z$enrichment)), function(i) {
      a <- z$enrichment[i, ]; retained <- setdiff(z$sets[[a$pathway]], z$masks[[a$scenario]])
      le <- strsplit(a$leadingEdge, ';', fixed = TRUE)[[1]]
      length(retained) == a$size && all(le %in% retained) &&
        a$universe_genes == nrow(z$common) - length(z$masks[[a$scenario]])
    }, logical(1))))
  })
  testthat::test_that('new enrichment families, numerical limitations and missing readouts are explicit', {
    en <- z$enrichment
    testthat::expect_equal(nrow(en), 648L)
    testthat::expect_false(anyNA(en$NES) || anyNA(en$pval) || anyNA(en$camera_fixed_P) || anyNA(en$camera_estimated_P))
    testthat::expect_true(all(en$family_size[en$scenario == 'original'] == 54))
    testthat::expect_true(all(is.na(en$family_size[en$scenario != 'original'])))
    valid <- vapply(c('fgsea', 'camera_fixed', 'camera_estimated'), function(method) {
      qcol <- paste0(method, '_secondary_BH'); pcol <- if (method == 'fgsea') 'pval' else paste0(method, '_P')
      all(is.na(en[[qcol]][en$scenario != 'original'])) &&
        all(vapply(split(en[en$scenario == 'original', ], en$contrast[en$scenario == 'original']), function(a)
          isTRUE(all.equal(a[[qcol]], p.adjust(a[[pcol]], 'BH', n = 54))), logical(1)))
    }, logical(1))
    testthat::expect_true(all(valid))
    testthat::expect_equal(en$fgsea_precision_floor, en$pval <= cfg$fgsea_eps)
    testthat::expect_true(nrow(z$numerical_warnings) > 0)
    testthat::expect_false(any(grepl('pvalue|p_value|confidence', names(z$correlations), ignore.case = TRUE)))
    testthat::expect_equal(nrow(z$readouts), 216L)
    testthat::expect_true(all(is.na(z$readouts$external_logFC[z$readouts$HT1080_symbol == 'BEX1'])))
    testthat::expect_true(all(z$readouts$external_coverage[z$readouts$HT1080_symbol == 'BEX1'] == 'not_in_common_expression_universe'))
  })
})
