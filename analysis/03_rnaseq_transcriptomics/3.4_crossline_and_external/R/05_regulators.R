# CollecTRI target-pattern scores classified with the same shared rule.
regulators <- function(x, cfg) {
  src <- Reduce(intersect, lapply(x$tf, `[[`, 'source'))
  d <- data.frame(source = src)
  for (L in LINES) {
    t <- x$tf[[L]]
    d[[paste0('score_', L)]] <- t$score[match(src, t$source)]
    d[[paste0('BH_', L)]] <- t$BH_target_pattern[match(src, t$source)]
  }
  sc <- as.matrix(d[, paste0('score_', LINES)]); bh <- as.matrix(d[, paste0('BH_', LINES)])
  sig <- bh < cfg$alpha; sig[is.na(sig)] <- FALSE
  same <- apply(sign(sc), 1, function(v) length(unique(v)) == 1)
  opp <- sapply(seq_len(nrow(d)), function(i) length(unique(sign(sc[i, sig[i, ]]))) > 1)
  d$n_significant <- rowSums(sig)
  d$class <- ifelse(d$n_significant >= 2 & same, 'shared', ifelse(opp, 'divergent',
             ifelse(d$n_significant == 1, 'line_specific', ifelse(d$n_significant >= 2, 'agree_where_significant_third_opposite', 'not_significant'))))
  d$caveat <- ifelse(d$source %in% c('GRHL2', 'ZEB1'), 'score partly restates the GRHL2 transcript change (self/GRHL2 edge; see DU145 follow-up)', '')
  d <- d[order(match(d$class, c('shared', 'divergent')), -d$n_significant, d$source), ]
  tab(d, 'regulators_all_common'); tab(d[d$class %in% c('shared', 'divergent'), ], 'regulators_shared_divergent')
  d
}
