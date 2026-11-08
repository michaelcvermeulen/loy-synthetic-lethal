enrichment_adjust <- function(z, p_column, family_columns = 'contrast') {
  z$BH_across_collections <- z$BH_within_collection <- z$BH_within_subcollection <- NA_real_
  adjust <- function(columns) {
    groups <- do.call(paste, c(z[, columns, drop = FALSE], sep = '\r'))
    ans <- rep(NA_real_, nrow(z))
    for (j in split(seq_len(nrow(z)), groups)) ans[j] <- p.adjust(z[[p_column]][j], 'BH', n = length(j))
    ans
  }
  z$BH_across_collections <- adjust(family_columns)
  z$BH_within_collection <- adjust(c(family_columns, 'collection'))
  z$BH_within_subcollection <- adjust(c(family_columns, 'collection', 'subcollection'))
  z
}

