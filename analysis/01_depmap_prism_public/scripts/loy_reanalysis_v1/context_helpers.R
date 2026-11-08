# Helpers for the context-interaction stages. Source after scan_helpers.R.
#
# For a context C and Y status Y the model is
#   response ~ lineage + technical + Y + C + Y:C
# and the Y:C coefficient is the test of context dependence. For a binary
# context it is (Y effect inside C) minus (Y effect outside C). A continuous
# context is centred and scaled, so Y:C is the change in the Y effect per SD.

power_multiplier <- qnorm(0.975) + qnorm(0.8)

# Lines of each Y-by-context cell among the given rows; for a continuous
# context, the number of lines per Y group with the context measured.
context_cells <- function(exposed, value, type) {
  if (type == "continuous") {
    return(c(neg_in = sum(exposed == 1 & !is.na(value)), pos_in = sum(exposed == 0 & !is.na(value)), neg_out = NA, pos_out = NA))
  }
  c(
    neg_in = sum(exposed == 1 & value %in% 1), pos_in = sum(exposed == 0 & value %in% 1),
    neg_out = sum(exposed == 1 & value %in% 0), pos_out = sum(exposed == 0 & value %in% 0)
  )
}

context_supported <- function(cells, type) {
  if (type == "continuous") {
    return(all(cells[c("neg_in", "pos_in")] >= params$context$min_group_continuous))
  }
  all(cells >= params$context$min_cell)
}

# Add the context columns to an analysis frame: `ctx` (0/1, or z-score for a
# continuous context) and `interaction` = exposed * ctx. Lines without the
# context are dropped.
with_context <- function(frame, value, type) {
  out <- copy(frame)[, ctx := as.numeric(value)]
  out <- out[!is.na(ctx)]
  if (type == "continuous") out[, ctx := (ctx - mean(ctx)) / sd(ctx)]
  out[, interaction := exposed * ctx]
  out
}

# BH within each context class and across the whole family.
add_family_fdr <- function(dt, p_col, prefix) {
  dt[, (paste0(prefix, "_fdr_class")) := bh(get(p_col)), by = class]
  dt[, (paste0(prefix, "_fdr_global")) := bh(get(p_col))]
  dt[]
}
