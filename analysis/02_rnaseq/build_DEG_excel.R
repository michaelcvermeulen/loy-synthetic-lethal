#!/usr/bin/env Rscript
# Build a multi-sheet Excel workbook of significant DEGs across all three
# cell lines and every available contrast. Run from the repo root or from
# this directory.
#
# Significance threshold: FDR < 0.05 AND |logFC| > 0.75 (project defaults).
# Output: analysis/02_rnaseq/significant_DEGs.xlsx

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(openxlsx)
})

# FDR threshold can be overridden from the command line:
#   Rscript build_DEG_excel.R 0.05
#   Rscript build_DEG_excel.R 0.1
args  <- commandArgs(trailingOnly = TRUE)
ALPHA <- if (length(args) >= 1) as.numeric(args[[1]]) else 0.05
stopifnot(is.finite(ALPHA), ALPHA > 0, ALPHA < 1)
LFC_THR <- 0.75

BASE_DIR <- "a:/Dropbox/LOY_paralogs/projects/github/loy-synthetic-lethality/analysis/02_rnaseq"
# Encode the FDR in the filename so multiple thresholds coexist
fdr_tag  <- sub("\\.", "p", sprintf("%.2g", ALPHA))   # 0.05 -> "0p05", 0.1 -> "0p1"
OUT_XLSX <- file.path(BASE_DIR,
                      sprintf("significant_DEGs_FDR%s.xlsx", fdr_tag))

# (cell_line, contrast_label, file_path) tuples to ingest. LNCaP files were
# written with a "ht1080_" prefix due to a template-string carryover in the
# LNCaP Rmd — the content is LNCaP, only the filenames are mislabeled.
de_jobs <- bind_rows(
  data.frame(cell_line = "HT1080",
             contrast  = c("LOY_v_sgCTRL_paired", "LOY_v_sgCTRL_unpaired",
                           "LOY_v_pooled",        "LOY_v_PARENTAL",
                           "sgCTRL_v_PARENTAL"),
             file      = c("ht1080_LOY_v_sgCTRL_DE_paired.csv",
                           "ht1080_LOY_v_sgCTRL_DE_unpaired.csv",
                           "ht1080_LOY_v_pooled_DE.csv",
                           "ht1080_LOY_v_PARENTAL_DE.csv",
                           "ht1080_sgCTRL_v_PARENTAL_DE.csv"),
             stringsAsFactors = FALSE),
  data.frame(cell_line = "A549",
             contrast  = c("LOY_v_sgCTRL_paired", "LOY_v_sgCTRL_unpaired",
                           "LOY_v_pooled",        "LOY_v_PARENTAL",
                           "sgCTRL_v_PARENTAL"),
             file      = c("a549_LOY_v_sgCTRL_DE_paired.csv",
                           "a549_LOY_v_sgCTRL_DE_unpaired.csv",
                           "a549_LOY_v_pooled_DE.csv",
                           "a549_LOY_v_PARENTAL_DE.csv",
                           "a549_sgCTRL_v_PARENTAL_DE.csv"),
             stringsAsFactors = FALSE),
  data.frame(cell_line = "LNCaP",
             contrast  = c("LOY_v_sgCTRL_paired", "LOY_v_sgCTRL_unpaired"),
             file      = c("ht1080_LOY_v_sgCTRL_DE_paired.csv",
                           "ht1080_LOY_v_sgCTRL_DE_unpaired.csv"),
             stringsAsFactors = FALSE))

de_jobs$path <- file.path(BASE_DIR, de_jobs$cell_line, "DE", de_jobs$file)
stopifnot(all(file.exists(de_jobs$path)))

# Load + filter each DE table; tag direction.
load_sig <- function(path, cell_line, contrast) {
  d <- fread(path)
  # Make sure we have logFC + FDR + gene_name columns regardless of which
  # contrast wrote the file. Different write paths exist (run_de_contrast
  # writes logFC/FDR/PValue/t/logCPM; some pooled writers may differ).
  if (!"FDR" %in% names(d)) {
    if ("adj.P.Val" %in% names(d)) setnames(d, "adj.P.Val", "FDR")
  }
  if (!"PValue" %in% names(d) && "P.Value" %in% names(d)) {
    setnames(d, "P.Value", "PValue")
  }
  if (!"gene_name" %in% names(d) && "gene" %in% names(d)) {
    setnames(d, "gene", "gene_name")
  }
  needed <- c("gene_name", "logFC", "FDR")
  if (!all(needed %in% names(d))) {
    stop("Missing columns in ", path, ": ",
         paste(setdiff(needed, names(d)), collapse = ", "))
  }
  d <- d[!is.na(FDR) & !is.na(logFC) & FDR < ALPHA & abs(logFC) > LFC_THR]
  d[, direction := ifelse(logFC > 0, "UP", "DOWN")]
  d[, cell_line := cell_line]
  d[, contrast  := contrast]
  # Keep a tidy column subset, in a consistent order.
  keep_cols <- c("cell_line", "contrast", "direction", "gene_name",
                 "logFC", "FDR")
  if ("PValue" %in% names(d)) keep_cols <- c(keep_cols, "PValue")
  if ("logCPM" %in% names(d)) keep_cols <- c(keep_cols, "logCPM")
  if ("t"      %in% names(d)) keep_cols <- c(keep_cols, "t")
  d <- d[, ..keep_cols]
  # setorder() doesn't accept function-of-column expressions; use base order()
  d <- d[order(direction, -abs(logFC))]
  d
}

all_sig <- mapply(load_sig,
                  path      = de_jobs$path,
                  cell_line = de_jobs$cell_line,
                  contrast  = de_jobs$contrast,
                  SIMPLIFY  = FALSE)
names(all_sig) <- paste(de_jobs$cell_line, de_jobs$contrast, sep = "_")

# Long-format master table
master <- rbindlist(all_sig, fill = TRUE)
cat(sprintf("Master long-format table: %d rows across %d (cell x contrast) combos\n",
            nrow(master), length(all_sig)))

# Summary counts. Use a full (cell, contrast) skeleton so combos with 0 sig
# DEGs (e.g. A549 LOY-vs-sgCTRL where the sgRNA control is non-quiet) still
# show up as explicit zeros instead of being dropped silently.
skeleton <- as.data.table(de_jobs)[, .(cell_line, contrast)]
summary_tbl <- master[, .(n_genes = .N),
                      by = .(cell_line, contrast, direction)]
summary_wide <- dcast(summary_tbl, cell_line + contrast ~ direction,
                      value.var = "n_genes", fill = 0L)
summary_wide <- merge(skeleton, summary_wide,
                      by = c("cell_line", "contrast"), all.x = TRUE)
for (col in c("UP", "DOWN")) {
  if (!col %in% names(summary_wide)) summary_wide[[col]] <- 0L
  summary_wide[is.na(get(col)), (col) := 0L]
}
summary_wide[, total := UP + DOWN]
setorder(summary_wide, cell_line, contrast)

# Write workbook
wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold",
                            fgFill = "#11335d", fontColour = "#FFFFFF",
                            halign = "center", border = "TopBottom",
                            borderColour = "#11335d")
num_style <- createStyle(numFmt = "0.000")

# Sheet 1: Summary counts
addWorksheet(wb, "Summary", tabColour = "#11335d")
writeData(wb, "Summary", summary_wide, headerStyle = header_style)
freezePane(wb, "Summary", firstRow = TRUE)
setColWidths(wb, "Summary", cols = 1:ncol(summary_wide),
             widths = c(12, 28, 8, 8, 8))

# Sheet 2: Master long format
addWorksheet(wb, "All_genes_long", tabColour = "#b89d29")
writeData(wb, "All_genes_long", master, headerStyle = header_style)
freezePane(wb, "All_genes_long", firstRow = TRUE)
addFilter(wb, "All_genes_long", row = 1, cols = seq_len(ncol(master)))
setColWidths(wb, "All_genes_long", cols = seq_len(ncol(master)),
             widths = c(12, 28, 10, 14, 10, 12, 12,
                        rep(10, max(0, ncol(master) - 7))))
# Numeric formatting for logFC, FDR, PValue columns
num_cols <- which(names(master) %in% c("logFC", "FDR", "PValue",
                                       "logCPM", "t"))
if (length(num_cols) > 0) {
  addStyle(wb, "All_genes_long", num_style,
           rows = 2:(nrow(master) + 1), cols = num_cols,
           gridExpand = TRUE)
}

# Sheets 3..N: one per (cell x contrast). Tab name <= 31 chars (Excel limit).
shorten <- function(s) substr(gsub("_", "-", s), 1, 31)
for (i in seq_along(all_sig)) {
  nm   <- names(all_sig)[i]
  tab  <- shorten(nm)
  df   <- all_sig[[i]]
  cell <- de_jobs$cell_line[i]
  colour <- switch(cell,
                   "HT1080" = "#FFD580",
                   "A549"   = "#A8D5BA",
                   "LNCaP"  = "#D4B7E0",
                   "#cccccc")
  addWorksheet(wb, tab, tabColour = colour)
  writeData(wb, tab, df, headerStyle = header_style)
  freezePane(wb, tab, firstRow = TRUE)
  addFilter(wb, tab, row = 1, cols = seq_len(ncol(df)))
  setColWidths(wb, tab, cols = seq_len(ncol(df)),
               widths = c(12, 28, 10, 14, 10, 12, 12,
                          rep(10, max(0, ncol(df) - 7))))
  if (nrow(df) > 0) {
    nc_idx <- which(names(df) %in% c("logFC", "FDR", "PValue",
                                     "logCPM", "t"))
    if (length(nc_idx) > 0) {
      addStyle(wb, tab, num_style,
               rows = 2:(nrow(df) + 1), cols = nc_idx,
               gridExpand = TRUE)
    }
  }
}

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("Wrote %s\n  %d sheets, %d total significant DEG rows\n",
            OUT_XLSX, length(all_sig) + 2L, nrow(master)))
print(summary_wide)
