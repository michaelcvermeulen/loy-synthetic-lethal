# ---------------------------------------------------------------------------
# read_xlsx_cells() - recover cell FILL COLOURS and FORMULAS from an .xlsx
#
# Needed because data/LOY_clones.xlsx carries two things readxl cannot see:
#
#   1. the Gel, DDX3X and Call columns record their result as a cell background
#      colour with no text in the cell, so those columns come back empty;
#   2. the ddPCR columns are formulas of the form "=147.9/307.5" - the two
#      ddPCR copy numbers (target, reference) behind each ratio. readxl returns
#      only the computed ratio, which throws away the counts the ratio's
#      uncertainty depends on.
#
# tidyxl would do both in one call but is not installed, and an .xlsx is just a
# zip of XML, so read it directly.
#
# Returns a tibble: sheet, cell, row, col, argb, formula
#   argb    - NA where the cell has no solid rgb fill
#   formula - NA where the cell holds a literal rather than a formula
# ---------------------------------------------------------------------------

library(xml2)
library(dplyr)

read_xlsx_cells <- function(path, sheet = NULL) {

  tmp <- file.path(tempdir(), paste0("xlsxcells_", basename(path)))
  unlink(tmp, recursive = TRUE)
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(path, exdir = tmp)

  strip <- function(f) xml_ns_strip(read_xml(file.path(tmp, f)))

  # --- style table: cell style index -> fill index -> ARGB -------------------
  sty <- strip("xl/styles.xml")
  fill_argb <- xml_find_all(sty, ".//fills/fill") %>%
    vapply(function(f) {
      fg <- xml_find_first(f, "./patternFill/fgColor")
      pt <- xml_attr(xml_find_first(f, "./patternFill"), "patternType")
      # Only solid rgb fills carry a result here. "none"/"gray125" are the
      # default empty styles, and theme colours are used only on the blank
      # spacer columns - both are left as NA rather than guessed at.
      if (is.na(pt) || pt != "solid" || inherits(fg, "xml_missing")) return(NA_character_)
      xml_attr(fg, "rgb")
    }, character(1))

  xf_fill <- xml_find_all(sty, ".//cellXfs/xf") %>%
    vapply(function(x) as.integer(xml_attr(x, "fillId")), integer(1))

  # --- sheet name -> worksheet part -----------------------------------------
  wb   <- strip("xl/workbook.xml")
  rels <- strip("xl/_rels/workbook.xml.rels")
  rel_target <- setNames(xml_attr(xml_find_all(rels, ".//Relationship"), "Target"),
                         xml_attr(xml_find_all(rels, ".//Relationship"), "Id"))

  sheet_nodes <- xml_find_all(wb, ".//sheets/sheet")
  sheet_name  <- xml_attr(sheet_nodes, "name")
  # the r:id attribute survives ns-stripping as plain "id"
  sheet_rid   <- xml_attr(sheet_nodes, "id")
  sheet_part  <- unname(rel_target[sheet_rid])
  sheet_part  <- ifelse(startsWith(sheet_part, "/"),
                        sub("^/", "", sheet_part),
                        file.path("xl", sheet_part))

  keep <- if (is.null(sheet)) seq_along(sheet_name) else which(sheet_name %in% sheet)
  stopifnot(length(keep) > 0)

  out <- lapply(keep, function(i) {
    ws <- strip(sheet_part[i])
    cells <- xml_find_all(ws, ".//sheetData/row/c")
    if (length(cells) == 0) return(NULL)
    ref <- xml_attr(cells, "r")
    s   <- suppressWarnings(as.integer(xml_attr(cells, "s")))
    s[is.na(s)] <- 0L
    tibble(
      sheet   = sheet_name[i],
      cell    = ref,
      row     = as.integer(sub("^[A-Z]+", "", ref)),
      col     = vapply(sub("[0-9]+$", "", ref), function(L) {
        # "A"->1, "AA"->27
        chars <- utils::head(strsplit(L, "")[[1]], 3)
        Reduce(function(a, b) a * 26L + (match(b, LETTERS)), chars, 0L)
      }, integer(1), USE.NAMES = FALSE),
      argb    = fill_argb[xf_fill[s + 1L] + 1L],
      formula = vapply(cells, function(cl) {
        f <- xml_find_first(cl, "./f")
        if (inherits(f, "xml_missing")) NA_character_ else xml_text(f)
      }, character(1))
    )
  })

  bind_rows(out)
}

# ---------------------------------------------------------------------------
# split_ratio_formula() - "=147.9/307.5" -> target 147.9, reference 307.5
#
# Every ddPCR formula in this workbook is a bare two-number division: the
# target copies/uL over the reference copies/uL. Anything else (a literal, a
# cell reference, a longer expression) returns NA rather than being guessed at,
# so a change of habit in the sheet shows up as missing counts instead of
# silently wrong ones.
# ---------------------------------------------------------------------------

split_ratio_formula <- function(x) {
  x  <- stringr::str_remove(stringr::str_trim(x), "^=")
  ok <- stringr::str_detect(x, "^[0-9.]+\\s*/\\s*[0-9.]+$") & !is.na(x)
  parts <- stringr::str_split_fixed(x, "\\s*/\\s*", 2)
  tibble(
    target_copies = ifelse(ok, suppressWarnings(as.numeric(parts[, 1])), NA_real_),
    ref_copies    = ifelse(ok, suppressWarnings(as.numeric(parts[, 2])), NA_real_)
  )
}
