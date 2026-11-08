# Stage 00: fail fast if the environment cannot do what later stages need.
# - snapshot the msigdbr gene sets into derived/ so jobs do not depend on the user cache
# - check that figures render with the house style in PDF, SVG and PNG
# - check that an Rmd renders to HTML and to PDF
# - check that an .rds written here reads back correctly in the main (forest) environment

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "plot_style.R"))

reference_dir <- file.path(derived_dir, "reference")
smoke_dir <- file.path(results_dir, "smoke_test")
dir.create(reference_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(smoke_dir, showWarnings = FALSE, recursive = TRUE)

# --- gene sets ---------------------------------------------------------------
collections <- list(
  hallmark = list(collection = "H", subcollection = NULL),
  reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  gobp = list(collection = "C5", subcollection = "GO:BP"),
  positional = list(collection = "C1", subcollection = NULL)
)
gene_sets <- lapply(collections, function(spec) {
  tab <- as.data.table(suppressWarnings(msigdbr::msigdbr(
    species = "Homo sapiens", collection = spec$collection, subcollection = spec$subcollection
  )))
  split(tab$gene_symbol, tab$gs_name)
})
set_counts <- vapply(gene_sets, length, integer(1))
log_msg("gene sets: %s", paste(names(set_counts), set_counts, sep = "=", collapse = ", "))
stopifnot(set_counts[["hallmark"]] == 50, all(set_counts > 0))
saveRDS(
  list(sets = gene_sets, msigdbr_version = as.character(utils::packageVersion("msigdbr")), created = Sys.time()),
  file.path(reference_dir, "gene_sets.rds")
)

# --- figure rendering --------------------------------------------------------
toy <- data.frame(
  group = factor(rep(c("Y-positive", "Y-negative"), each = 20), levels = c("Y-positive", "Y-negative")),
  value = c(rnorm(20), rnorm(20, 1))
)
figure <- ggplot2::ggplot(toy, ggplot2::aes(group, value, colour = group)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.15, size = 0.8) +
  scale_colour_y_status() +
  theme_house()
paths <- save_figure(figure, "smoke_test", width = 3, height = 2.5, dir = smoke_dir)
stopifnot(all(file.exists(paths)), all(file.size(paths) > 0))
log_msg("figure rendering ok: %s", paste(basename(paths), collapse = ", "))

# --- report rendering --------------------------------------------------------
# Named differently from the figure: LaTeX resolves an image without its
# extension and would otherwise pick up the PDF it is busy writing.
rmd <- file.path(smoke_dir, "render_check.Rmd")
writeLines(c(
  "---", "title: \"Smoke test\"", "output:",
  "  html_document:", "    self_contained: true",
  "  pdf_document:", "    latex_engine: xelatex",
  "---", "", "A table and a figure.", "",
  "```{r, echo=FALSE}", "knitr::kable(head(mtcars, 3))", "```", "",
  "```{r, echo=FALSE, out.width='50%'}", "knitr::include_graphics('smoke_test.png')", "```"
), rmd)
rendered <- rmarkdown::render(rmd, output_format = "all", quiet = TRUE)
stopifnot(all(file.exists(rendered)), any(grepl("\\.pdf$", rendered)), any(grepl("\\.html$", rendered)))
log_msg("report rendering ok: %s", paste(basename(rendered), collapse = ", "))

# --- cross-version .rds ------------------------------------------------------
probe <- list(
  matrix = matrix(c(1.5, NA, 3, 4), 2, 2, dimnames = list(c("ACH-000001", "ACH-000002"), c("A (1)", "B (2)"))),
  frame = data.frame(ModelID = c("ACH-000001", "ACH-000002"), lineage = c("Lung", "Skin"), stringsAsFactors = FALSE)
)
probe_path <- file.path(smoke_dir, "probe.rds")
saveRDS(probe, probe_path)
main_rscript <- "/global/home/hpc5173/miniconda3/envs/main/bin/Rscript"
check <- system2(main_rscript, c("--vanilla", "-e", shQuote(sprintf(
  "p <- readRDS('%s'); stopifnot(is.matrix(p$matrix), is.na(p$matrix[2, 1]), p$matrix[1, 2] == 3, identical(rownames(p$matrix), c('ACH-000001', 'ACH-000002')), nrow(p$frame) == 2); cat('rds ok in', R.version.string)",
  probe_path
))), stdout = TRUE, stderr = TRUE)
log_msg("cross-version rds: %s", paste(check, collapse = " "))
stopifnot(any(grepl("rds ok", check)))

writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info_scRNA.txt"))
log_msg("smoke test passed")
