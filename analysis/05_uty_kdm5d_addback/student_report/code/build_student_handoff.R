# Refresh the exploratory analysis and portable report using existing CAMERA outputs.
roots <- c(getwd(), file.path(getwd(), 'THESIS_ADDBACK_2026'))
root <- normalizePath(roots[file.exists(file.path(roots, 'outputs/fit_primary.rds'))][1], winslash = '/')
setwd(root)
source('code/run_responsive_genes.R', local = new.env(parent = globalenv()))
for (d in c('student_report/data', 'student_report/code', 'student_report/figures', 'report/figures'))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
sources <- list.files('outputs/responsive_genes', full.names = TRUE)
stopifnot(all(file.copy(sources, file.path('student_report/data', gsub('[/\\\\]', '__', sources)), overwrite = TRUE)))
figures <- list.files('figures', pattern = '^student_responsive_.*\\.(png|pdf|svg)$', full.names = TRUE)
stopifnot(length(figures) == 24,
  all(file.copy(figures, file.path('student_report/figures', basename(figures)), overwrite = TRUE)),
  all(file.copy(figures, file.path('report/figures', basename(figures)), overwrite = TRUE)))
file.copy('figures/student_figure_index.csv', 'student_report/figures/student_figure_index.csv', overwrite = TRUE)
signature_files <- c(list.files('figures', pattern = '^ch03_main_ht1080_addback_baseline_reversal[.](png|pdf|svg)$', full.names = TRUE), list.files('figures/student_panels', pattern = '^panel_[fghi]_.*[.](png|pdf|svg)$', full.names = TRUE))
stopifnot(all(file.copy(signature_files, file.path('student_report/figures', basename(signature_files)), overwrite = TRUE)))
file.copy('figures/ch03_main_ht1080_addback_baseline_reversal.png', 'report/figures/ch03_main_ht1080_addback_baseline_reversal.png', overwrite = TRUE)
scripts <- c('code/export_signature_pdfs.py', 'code/plot_signature_panels.R', 'code/restyle_signature_figure.R', 'code/run_analysis.R', 'code/run_responsive_genes.R', 'code/build_student_handoff.R', 'code/run_camera_report.R', 'report/student_report.Rmd')
stopifnot(all(file.copy(scripts, file.path('student_report/code', basename(scripts)), overwrite = TRUE)))
rmarkdown::render('report/student_report.Rmd', output_file = 'HT1080_addback_student_report_2026-09-21.html',
  output_dir = 'student_report', knit_root_dir = root, envir = new.env(parent = globalenv()), quiet = TRUE)
readme <- if (file.exists('student_report/README.md')) readLines('student_report/README.md', warn = FALSE) else
  c('# HT1080 add-back student handoff', '', 'Open HT1080_addback_student_report_2026-09-21.html in a browser.')
note <- 'Exploratory response gene lists are in data/outputs__responsive_genes__*.csv. The report includes selection rules, TREAT sensitivity and three standalone figures. Run code/build_student_handoff.R from the original analysis directory to refresh these outputs and the ZIP using existing CAMERA results.'
if (!note %in% readme) writeLines(c(readme, '', note), 'student_report/README.md')
finalize_handoff <- function() {
files <- list.files('student_report', recursive = TRUE, full.names = TRUE)
files <- files[!grepl('FILE_MANIFEST.csv$', files)]
data.table::fwrite(data.frame(file = sub('^student_report/', '', files), bytes = file.size(files),
  sha256 = vapply(files, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'student_report/FILE_MANIFEST.csv')
zipfile <- file.path(root, 'HT1080_addback_student_handoff_2026-09-21.zip')
zip::zipr(zipfile, files = 'student_report', root = root, include_directories = TRUE)
stopifnot(file.size(zipfile) > 0)
}
finalize_handoff()
cat('Rendered HTML and rebuilt student ZIP.\n')
