# Plot-only rebuild from retained numerical results. No model or test is rerun.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(ggplot2); library(patchwork)})
roots <- c(getwd(), file.path(getwd(), 'THESIS_ADDBACK_2026'))
root <- roots[file.exists(file.path(roots, 'outputs/dge_primary.rds'))][1]
setwd(root)
STYLE <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style'
source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))
update_geom_defaults('text', list(family = 'Open Sans', colour = '#222222'))
group_levels <- c('CTRL', 'LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
cols <- c(CTRL = '#3C5488', LOY = '#E64B35', LOY_LXHA = '#7E6148',
          LOY_UTY_AB = '#00A087', LOY_KDM5D_AB = '#F39B7F', LOY_KDM5D_UTY_AB = '#4DBBD5')
labs_group <- c(CTRL = 'sgCTRL', LOY = 'Y-depleted', LOY_LXHA = 'Empty vector',
                LOY_UTY_AB = 'UTY', LOY_KDM5D_AB = 'KDM5D', LOY_KDM5D_UTY_AB = 'UTY + KDM5D')
save_multi <- function(p, stem, width, height) {
  ggsave(file.path('figures', paste0(stem, '.png')), p, width = width, height = height, dpi = 300, bg = 'white')
  ggsave(file.path('figures', paste0(stem, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path('figures', paste0(stem, '.svg')), p, width = width, height = height, device = svglite::svglite)
}
score_df <- fread('outputs/signature/frozen_signature_scores.csv', data.table = FALSE)
score_contrasts <- fread('outputs/signature/score_contrasts.csv', data.table = FALSE)
sig <- fread('outputs/signature/frozen_signature_mapping.csv', data.table = FALSE)
fr <- fread('outputs/signature/fry_directional.csv', data.table = FALSE)
fry_names <- c('UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA', 'DUAL_v_UTY', 'DUAL_v_KDM5D', 'LOY_v_LXHA')
dge <- readRDS('outputs/dge_primary.rds')
samples <- fread('outputs/sample_manifest.csv', data.table = FALSE)
stopifnot(identical(colnames(dge), samples$sample_id))
lcpm <- cpm(dge, log = TRUE, prior.count = 2)
sig_ids <- sig$gene_id_stable[sig$represented & sig$gene_id_stable %in% rownames(lcpm)]
gm <- sapply(group_levels, function(g) rowMeans(lcpm[sig_ids, samples$final_group == g, drop = FALSE]))
source('code/plot_signature_panels.R')
stopifnot(nrow(fr_plot) == 12, identical(fr_plot$FDR_12, fr$FDR_12), identical(fr_plot$Direction, fr$Direction))
files <- c('figures/ch03_main_ht1080_addback_baseline_reversal.png',
           'figures/ch03_main_ht1080_addback_baseline_reversal.pdf',
           'figures/ch03_main_ht1080_addback_baseline_reversal.svg',
           list.files('figures/student_panels', pattern = '^panel_[fghi]_.*[.](png|pdf|svg)$', full.names = TRUE))
stopifnot(all(file.copy(files, file.path('student_report/figures', basename(files)), overwrite = TRUE)))
file.copy('figures/ch03_main_ht1080_addback_baseline_reversal.png', 'report/figures/ch03_main_ht1080_addback_baseline_reversal.png', overwrite = TRUE)
file.copy('outputs/signature/fry_panel_display.csv', 'student_report/data/outputs__signature__fry_panel_display.csv', overwrite = TRUE)
scripts <- c('code/run_analysis.R', 'code/plot_signature_panels.R', 'code/restyle_signature_figure.R')
file.copy(scripts, file.path('student_report/code', basename(scripts)), overwrite = TRUE)
writeLines(c('PASS: 12 retained fry results, directions and FDR values unchanged.',
             'Open Sans Regular/SemiBold from the thesis figure skill; SVG exported by svglite, PDF printed from SVG with embedded fonts using headless Chrome.',
             'Panel D: all FDR values shown; blue=decrease, orange=increase; bold indicates FDR<0.05.'),
           'provenance/signature_restyle_verification.txt')

