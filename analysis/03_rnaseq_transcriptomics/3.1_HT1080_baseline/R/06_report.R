finalize_run <- function(x, e, p, d, cfg, repo) {
  message('Writing provenance, validation results and report')
  # Older symbol-only exports cannot safely be joined at ambiguous loci.
  old_path <- file.path(repo,'analysis/HT1080/reannotate_v105/HT1080_DE_with_coords_paired.csv')
  old <- as.data.frame(data.table::fread(old_path))
  old$symbol_unambiguous <- !(duplicated(old$gene)|duplicated(old$gene,fromLast=TRUE))
  a <- e$effects[e$effects$primary_eligible, ]
  a$symbol_unambiguous <- !(duplicated(a$gene_name)|duplicated(a$gene_name,fromLast=TRUE))
  cmp <- merge(a[a$symbol_unambiguous,c('gene_id','gene_name','chrom','AUG_2025','mean_difference')],
                old[old$symbol_unambiguous,c('gene','chr','logFC')],by.x='gene_name',by.y='gene')
  cmp <- cmp[cmp$chrom == cmp$chr, ]; names(cmp)[names(cmp)=='logFC'] <- 'old_paired_logFC'
  cmp$comparison_scope <- 'unambiguous_symbol_and_chromosome_only_old_export_lacks_exact_gene_id'
  tab(cmp,'legacy_August_effect_comparison')
  primary <- e$effects[e$effects$primary_eligible, ]
  rr <- do.call(rbind,lapply(split(d$regional,paste(d$regional$batch,d$regional$region)),function(z) {
    data.frame(batch=z$batch[1],region=z$region[1],
      LOY_mean_CPM=mean(z$sum_CPM[z$treatment=='LOY']),CTRL_mean_CPM=mean(z$sum_CPM[z$treatment=='CTRL']))
  }))
  rr$ratio_of_period_mean_CPM <- ifelse(rr$CTRL_mean_CPM>0,rr$LOY_mean_CPM/rr$CTRL_mean_CPM,NA_real_)
  tab(rr,'Y_PAR_MT_period_summary')
  mainfg <- p$fgsea[p$fgsea$contrast=='combined', ]
  summary <- list(raw_genes=nrow(x$counts),libraries=ncol(x$counts),core_libraries=nrow(x$core),
    primary_genes=nrow(primary),unique_symbols=nrow(e$mapping),candidates=sum(primary$descriptive_candidate),
    period_pearson=cor(primary$AUG_2025,primary$MAY_2026),
    period_spearman=cor(primary$AUG_2025,primary$MAY_2026,method='spearman'),
    period_direction_fraction=mean(primary$period_direction_agrees),
    RLE_correlation=cor(primary$mean_difference,primary$RLE_mean_difference),
    prior_correlation=cor(primary$mean_difference,primary$prior025_mean_difference),
    median_LOO_change=median(primary$LOO_max_abs_change),
    genes_pass_total50=sum(primary$eligible_min_total_50),genes_pass_total100=sum(primary$eligible_min_total_100),
    assignment_min=min(e$qc$assigned_fraction),assignment_max=max(e$qc$assigned_fraction),
    old_comparable_genes=nrow(cmp),old_correlation=cor(cmp$AUG_2025,cmp$old_paired_logFC),
    pathways_tested=nrow(mainfg),pathways_BH05=sum(mainfg$BH_across_collections<0.05,na.rm=TRUE),
    msigdb_resource_sets=nrow(p$set_coverage),msigdb_collections=length(unique(p$set_coverage$collection)),
    camera_tested=sum(p$camera$contrast=='combined'),
    camera_BH05=sum(p$camera$contrast=='combined' & p$camera$BH_across_collections<0.05,na.rm=TRUE),
    camera_sensitivity_BH05=sum(p$camera_sensitivity$contrast=='combined' & p$camera_sensitivity$BH_across_collections<0.05,na.rm=TRUE),
    TFs_tested=sum(p$tf$contrast=='combined'),
    msigdb_version=unique(data.table::fread('resources/msigdb_all_set_metadata.csv')$db_version),
    generated_at=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'))
  saveRDS(summary,'outputs/models/run_summary.rds')
  yaml::write_yaml(summary,'provenance/run_summary.yml')
  pkgs <- c('edgeR','limma','data.table','rtracklayer','SummarizedExperiment','GenomicRanges','S4Vectors',
    'fgsea','BiocParallel','msigdbr','decoupleR','progeny','ggplot2','patchwork','ggrepel','scales',
    'ragg','svglite','colorspace','systemfonts','digest','yaml','rmarkdown','knitr','renv','testthat','png')
  tab(data.frame(package=pkgs,version=vapply(pkgs,function(z)as.character(packageVersion(z)),character(1))), 'package_versions')
  capture.output(sessionInfo(),file='provenance/sessionInfo.txt')
  sink('provenance/scientific_checks.txt')
  tryCatch({
    source('tests/check_scientific_invariants.R')
    source('tests/check_hdac12_invariants.R')
    source('tests/check_hdac_external_invariants.R')
    source('tests/check_epigenetic_invariants.R')
  }, finally = sink())
  # Include all actual raw tables and assignment summaries, not only the assembled matrix.
  paths <- unique(c(file.path(repo,cfg$raw_counts),file.path(repo,cfg$gtf),x$meta$source_count_file,
    paste0(x$meta$source_count_file,'.summary'),cfg$cn_source,old_path,
    'A:/Dropbox/LOY_paralogs/projects/proposals_summaries/deep_research/RNAseq_methods_review.md',
    list.files('config',full.names=TRUE),list.files('R',full.names=TRUE),
    list.files('resources',full.names=TRUE,recursive=TRUE),list.files('tests',full.names=TRUE),
    'run.R','run_hdac12_step2.R','run_hdac12_step3.R','prepare_resources.R','prepare_chromatin_resources.R',
    'audit_hdac12_sources.py','prepare_hdac12_external.py','audit_hdac12_external.py','render_report.R','review_figures.R',
    'run_epigenetic_deep_dive.R','prepare_epigenetic_resources.R','audit_epigenetic_sources.py','EPIGENETIC_DEEP_DIVE_PLAN.md',
    'run_chromatin_GSEA_plots.R','finalize_chromatin_GSEA_plots.R',
    list.files('report',pattern='\\.Rmd$',full.names=TRUE),
    'README.md','PLAN.md','HDAC12_PLAN.md','FIGURES.md','provenance/DECISIONS.md','provenance/ZMYND8_interpretation.md',
    'provenance/chromatin_interpretation.md','provenance/HDAC12_step2_interpretation.md','provenance/HDACEXT_interpretation.md',
    'FIGURE_CAPTIONS.md'))
  assert(all(file.exists(paths)),'Missing declared provenance dependency')
  provenance <- data.frame(path=paths,sha256=vapply(paths,hash_file,character(1)),
                           bytes=file.info(paths)$size,recorded_at=summary$generated_at)
  data.table::fwrite(provenance,'provenance/input_code_resource_hashes.csv')
  # Snapshot records the existing compatible library; it does not update installed packages.
  needs_snapshot <- !file.exists('renv.lock')
  if (!needs_snapshot) {
    locked <- jsonlite::read_json('renv.lock')$Packages
    same <- vapply(locked,function(z) identical(packageDescription(z$Package)$Version,z$Version),logical(1))
    needs_snapshot <- !all(same) || !all(pkgs %in% names(locked))
  }
  if (needs_snapshot) renv::snapshot(project=getwd(),lockfile='renv.lock',packages=pkgs,prompt=FALSE,force=TRUE)
  provenance <- rbind(provenance,data.frame(path='renv.lock',sha256=hash_file('renv.lock'),
    bytes=file.info('renv.lock')$size,recorded_at=summary$generated_at))
  data.table::fwrite(provenance,'provenance/input_code_resource_hashes.csv')
  rmarkdown::pandoc_convert(normalizePath('outputs/captions/figure_captions.md',winslash='/'),to='docx',
    output=file.path(normalizePath('outputs/captions',winslash='/'),'figure_captions.docx'))
  rmarkdown::render('report/analysis.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),
                     envir=new.env(parent=globalenv()),quiet=TRUE)
  rmarkdown::render('report/HDAC12_gene_drivers.Rmd', output_dir=normalizePath('outputs/report',winslash='/'),
                     envir=new.env(parent=globalenv()), quiet=TRUE)
  rmarkdown::render('report/HDAC12_external_comparison.Rmd', output_dir=normalizePath('outputs/report',winslash='/'),
                     envir=new.env(parent=globalenv()), quiet=TRUE)
  rmarkdown::render('report/epigenetic_deep_dive.Rmd', output_dir=normalizePath('outputs/report',winslash='/'),
                     envir=new.env(parent=globalenv()), quiet=TRUE)
  rmarkdown::render('report/chromatin_GSEA_curves.Rmd', output_dir=normalizePath('outputs/report',winslash='/'),
                     envir=new.env(parent=globalenv()), quiet=TRUE)
  message('Core analysis and report complete. Figure visual review is recorded separately.')
  invisible(summary)
}
