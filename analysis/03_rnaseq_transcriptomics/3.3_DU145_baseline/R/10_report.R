html_escape <- function(s) {s<-gsub('&','&amp;',as.character(s),fixed=TRUE);s<-gsub('<','&lt;',s,fixed=TRUE);gsub('>','&gt;',s,fixed=TRUE)}
html_table <- function(d) {
  d <- as.data.frame(d)
  for (j in seq_along(d)) if (is.numeric(d[[j]])) d[[j]] <- format(signif(d[[j]],4),trim=TRUE)
  paste0('<div class="table"><table><thead><tr>',paste0('<th>',html_escape(names(d)),'</th>',collapse=''),'</tr></thead><tbody>',
         paste(apply(d,1,function(row)paste0('<tr>',paste0('<td>',html_escape(row),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table></div>')
}
write_report <- function(x,e,id,y,cn,p,cfg) {
  t <- e$de$table; ds <- e$de_summary; yt <- y$table; cc <- cn$concordance
  du <- id$ranked[id$ranked$StrippedCellLineName=='DU145' & id$ranked$profile %in% c('sgCTRL_mean','LOY_mean'),]
  ln <- id$ranked[id$ranked$StrippedCellLineName=='LNCAPCLONEFGC' & id$ranked$profile %in% c('sgCTRL_mean','LOY_mean'),]
  pct <- function(g) yt$percent_of_sgCTRL[match(g,yt$gene_name)]
  dep <- yt[yt$expressed & yt$lpWGS_interval %in% 'depleted MSY' & yt$CTRL_mean_CPM>=5,]
  summary <- list(generated_at=format(Sys.time(),'%Y-%m-%d %H:%M:%S %z'),libraries=ncol(x$counts),primary_genes=sum(e$effects$primary_eligible),
    PC1_percent=100*e$pca_variance[1],PC2_percent=100*e$pca_variance[2],PC1_passage_R2=e$pc_assoc$passage_R2[1],PC2_treatment_R2=e$pc_assoc$treatment_R2[2],
    paired_FDR05=ds$genes[1],paired_FDR05_higher=ds$higher[1],paired_FDR05_lower=ds$lower[1],FDR05_LFC1=ds$genes[2],thesis_DE=ds$genes[3],thesis_DE_higher=ds$higher[3],thesis_DE_lower=ds$lower[3],unpaired_FDR05=ds$genes[4],
    DU145_rank_sgCTRL=du$rank[du$profile=='sgCTRL_mean'],DU145_rho_sgCTRL=du$spearman[du$profile=='sgCTRL_mean'],LNCaP_rank_sgCTRL=ln$rank[ln$profile=='sgCTRL_mean'],
    RPS4Y1_percent=pct('RPS4Y1'),ZFY_percent=pct('ZFY'),depleted_MSY_percent_min=min(dep$percent_of_sgCTRL),depleted_MSY_percent_max=max(dep$percent_of_sgCTRL),
    CN_r_gene=cc$Pearson[cc$level=='gene'],CN_r_window=cc$Pearson[cc$level=='window_10Mb'],CN_r_arm=cc$Pearson[cc$level=='arm'],
    thesis_DE_CN_stable=cn$de_cn$CN_stable[1],FDR05_CN_concordant=cn$de_cn$CN_changed_concordant[2],enrichment_complete=!is.null(p))
  if (!is.null(p)) summary <- c(summary,list(pathways_tested=nrow(p$fgsea),fgsea_global_BH05=sum(p$fgsea$BH_across_collections<.05,na.rm=TRUE),
    fgsea_survive_CN=sum(p$fgsea$survives_CN_screen,na.rm=TRUE),fgsea_interpretable=sum(p$fgsea$interpretable,na.rm=TRUE),camera_global_BH05=sum(p$camera$BH_across_collections<.05,na.rm=TRUE),
    camera_estimated_BH05=sum(p$camera_sensitivity$BH_across_collections<.05,na.rm=TRUE)))
  jsonlite::write_json(summary,'provenance/run_summary.json',pretty=TRUE,auto_unbox=TRUE)
  s <- summary
  lines <- c('# DU145 LOY versus sgCTRL: thesis reanalysis','',paste('Generated:',s$generated_at),'',
    '## Design',sprintf('Six libraries from the May 2026 sampling window: three sgCTRL and three LOY at matched passages P12, P13 and P14. The LOY population is the CP + Q53 + P47 2x population (author, 23 September 2026). Source records label the line LNCaP; the RNA identifies it as DU145. Passages are repeated sampling of one population per arm, not independent editing events. The primary model is paired voom-limma, ~ passage + condition (residual df 2), with robust empirical Bayes. PC1 (%.1f%%) tracks passage (R2 = %.2f) and PC2 (%.1f%%) tracks condition (R2 = %.2f).',s$PC1_percent,s$PC1_passage_R2,s$PC2_percent,s$PC2_treatment_R2),'',
    '## Identity',sprintf('Against %s default DepMap 26Q1 profiles (8,000 most variable shared genes, Spearman), the sgCTRL mean ranks DU145 first (rho = %.3f); LNCaP clone FGC ranks %s. AR, KLK3, KLK2, FOLH1 and HOXB13 have zero reads in all six libraries. STR profiling was not performed in this analysis.',format(length(unique(id$ranked$ModelID)),big.mark=','),s$DU145_rho_sgCTRL,s$LNCaP_rank_sgCTRL),'',
    '## Chromosome Y',sprintf('Expression follows the lpWGS segment structure. The two expressed MSY genes in the retained distal Yp interval stay at %.0f%% (RPS4Y1) and %.0f%% (ZFY) of sgCTRL; expressed genes in the depleted MSY interval (sgCTRL CPM >= 5) are at %.1f-%.1f%%. PAR RNA is chrX-assigned and cannot be attributed to X or Y.',s$RPS4Y1_percent,s$ZFY_percent,s$depleted_MSY_percent_min,s$depleted_MSY_percent_max),'',
    '## Autosomal and X expression',sprintf('Across %s primary genes, %s pass FDR < 0.05 (%s higher, %s lower); %s also have |log2FC| >= 1, and %s meet the thesis rule (same direction in all three passage pairs; %s higher, %s lower). The unpaired sensitivity model gives %s at FDR < 0.05.',format(s$primary_genes,big.mark=','),s$paired_FDR05,s$paired_FDR05_higher,s$paired_FDR05_lower,s$FDR05_LFC1,s$thesis_DE,s$thesis_DE_higher,s$thesis_DE_lower,s$unpaired_FDR05),'',
    '## Copy number',sprintf('RNA effects track lpWGS LOY-2X minus sgCTRL segment differences increasingly with genomic scale: Pearson r = %.2f (gene), %.2f (10-Mb window), %.2f (arm). %s of the FDR < 0.05 genes are copy-number-concordant; %s of %s thesis DE genes are copy-number-stable.',s$CN_r_gene,s$CN_r_window,s$CN_r_arm,s$FDR05_CN_concordant,s$thesis_DE_CN_stable,s$thesis_DE),'')
  if (!is.null(p)) lines <- c(lines,'## Pathways',sprintf('Of %s eligible MSigDB 2026.1.Hs sets, %s pass global BH < 0.05 by effect-ranked fgsea and %s also pass after removing copy-number-changed genes, and %s of those are not positionally concentrated (C1 sets or >= 50%% of the leading edge on one arm). Fixed-correlation CAMERA gives %s; estimated-correlation CAMERA gives %s.',format(s$pathways_tested,big.mark=','),s$fgsea_global_BH05,s$fgsea_survive_CN,s$fgsea_interpretable,s$camera_global_BH05,s$camera_estimated_BH05),'')
  lines <- c(lines,'Reproduce with R 4.5.2 from this directory: `Rscript --vanilla run.R` (`--core-only`, `--enrichment-only`, `--report-only`). Input, resource, DepMap and code hashes and validation are in `provenance/` and `outputs/tables/`.')
  writeLines(lines,'outputs/report/analysis.md')
  body <- c('<h1>DU145 LOY versus sgCTRL</h1>')
  for (l in lines[-(1:3)]) body <- c(body,if (grepl('^## ',l)) paste0('<h2>',html_escape(sub('^## ','',l)),'</h2>') else if (nzchar(l)) paste0('<p>',html_escape(l),'</p>'))
  body <- c(body,'<h2>DE rule summary</h2>',html_table(ds))
  body <- c(body,'<h2>Thesis DE genes</h2>',html_table(merge(t[t$thesis_DE,c('gene_id','gene_name','chrom','logFC','FDR','pairs_agreeing_with_logFC')],cn$genes[,c('gene_id','t2t_arm','CN_delta_LOY2X','CN_class')],by='gene_id')[,-1]))
  body <- c(body,'<h2>Chromosome Y and PAR</h2>',html_table(yt[yt$expressed,c('gene_name','region','t2t_mid_Mb','lpWGS_interval','CTRL_mean_CPM','percent_of_sgCTRL','logFC','BH_Y_PAR_family','lpWGS_LOY2X_over_sgCTRL')]))
  body <- c(body,'<h2>X-Y paralogs</h2>',html_table(y$pairs[,c('X','Y','X_mean_difference','X_BH_inclusive','Y_CTRL_CPM','Y_mean_difference','summed_XY_percent_of_sgCTRL')]))
  body <- c(body,'<h2>Copy-number concordance</h2>',html_table(cc),html_table(cn$de_cn))
  if (!is.null(p)) {
    top <- p$fgsea[p$fgsea$BH_across_collections<.05,]; top <- head(top[order(top$pval),c('pathway','collection','NES','BH_across_collections','NES_CN_stable_only','survives_CN_screen','positionally_concentrated','interpretable','leading_edge_top_arm','leading_edge_top_arm_fraction')],40)
    body <- c(body,'<h2>Enrichment</h2>',html_table(p$summary),html_table(top))
  }
  for (f in list.files('outputs/figures',pattern='\\.png$',full.names=TRUE)) {
    stem <- sub('\\.png$','',basename(f)); cap <- readLines(paste0('outputs/captions/',stem,'.md'),warn=FALSE)
    body <- c(body,paste0('<h2>',html_escape(gsub('_',' ',stem)),'</h2><img src="data:image/png;base64,',base64enc::base64encode(f),'" alt="',html_escape(stem),'"/><p>',html_escape(paste(cap[-c(1,2)],collapse=' ')),'</p>'))
  }
  html <- c('<!doctype html><html><head><meta charset="utf-8"><title>DU145 RNA-seq reanalysis</title><style>body{max-width:1100px;margin:40px auto;padding:0 24px;font-family:Arial,sans-serif;line-height:1.55;color:#22303b}h1,h2{color:#175b80}img{max-width:100%;height:auto}.table{overflow-x:auto}table{border-collapse:collapse;font-size:13px}th,td{padding:6px;border-bottom:1px solid #ddd;text-align:left}th{background:#edf3f7}</style></head><body>',body,'</body></html>')
  writeLines(html,'outputs/report/analysis.html')
  invisible(summary)
}
