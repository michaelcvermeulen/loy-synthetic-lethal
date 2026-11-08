par_branch<-function(x,e,cfg) {
  a<-e$effects
  a$PAR_region<-NA_character_;a$PAR_boundary_overlap<-FALSE
  bed<-data.table::fread('config/hg38_PAR.bed',header=FALSE)
  for(j in seq_len(nrow(bed))) {
    ix<-which(a$chrom==bed[[1]][j]&a$start<=bed[[3]][j]&a$end>bed[[2]][j])
    a$PAR_region[ix]<-bed[[4]][j]
    a$PAR_boundary_overlap[ix]<-a$start[ix]<=bed[[2]][j]|a$end[ix]>bed[[3]][j]
  }
  # Secondary inclusive fit estimates the voom trend over all expressed nuclear genes.
  # The primary non-PAR analysis, its normalization and its test family are unchanged.
  k<-a$expressed & a$chrom%in%c(paste0('chr',1:22),'chrX','chrY')
  anno<-a[k,c('gene_id','gene_name','chrom','start','end','PAR','PAR_region','PAR_boundary_overlap','mean_difference','CTRL_mean_CPM','LOY_mean_CPM','LOO_all_same_direction')]
  fit<-fit_unpaired(x$counts[k,x$core$sample_key,drop=FALSE],anno,e$norm$effective,x$core)
  all<-fit$table;names(all)[names(all)=='FDR']<-'BH_all_expressed_nuclear_genes'
  tab(all,'DE_secondary_inclusive_nuclear_voom')
  par<-all[all$PAR,];par$BH_PAR_family<-p.adjust(par$PValue,'BH')
  par$CPM_percent_change<-100*(par$LOY_mean_CPM/par$CTRL_mean_CPM-1)
  par$analysis_role<-'secondary_PAR_followup_requested_after_primary_results;unpaired_exploratory'
  tab(par,'DE_PAR_secondary_voom')
  tab(a[a$PAR,],'PAR_all_gene_descriptive_effects')
  saveRDS(fit,'outputs/models/secondary_inclusive_nuclear_voom.rds')
  gene_ids<-par$gene_id[par$chrom=='chrX'&par$gene_name%in%c('PLCXD1','GTPBP6','PPP2R3B','CRLF2','CSF2RA','SLC25A6','ASMTL','AKAP17A','DHRSX','ZBED1','CD99','SPRY3','VAMP7')]
  abundance<-do.call(rbind,lapply(gene_ids,function(g) {
    i<-match(g,a$gene_id)
    data.frame(gene_id=g,gene_name=a$gene_name[i],PAR_region=a$PAR_region[i],x$core[,c('sample_key','treatment','sample_index')],CPM=as.numeric(e$norm$CPM[g,]))
  }))
  tab(abundance,'PAR_selected_gene_sample_abundance')
  suppressPackageStartupMessages(library(ggplot2))
  p<-ggplot(abundance,aes(treatment,CPM,color=treatment))+geom_point(position=position_jitter(width=.08,seed=42),size=1.7)+facet_wrap(~PAR_region+gene_name,ncol=4,scales='free_y')+scale_color_manual(values=COLORS)+theme_bw(base_size=11)+labs(title='A549 pseudoautosomal transcript abundance',x=NULL,y='CPM using primary TMM scaling')
  save_figure(p,'Fig09_PAR_expression',11,9,'Selected expressed PAR1 and PAR2 genes; all eight libraries shown using the primary normalization. Chromosome-X feature coordinates do not establish the chromosome of origin for homologous PAR transcripts. The complete PAR table contains every expressed PAR-overlapping feature, identifies boundary-spanning genes, and reports both PAR-family and inclusive nuclear BH. This secondary unpaired analysis remains conditional on uncertain derivation history.')
  invisible(par)
}
