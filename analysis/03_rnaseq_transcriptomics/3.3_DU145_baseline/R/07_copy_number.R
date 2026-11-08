# RNA effects against lpWGS dosage differences between the 2x LOY population and sgCTRL.
# Saved QDNAseq segment values (T2T-CHM13v2.0) are used without resegmentation.
t2t_arms <- function(path) {
  cb <- data.table::fread(proj_path(path),header=FALSE,col.names=c('chrom','start','end','band','stain'))
  do.call(rbind,lapply(split(cb,cb$chrom),function(z) {
    p <- z[startsWith(z$band,'p')]; q <- z[startsWith(z$band,'q')]
    if (!nrow(p) || !nrow(q)) return(NULL)
    data.frame(chrom=z$chrom[1],p_end=max(p$end),q_start=min(q$start))
  }))
}
assign_arm <- function(chrom,mid,arms) {
  i <- match(chrom,arms$chrom)
  ifelse(is.na(i),NA_character_,ifelse(mid<=arms$p_end[i],paste0(chrom,'p'),ifelse(mid>=arms$q_start[i],paste0(chrom,'q'),NA_character_)))
}
du145_copy_number <- function(x,e,cfg) {
  a <- e$effects[e$effects$primary_eligible,]
  spans <- t2t_gene_spans(cfg$t2t_genes)
  mp <- map_to_t2t(cbind(a[,c('gene_id','gene_name','chrom','PAR')],stable_gene_id=a$stable_gene_id),spans,y_target=FALSE)
  g <- cbind(a[,c('gene_id','gene_name','chrom','arm','mean_CPM','mean_difference')],mp[,-1])
  bins <- load_wgs_bins(cfg)
  for (s in c(cfg$wgs_ctrl,cfg$wgs_loy,cfg$wgs_loy_alt)) {
    z <- gene_bin_cn(g,bins,s)
    lab <- c('sgCTRL','LOY2X','DMSO_LOY1X')[match(s,c(cfg$wgs_ctrl,cfg$wgs_loy,cfg$wgs_loy_alt))]
    g[[paste0('log2_seg_',lab)]] <- z$log2_seg; g[[paste0('n_bins_',lab)]] <- z$n_bins
  }
  g$CN_delta_LOY2X <- g$log2_seg_LOY2X-g$log2_seg_sgCTRL
  g$CN_delta_DMSO_LOY1X <- g$log2_seg_DMSO_LOY1X-g$log2_seg_sgCTRL
  g$CN_changed <- is.finite(g$CN_delta_LOY2X) & abs(g$CN_delta_LOY2X)>=cfg$cn_abs_delta
  arms <- t2t_arms(cfg$t2t_cytobands)
  g$t2t_mid <- (g$t2t_start+g$t2t_end)/2
  g$t2t_arm <- assign_arm(g$t2t_chrom,g$t2t_mid,arms)
  g$t2t_window_10Mb <- ifelse(is.na(g$t2t_mid),NA,paste0(g$t2t_chrom,':',floor(g$t2t_mid/1e7)*10,'Mb'))
  de <- e$de$table
  g <- merge(g,de[,c('gene_id','logFC','FDR','thesis_DE','pairs_agreeing_with_logFC')],by='gene_id',sort=FALSE)
  g$CN_concordant <- g$CN_changed & sign(g$CN_delta_LOY2X)==sign(g$logFC)
  g$CN_class <- ifelse(!is.finite(g$CN_delta_LOY2X),'no_usable_bins',ifelse(!g$CN_changed,'CN_stable',ifelse(g$CN_concordant,'CN_changed_concordant','CN_changed_discordant')))
  tab(g,'CN_gene_level')
  tab(g[g$thesis_DE,],'CN_thesis_DE_annotation')
  tab(data.frame(t2t_mapping=names(table(g$t2t_mapping)),genes=as.integer(table(g$t2t_mapping))),'CN_t2t_mapping_summary')
  # Arm and window level: RNA mean effect versus DNA bin-median difference.
  b <- bins[bins$use & is.finite(bins$log2_seg),]
  wide <- data.table::dcast(data.table::as.data.table(b),chromosome+start+end~sample,value.var='log2_seg')
  wide <- as.data.frame(wide)
  wide$delta <- wide[[cfg$wgs_loy]]-wide[[cfg$wgs_ctrl]]
  wide$delta_alt <- wide[[cfg$wgs_loy_alt]]-wide[[cfg$wgs_ctrl]]
  wide$mid <- (wide$start+wide$end)/2
  wide$arm <- assign_arm(wide$chromosome,wide$mid,arms)
  wide$window <- paste0(wide$chromosome,':',floor(wide$mid/1e7)*10,'Mb')
  unit_table <- function(key_g,key_b,label,min_genes) {
    rna <- aggregate(mean_difference~key,data.frame(key=g[[key_g]],mean_difference=g$mean_difference),mean)
    ng <- aggregate(n~key,data.frame(key=g[[key_g]],n=1),sum)
    dna <- aggregate(cbind(delta,delta_alt)~key,data.frame(key=wide[[key_b]],delta=wide$delta,delta_alt=wide$delta_alt),median)
    nb <- aggregate(n~key,data.frame(key=wide[[key_b]],n=1),sum)
    z <- merge(merge(merge(rna,ng,by='key'),dna,by='key'),nb,by='key',suffixes=c('_genes','_bins'))
    names(z) <- c('unit','RNA_mean_difference','genes','DNA_delta_LOY2X','DNA_delta_DMSO_LOY1X','bins')
    z$level <- label; z[z$genes>=min_genes,]
  }
  arm <- unit_table('t2t_arm','arm','arm',20); win <- unit_table('t2t_window_10Mb','window','window_10Mb',10)
  auto <- function(z) !grepl('^chrY',z$unit)
  tab(arm,'CN_arm_level_RNA_vs_DNA'); tab(win,'CN_window_10Mb_RNA_vs_DNA')
  gl <- g[is.finite(g$CN_delta_LOY2X),]
  conc <- function(xv,yv,lab,n) data.frame(level=lab,units=n,Pearson=cor(xv,yv),slope=unname(coef(lm(yv~xv))[2]))
  conc_tab <- rbind(conc(gl$CN_delta_LOY2X,gl$mean_difference,'gene',nrow(gl)),
                    conc(win$DNA_delta_LOY2X[auto(win)],win$RNA_mean_difference[auto(win)],'window_10Mb',sum(auto(win))),
                    conc(arm$DNA_delta_LOY2X[auto(arm)],arm$RNA_mean_difference[auto(arm)],'arm',sum(auto(arm))),
                    conc(arm$DNA_delta_DMSO_LOY1X[auto(arm)],arm$RNA_mean_difference[auto(arm)],'arm_vs_DMSO_LOY1X_lpWGS',sum(auto(arm))))
  tab(conc_tab,'CN_RNA_DNA_concordance')
  de_cn <- data.frame(set=c('thesis_DE','FDR<0.05'),
    genes=c(sum(g$thesis_DE),sum(g$FDR<cfg$de_fdr)),
    CN_stable=c(sum(g$thesis_DE & g$CN_class=='CN_stable'),sum(g$FDR<cfg$de_fdr & g$CN_class=='CN_stable')),
    CN_changed_concordant=c(sum(g$thesis_DE & g$CN_class=='CN_changed_concordant'),sum(g$FDR<cfg$de_fdr & g$CN_class=='CN_changed_concordant')),
    CN_changed_discordant=c(sum(g$thesis_DE & g$CN_class=='CN_changed_discordant'),sum(g$FDR<cfg$de_fdr & g$CN_class=='CN_changed_discordant')),
    no_usable_bins=c(sum(g$thesis_DE & g$CN_class=='no_usable_bins'),sum(g$FDR<cfg$de_fdr & g$CN_class=='no_usable_bins')))
  tab(de_cn,'CN_DE_summary')
  fdr_arm <- as.data.frame(table(arm=g$t2t_arm[g$FDR<cfg$de_fdr]))
  tab(fdr_arm[order(-fdr_arm$Freq),],'CN_FDR05_genes_by_arm')
  list(genes=g,arm=arm,window=win,concordance=conc_tab,de_cn=de_cn,bins=wide)
}
