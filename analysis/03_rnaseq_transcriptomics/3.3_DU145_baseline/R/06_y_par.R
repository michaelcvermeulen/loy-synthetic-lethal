proj_path <- function(p) file.path(dirname(REPO),p)
# Shared genome helpers: GENCODE v49 (GRCh38) genes -> CHM13v2.0 spans, and saved QDNAseq bins.
t2t_gene_spans <- function(path) {
  g <- utils::read.delim(gzfile(proj_path(path)),header=TRUE,quote='',comment.char='',colClasses='character')[,c(1,2,3,13,21)]
  names(g) <- c('chrom','start','end','symbol','source_gene')
  g <- data.table::as.data.table(g); g$start <- as.integer(g$start); g$end <- as.integer(g$end)
  g <- g[!is.na(g$source_gene) & grepl('^ENSG',g$source_gene)]
  g$stable_gene_id <- sub('\\..*$','',g$source_gene)
  s <- g[,.(t2t_start=min(start),t2t_end=max(end),symbol=symbol[1]),by=.(chrom,stable_gene_id)]
  as.data.frame(s)
}
map_to_t2t <- function(a,spans,y_target=FALSE) {
  # PAR genes are chrX-assigned in the RNA; for chromosome-Y display map them to CHM13 chrY PAR spans.
  target <- ifelse(y_target & a$PAR,'chrY',a$chrom)
  k1 <- paste(target,a$stable_gene_id); k2 <- paste(spans$chrom,spans$stable_gene_id)
  i <- match(k1,k2); method <- ifelse(is.na(i),NA_character_,'stable_ENSG_same_chrom')
  sym_key <- paste(spans$chrom,spans$symbol)
  uniq <- !(duplicated(sym_key)|duplicated(sym_key,fromLast=TRUE))
  j <- match(paste(target,a$gene_name),sym_key[uniq])
  fb <- is.na(i) & !is.na(j) & !is.na(a$gene_name) & nzchar(a$gene_name)
  i[fb] <- which(uniq)[j[fb]]; method[fb] <- 'symbol_same_chrom_unique'
  data.frame(gene_id=a$gene_id,t2t_chrom=ifelse(is.na(i),NA,spans$chrom[i]),t2t_start=spans$t2t_start[i],
             t2t_end=spans$t2t_end[i],t2t_mapping=ifelse(is.na(method),'unmapped',method))
}
load_wgs_bins <- function(cfg) {
  d <- as.data.frame(readRDS(proj_path(cfg$wgs_bins)))
  d <- d[d$sample %in% c(cfg$wgs_ctrl,cfg$wgs_loy,cfg$wgs_loy_alt),c('sample','chromosome','start','end','log2','log2_seg','use')]
  d
}
gene_bin_cn <- function(map,bins,sample) {
  b <- bins[bins$sample==sample & bins$use & is.finite(bins$log2_seg),]
  out <- rep(NA_real_,nrow(map)); n <- rep(0L,nrow(map))
  for (ch in unique(stats::na.omit(map$t2t_chrom))) {
    bb <- b[b$chromosome==ch,]; if (!nrow(bb)) next
    ix <- which(map$t2t_chrom==ch)
    for (i in ix) {
      hit <- bb$start<=map$t2t_end[i] & bb$end>=map$t2t_start[i]
      if (any(hit)) {out[i] <- median(bb$log2_seg[hit]); n[i] <- sum(hit)}
    }
  }
  list(log2_seg=out,n_bins=n)
}

du145_y_par <- function(x,e,cfg) {
  a <- e$effects; m <- x$core; counts <- x$counts[,m$sample_key,drop=FALSE]
  spans <- t2t_gene_spans(cfg$t2t_genes)
  # Inclusive paired fit over all expressed nuclear genes (autosomes, X, Y incl. PAR), primary TMM exposures.
  k <- a$expressed & a$chrom %in% c(paste0('chr',1:22),'chrX','chrY')
  anno <- a[k,c('gene_id','gene_name','gene_type','chrom','arm','start','end','PAR','mean_difference','CTRL_mean_CPM','LOY_mean_CPM')]
  fit <- fit_voom(counts[k,],anno,e$norm$effective,m,paired=TRUE)
  inc <- fit$table; names(inc)[names(inc)=='FDR'] <- 'BH_all_expressed_nuclear'
  inc <- thesis_de_call(transform(inc,FDR=BH_all_expressed_nuclear),e$pair_diff[k,,drop=FALSE],cfg)
  tab(inc,'DE_inclusive_nuclear_paired_voom')
  saveRDS(fit,'outputs/models/inclusive_nuclear_paired_voom.rds')
  # Y and PAR family: every chrY or PAR-overlapping feature with any counts, tested where expressed.
  yp <- which((a$chrom=='chrY' | a$PAR) & rowSums(counts)>0)
  y <- a[yp,c('gene_id','gene_name','gene_type','chrom','arm','start','end','PAR','expressed','CTRL_mean_CPM','LOY_mean_CPM','mean_difference',grep('^pair_',names(a),value=TRUE))]
  y$region <- ifelse(y$PAR,ifelse(y$start<5e6,'PAR1','PAR2'),'MSY')
  st <- inc[match(y$gene_id,inc$gene_id),c('logFC','CI.L','CI.R','PValue','BH_all_expressed_nuclear')]
  y <- cbind(y,st)
  y$BH_Y_PAR_family <- NA_real_; tested <- is.finite(y$PValue)
  y$BH_Y_PAR_family[tested] <- p.adjust(y$PValue[tested],'BH')
  y$percent_of_sgCTRL <- ifelse(y$CTRL_mean_CPM>0,100*y$LOY_mean_CPM/y$CTRL_mean_CPM,NA_real_)
  cpm <- e$norm$CPM[y$gene_id,,drop=FALSE]
  for (b in x$pairs$block) {
    r <- x$pairs[x$pairs$block==b,]
    y[[paste0('percent_of_sgCTRL_',b)]] <- ifelse(cpm[,r$CTRL]>0,100*cpm[,r$LOY]/cpm[,r$CTRL],NA_real_)
  }
  mp <- map_to_t2t(cbind(a[yp,c('gene_id','gene_name','chrom','PAR')],stable_gene_id=x$annotation$stable_gene_id[yp]),spans,y_target=TRUE)
  y <- cbind(y,mp[,-1])
  y$t2t_mid_Mb <- (y$t2t_start+y$t2t_end)/2e6
  # lpWGS context along chromosome Y (saved segment values; no new segmentation).
  bins <- load_wgs_bins(cfg)
  ymap <- data.frame(t2t_chrom=ifelse(is.na(y$t2t_chrom),NA,'chrY'),t2t_start=y$t2t_start,t2t_end=y$t2t_end)
  for (s in c(cfg$wgs_ctrl,cfg$wgs_loy)) {
    z <- gene_bin_cn(ymap,bins,s); lab <- if (s==cfg$wgs_ctrl) 'sgCTRL' else 'LOY2X'
    y[[paste0('lpWGS_copy_ratio_',lab)]] <- 2^z$log2_seg; y[[paste0('lpWGS_bins_',lab)]] <- z$n_bins
  }
  y$lpWGS_LOY2X_over_sgCTRL <- y$lpWGS_copy_ratio_LOY2X/y$lpWGS_copy_ratio_sgCTRL
  seg_retained <- function(mb) ifelse(is.na(mb),NA_character_,ifelse(mb<2.45832,'PAR1',ifelse(mb<=4.8,'retained distal Yp (2.46-4.8 Mb)',ifelse(mb<62.122809,'depleted MSY','PAR2'))))
  y$lpWGS_interval <- seg_retained(y$t2t_mid_Mb)
  y <- y[order(y$t2t_mid_Mb,na.last=TRUE),]
  tab(y,'Y_PAR_gene_table')
  per <- do.call(rbind,lapply(which(y$expressed),function(i) data.frame(gene_id=y$gene_id[i],gene_name=y$gene_name[i],region=y$region[i],
    lpWGS_interval=y$lpWGS_interval[i],t2t_mid_Mb=y$t2t_mid_Mb[i],m[,c('sample_key','treatment','passage_block')],
    CPM=as.numeric(e$norm$CPM[y$gene_id[i],m$sample_key]),raw_count=as.integer(counts[y$gene_id[i],]))))
  tab(per,'Y_PAR_sample_abundance')
  # Region sums per library with the same core exposures.
  cls <- ifelse(a$PAR,ifelse(a$start<5e6,'PAR1','PAR2'),ifelse(a$chrom=='chrY','MSY',NA))
  msy_int <- setNames(y$lpWGS_interval,y$gene_id)[a$gene_id]
  cls[!is.na(cls) & cls=='MSY'] <- ifelse(!is.na(msy_int[!is.na(cls) & cls=='MSY']) & msy_int[!is.na(cls) & cls=='MSY']=='retained distal Yp (2.46-4.8 Mb)','MSY retained distal Yp','MSY depleted')
  reg <- do.call(rbind,lapply(stats::na.omit(unique(cls)),function(r) data.frame(region=r,m[,c('sample_key','treatment','passage_block')],
    sum_CPM=colSums(e$norm$CPM[which(cls==r),m$sample_key,drop=FALSE]))))
  tab(reg,'Y_PAR_region_abundance')
  # X-Y paralog pairs.
  pairs <- data.frame(X=c('KDM6A','KDM5C','DDX3X','EIF1AX','RPS4X','ZFX','USP9X','NLGN4X','TBL1X','PRKX','TMSB4X'),
                      Y=c('UTY','KDM5D','DDX3Y','EIF1AY','RPS4Y1','ZFY','USP9Y','NLGN4Y','TBL1Y','PRKY','TMSB4Y'))
  find <- function(sym,ch) {ix<-which(a$gene_name==sym & a$chrom==ch & !a$PAR); if(!length(ix)) NA_integer_ else ix[which.max(a$mean_CPM[ix])]}
  for (j in seq_len(nrow(pairs))) {
    for (s in c('X','Y')) {
      i <- find(pairs[[s]][j],if (s=='X') 'chrX' else 'chrY')
      pairs[j,paste0(s,'_CTRL_CPM')] <- if (is.na(i)) NA else a$CTRL_mean_CPM[i]
      pairs[j,paste0(s,'_LOY_CPM')] <- if (is.na(i)) NA else a$LOY_mean_CPM[i]
      pairs[j,paste0(s,'_mean_difference')] <- if (is.na(i)) NA else a$mean_difference[i]
      z <- if (is.na(i)) NULL else inc[inc$gene_id==a$gene_id[i],]
      pairs[j,paste0(s,'_logFC')] <- if (is.null(z) || !nrow(z)) NA else z$logFC
      pairs[j,paste0(s,'_BH_inclusive')] <- if (is.null(z) || !nrow(z)) NA else z$BH_all_expressed_nuclear
    }
  }
  pairs$Y_expressed_in_sgCTRL <- is.finite(pairs$Y_CTRL_CPM) & pairs$Y_CTRL_CPM>=1
  pairs$summed_XY_percent_of_sgCTRL <- ifelse(pairs$Y_expressed_in_sgCTRL,100*(pairs$X_LOY_CPM+pairs$Y_LOY_CPM)/(pairs$X_CTRL_CPM+pairs$Y_CTRL_CPM),NA)
  tab(pairs,'XY_paralog_effects')
  list(table=y,per_sample=per,regions=reg,pairs=pairs,inclusive=inc)
}
