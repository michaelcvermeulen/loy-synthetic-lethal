# Cell-line identity from the RNA itself: DepMap 26Q1 correlation and prostate lineage markers.
# The source records call this line LNCaP; the thesis reports it as DU145.
du145_identity <- function(x,e,cfg) {
  root <- cfg$depmap_root
  manifest <- data.table::fread(file.path(root,'docs/data_file_manifest.csv'))
  verify <- function(rel) {
    path <- file.path(root,rel); expected <- manifest$sha256[manifest$new_relative_path==rel]
    assert(length(expected)==1,paste('DepMap manifest entry missing:',rel))
    observed <- hash_file(path); assert(identical(observed,expected),paste('DepMap hash mismatch:',rel))
    data.frame(path=path,sha256=observed)
  }
  hashes <- rbind(verify(cfg$depmap_expression),verify(cfg$depmap_model))
  tab(hashes,'identity_depmap_input_hashes')
  model <- as.data.frame(data.table::fread(file.path(root,cfg$depmap_model)))
  expr <- data.table::fread(file.path(root,cfg$depmap_expression))
  expr <- expr[expr$IsDefaultEntryForModel %in% c('Yes',TRUE,'True')]
  ids <- expr$ModelID
  gene_cols <- setdiff(names(expr),c('V1','SequencingID','ModelConditionID','ModelID','IsDefaultEntryForMC','IsDefaultEntryForModel'))
  dm <- as.matrix(expr[,gene_cols,with=FALSE]); rownames(dm) <- ids
  colnames(dm) <- sub(' \\(.*\\)$','',gene_cols)
  dm <- dm[,!duplicated(colnames(dm)),drop=FALSE]
  # Our libraries on the DepMap scale: log2(TPM + 1) from counting-feature lengths.
  a <- x$annotation; pc <- a$gene_type=='protein_coding' & a$chrom!='chrM'
  rpk <- x$counts[pc,,drop=FALSE]/(a$feature_length_bp[pc]/1e3)
  tpm <- sweep(rpk,2,colSums(rpk),'/')*1e6
  sym <- a$gene_name[pc]
  keep <- !is.na(sym) & nzchar(sym) & !duplicated(sym)
  ours <- log2(tpm[keep,,drop=FALSE]+1); rownames(ours) <- sym[keep]
  common <- intersect(colnames(dm),rownames(ours))
  sds <- apply(dm[,common,drop=FALSE],2,sd,na.rm=TRUE)
  top <- names(sort(sds,decreasing=TRUE))[seq_len(min(cfg$depmap_top_variable,length(sds)))]
  profiles <- cbind(ours[top,,drop=FALSE],
    sgCTRL_mean=rowMeans(ours[top,x$meta$treatment=='CTRL',drop=FALSE]),
    LOY_mean=rowMeans(ours[top,x$meta$treatment=='LOY',drop=FALSE]))
  rho <- cor(t(dm[,top,drop=FALSE]),profiles,method='spearman',use='pairwise.complete.obs')
  info <- model[match(rownames(rho),model$ModelID),c('ModelID','StrippedCellLineName','OncotreeLineage','OncotreePrimaryDisease')]
  ranked <- lapply(colnames(rho),function(s) {
    d <- cbind(info,profile=s,spearman=rho[,s]); d <- d[order(-d$spearman),]; d$rank <- seq_len(nrow(d)); d
  })
  ranked <- do.call(rbind,ranked)
  tab(ranked,'identity_depmap_all_ranks')
  focus <- ranked[ranked$StrippedCellLineName %in% c('DU145','LNCAPCLONEFGC','PC3','VCAP','22RV1','C42B') | ranked$rank<=5,]
  tab(focus[order(focus$profile,focus$rank),],'identity_depmap_summary')
  markers <- c('AR','KLK3','KLK2','FOLH1','HOXB13','NKX3-1','TMPRSS2','KRT7','KRT8','KRT18','VIM')
  mk <- do.call(rbind,lapply(markers,function(g) {
    i <- which(a$gene_name==g & a$chrom!='chrM'); if (!length(i)) return(NULL); i <- i[which.max(rowSums(x$counts[i,,drop=FALSE]))]
    data.frame(gene_name=g,gene_id=a$gene_id[i],chrom=a$chrom[i],sample_key=colnames(x$counts),raw_count=as.integer(x$counts[i,]),TPM=as.numeric(tpm[match(a$gene_id[i],a$gene_id[pc]),]))
  }))
  tab(mk,'identity_lineage_markers')
  list(ranked=ranked,focus=focus,markers=mk,genes_used=length(top))
}
