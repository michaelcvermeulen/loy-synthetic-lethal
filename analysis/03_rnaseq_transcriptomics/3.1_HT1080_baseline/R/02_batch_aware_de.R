# Exploratory gene-level inference conditional on fixed populations/passages.
batch_aware_de <- function(x, e) {
  suppressPackageStartupMessages({library(limma); library(edgeR)})
  cfg <- jsonlite::fromJSON("config/batch_aware_de.json")
  m <- x$core
  pair_indices(m)
  keep <- e$effects$primary_eligible
  a <- e$effects[keep, , drop = FALSE]
  counts <- x$counts[a$gene_id, m$sample_key, drop = FALSE]
  effective <- setNames(e$norm$effective,colnames(e$logCPM))[m$sample_key]
  assert(nrow(m)==12 && nrow(counts)==25526, "Unexpected sample/gene universe")
  assert(identical(rownames(counts), a$gene_id) && all(is.finite(effective)), "Gene/normalisation mismatch")
  assert(all(counts>=0) && all(counts==round(counts)), "Expected original integer counts")
  primary <- camera_design(m)
  design <- primary$design
  contrasts <- primary$contrasts
  interaction <- setNames(rep(0,ncol(design)), colnames(design))
  interaction[c("LOY_AUG_2025","LOY_MAY_2026")] <- c(-1,1)
  contrasts <- cbind(contrasts, MAY_minus_AUG=interaction)
  m$treatment <- factor(m$treatment,levels=c("CTRL","LOY"))
  m$batch <- factor(m$batch,levels=c("AUG_2025","MAY_2026"))
  simple <- model.matrix(~batch+treatment,m)
  rownames(simple) <- m$sample_key
  simple_con <- matrix(0,ncol(simple),1,dimnames=list(colnames(simple),"combined"))
  simple_con["treatmentLOY",1] <- 1
  designs <- list(passage_period=design,batch_condition=simple)
  contrast_list <- list(passage_period=contrasts,batch_condition=simple_con)
  results <- fits <- summary <- list()
  annotation_cols <- c("gene_id","gene_name","gene_type","chrom","start","end","mean_CPM","mean_logCPM","mean_difference","AUG_2025","MAY_2026","descriptive_candidate","period_direction_agrees","LOO_all_same_direction")
  for (model in names(designs)) {
    D <- designs[[model]]; C <- contrast_list[[model]]
    assert(qr(D)$rank==ncol(D),paste("Rank-deficient model",model))
    v <- voom(counts,design=D,lib.size=effective,normalize.method="none",plot=FALSE,save.plot=TRUE)
    assert(all(is.finite(v$weights)) && all(v$weights>0), "Invalid voom weights")
    basefit <- lmFit(v,D)
    fit <- eBayes(contrasts.fit(basefit,C),trend=FALSE,robust=TRUE)
    tab(data.frame(sample_key=rownames(D),D,check.names=FALSE),paste0("batchDE_design_",model))
    tab(data.frame(coefficient=rownames(C),C,check.names=FALSE),paste0("batchDE_contrasts_",model))
    tab(data.frame(mean_log2_count=v$voom.xy$x,sqrt_residual_SD=v$voom.xy$y),paste0("batchDE_voom_points_",model))
    tab(data.frame(mean_log2_count=v$voom.line$x,sqrt_residual_SD=v$voom.line$y),paste0("batchDE_voom_trend_",model))
    tab(data.frame(sample_key=m$sample_key,effective_library_reads=effective,mean_weight=colMeans(v$weights),min_weight=apply(v$weights,2,min),max_weight=apply(v$weights,2,max)),paste0("batchDE_weights_",model))
    for (contrast in colnames(C)) {
      z <- topTable(fit,coef=contrast,number=Inf,sort.by="none",adjust.method="BH",confint=TRUE)
      z <- z[match(a$gene_id,rownames(z)), , drop=FALSE]
      names(z)[names(z)=="P.Value"] <- "PValue"
      names(z)[names(z)=="adj.P.Val"] <- "FDR"
      assert(all(is.finite(z$PValue)) && all(is.finite(z$FDR)), "Missing inferential results")
      assert(max(abs(z$FDR-p.adjust(z$PValue,"BH")))<1e-12,"BH mismatch")
      z <- cbind(a[,annotation_cols],model=model,contrast=contrast,z)
      z$FDR_below_0_05 <- z$FDR<cfg$alpha
      z$direction <- ifelse(z$logFC>0,"higher","lower")
      z$residual_df <- nrow(D)-ncol(D)
      z$statistical_scope <- "conditional_fixed_population_model_independent_residuals_not_generation_replication"
      key <- paste(model,contrast,sep="__")
      tab(z,paste0("batchDE_",key))
      results[[key]] <- z
      summary[[key]] <- data.frame(model=model,contrast=contrast,genes=nrow(z),residual_df=nrow(D)-ncol(D),FDR_below_0_05=sum(z$FDR_below_0_05),higher=sum(z$FDR_below_0_05 & z$logFC>0),lower=sum(z$FDR_below_0_05 & z$logFC<0),min_P=min(z$PValue),min_FDR=min(z$FDR))
    }
    fits[[model]] <- list(voom=v,fit=fit,basefit=basefit,design=D,contrasts=C)
  }
  tab(do.call(rbind,summary),"batchDE_summary")
  primary_result <- results$passage_period__combined
  sensitivity <- results$batch_condition__combined
  assert(identical(primary_result$gene_id,sensitivity$gene_id),"Model gene mismatch")
  comparison <- data.frame(gene_id=a$gene_id,gene_name=a$gene_name,
    descriptive_effect=a$mean_difference,passage_period_logFC=primary_result$logFC,
    passage_period_FDR=primary_result$FDR,batch_condition_logFC=sensitivity$logFC,
    batch_condition_FDR=sensitivity$FDR,descriptive_candidate=a$descriptive_candidate)
  tab(comparison,"batchDE_model_comparison")
  metrics <- data.frame(comparison=c("passage_period_vs_descriptive","passage_period_vs_batch_condition"),
    Pearson_r=c(cor(primary_result$logFC,a$mean_difference),cor(primary_result$logFC,sensitivity$logFC)),
    Spearman_r=c(cor(primary_result$logFC,a$mean_difference,method="spearman"),cor(primary_result$logFC,sensitivity$logFC,method="spearman")))
  tab(metrics,"batchDE_effect_concordance")
  tab(data.frame(gene_id=a$gene_id,gene_name=a$gene_name),"batchDE_gene_universe")
  # Confirm primary contrast is exactly the equal-weight combination of fitted effects.
  expected <- .5*(fits$passage_period$basefit$coefficients[,"LOY_AUG_2025"]+fits$passage_period$basefit$coefficients[,"LOY_MAY_2026"])
  assert(max(abs(expected-primary_result$logFC))<1e-10,"Combined contrast is not equally weighted")
  list(results=results,fits=fits,summary=do.call(rbind,summary),comparison=comparison,metrics=metrics,config=cfg)
}
