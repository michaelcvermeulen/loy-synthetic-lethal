root <- normalizePath("SEPT_11_2026_THESIS", winslash="/")
base <- file.path(root,"research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-05/derived/loy_reanalysis_v1")
m <- readRDS(file.path(base,"matrices/prism_REP1M.rds"))
c <- read.csv(file.path(base,"cohort.csv"),check.names=FALSE)
s <- read.csv(file.path(root,"figures/chapter3/source_data/prism_rep1m_volcano/plotted_results.csv"),check.names=FALSE)
s <- s[tolower(s$labelled)=="true",]
c <- c[c$ModelID %in% rownames(m) & !is.na(c$y_call) & c$y_call!="",]
long <- do.call(rbind,lapply(seq_len(nrow(s)),function(i) {
    z <- c[,c("ModelID","PatientID","StrippedCellLineName","OncotreeLineage","y_call")]
    z$feature <- s$feature[i]
    z$compound_name <- s$compound_name[i]
    z$response <- m[z$ModelID,s$feature[i]]
    z <- z[is.finite(z$response),]
    stopifnot(sum(z$y_call=="Y_negative")==s$n_negative[i],sum(z$y_call=="Y_positive")==s$n_positive[i])
    stopifnot(abs(mean(z$response[z$y_call=="Y_negative"])-s$mean_negative[i])<1e-8)
    stopifnot(abs(mean(z$response[z$y_call=="Y_positive"])-s$mean_positive[i])<1e-8)
    z
}))
out <- file.path(root,"figures/chapter3/source_data/prism_rep1m_boxplots")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(long,file.path(out,"cell_line_responses.csv"),row.names=FALSE)
write.csv(s,file.path(out,"compound_statistics.csv"),row.names=FALSE)
cat(nrow(long),"measurements; counts and means match saved analysis\n")
