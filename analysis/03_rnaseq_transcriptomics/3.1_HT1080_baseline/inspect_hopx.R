library(SummarizedExperiment)
r <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE'
x <- readRDS(file.path(r,'outputs/models/raw_counts_SE.rds'))
print(class(x));print(names(rowData(x)));print(as.data.frame(rowData(x))[grepl('HOPX',as.character(rowData(x)$gene_name)),])
i <- which(rowData(x)$gene_name=='HOPX');print(assay(x)[i,,drop=FALSE])
d <- readRDS(file.path(r,'outputs/programme_analysis_v1/models/DE.rds'));print(names(d))
