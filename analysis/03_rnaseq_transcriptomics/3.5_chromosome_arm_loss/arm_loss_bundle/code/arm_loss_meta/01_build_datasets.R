#!/usr/bin/env Rscript
# Arm-loss meta-analysis: load the three deposited count sets into one object.
#
#   GSE222379  ReDACT engineered 1q / 7p / 8q loss, 6 isogenic backgrounds (62 libs)
#   GSE336576  RPTEC 3p-loss clones, early (D18) + late (D90) passage (21 libs)
#   GSE150686  RPE1 monosomy 10;18 / 13 / 19p + inducible-p53 arms (33 libs)
#
# Replicate unit: technical replicates of the SAME clone are summed, so the unit
# of replication for a loss arm is the independent clone, not the library. Control
# libraries are replicates of one parental genotype and are kept as-is.
#
# Karyotypes come from the GEO characteristics, NOT the file names: GSE150686's
# "RM_10_ip53" is annotated "Monosomy Chr. 10 and 18", so naming it chr10-only
# would mis-specify which genes are cis.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "/global/home/hpc5173/projects/LOY/RNA/arm_loss_meta"
DATA <- file.path(ROOT, "data"); OUT <- file.path(ROOT, "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

strip_ver <- function(x) sub("\\..*$", "", x)

# GENCODE PAR_Y rows duplicate their chrX counterpart once the version suffix is
# stripped (44 of them in GSE336576). The masked-reference convention in this
# project keeps PAR on chrX, so the _PAR_Y copies are dropped.
# Files are NOT in a common gene order -- each deposited file is sorted by its own
# expression -- so every matrix is built by MATCHING on gene id, never by cbind.
read_counts <- function(f, header, id_col = 1, count_col = 2) {
  d <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = header)
  setnames(d, c(id_col, count_col), c("gene_id", "count"))
  d <- d[!grepl("_PAR_Y", gene_id)]
  d[, ensg := strip_ver(gene_id)]
  stopifnot(!anyDuplicated(d$ensg))
  setNames(as.numeric(d$count), d$ensg)
}
digest_chr <- function(x) {
  y <- sort(x); paste0(length(y), "_", substr(y[1], 1, 6), "_", substr(y[length(y)], 1, 6))
}
build_matrix <- function(vlist) {
  ref <- sort(names(vlist[[1]]))
  stopifnot(all(vapply(vlist, function(v) setequal(names(v), ref), logical(1))))
  m <- vapply(vlist, function(v) as.numeric(v[ref]), numeric(length(ref)))
  rownames(m) <- ref; storage.mode(m) <- "double"; m
}

# ================= GSE222379 =====================================================
f222 <- list.files(file.path(DATA, "GSE222379"), "\\.txt\\.gz$", full.names = TRUE)
l222 <- lapply(f222, read_counts, header = TRUE)
names(l222) <- sub("\\.txt\\.gz$", "", sub("^GSM[0-9]+_", "", basename(f222)))
# The A2780 libraries were quantified against a DIFFERENT annotation (48,162 genes
# vs 58,735; only 45,342 shared). Every contrast is within a cell line, so the
# blocks are kept separate rather than intersecting and throwing ~13k genes away.
blk222 <- vapply(l222, function(v) substr(digest_chr(names(v)), 1, 8), character(1))
cat("GSE222379 annotation blocks:\n"); print(table(blk222))
m222 <- lapply(split(names(l222), blk222), function(ids) build_matrix(l222[ids]))
names(m222) <- vapply(m222, function(m) if (ncol(m) == 10) "GSE222379_A2780" else "GSE222379_main",
                      character(1))

c222 <- data.table(sample = unlist(lapply(m222, colnames), use.names = FALSE))
c222[, base := sub("_rep_?[0-9]+$", "", sample)]                 # drops rep suffix only
c222[, cell_line := fcase(grepl("^MCF10A_hRas", base), "MCF10A_hRas",
                          default = sub("_.*$", "", base))]
c222[, tag := sub("^(MCF10A_hRas|[A-Za-z0-9]+)_", "", base)]      # e.g. 1q-loss_c1, WT, TK, TK_c1
c222[, condition := fcase(grepl("-loss", tag), "loss", grepl("^TK", tag), "TK",
                          default = "WT")]
c222[, arm_lost := fifelse(condition == "loss", sub("-loss.*$", "", tag), NA_character_)]
c222[, clone := fifelse(grepl("_c[0-9]+$", tag), sub(".*_(c[0-9]+)$", "\\1", tag), NA_character_)]
c222[, unit := fifelse(is.na(clone), sample, paste(base))]        # collapse reps within clone
c222[, dataset := "GSE222379"]

# ================= GSE336576 =====================================================
f336 <- list.files(file.path(DATA, "GSE336576"), "raw_counts", full.names = TRUE)
h336 <- fread(cmd = paste("zcat", shQuote(f336)), sep = "\t", header = TRUE, nrows = 0)
s336 <- names(h336)[-1]
l336 <- lapply(seq_along(s336), function(i) read_counts(f336, header = TRUE, count_col = i + 1))
names(l336) <- s336
m336 <- build_matrix(l336)

c336 <- data.table(sample = colnames(m336))
c336[, clone := sub("_.*$", "", sample)]
c336[, passage := fifelse(clone == "Ctrl", NA_character_, sub("^[^_]+_([^_]+)_.*$", "\\1", sample))]
c336[, condition := fifelse(clone == "Ctrl", "WT", "loss")]
c336[, arm_lost := fifelse(condition == "loss", "3p", NA_character_)]
c336[, cell_line := "RPTEC"]
c336[, unit := fifelse(condition == "loss", paste(clone, passage, sep = "_"), sample)]
c336[, dataset := "GSE336576"]

# ================= GSE150686 =====================================================
f150 <- list.files(file.path(DATA, "GSE150686"), "geneCounts\\.txt\\.gz$", full.names = TRUE)
l150 <- lapply(f150, read_counts, header = FALSE)
names(l150) <- sub("_geneCounts\\.txt\\.gz$", "", sub("^GSM[0-9]+_", "", basename(f150)))
m150 <- build_matrix(l150)

c150 <- data.table(sample = colnames(m150))
c150[, base := sub("_Rep[0-9]+$", "", sample)]
# karyotype + p53 background straight from the GEO characteristics
KARYO <- c(RPE1_p53KO = "diploid", RPE1_p53_KD = "diploid",
           RM_1018 = "chr10,chr18", RM_13 = "chr13", RM_19p = "19p",
           `RPE1_ip53_-D` = "diploid", `RPE1_ip53_+D` = "diploid",
           `RM_10_ip53_-D` = "chr10,chr18", `RM_10_ip53_+D` = "chr10,chr18",  # NOT chr10 alone
           `RM_13_ip53_-D` = "chr13", `RM_13_ip53_+D` = "chr13")
P53BG <- c(RPE1_p53KO = "p53KO", RM_1018 = "p53KO", RM_13 = "p53KO",
           RPE1_p53_KD = "p53KD", RM_19p = "p53KD",
           `RPE1_ip53_-D` = "ip53_off", `RM_10_ip53_-D` = "ip53_off", `RM_13_ip53_-D` = "ip53_off",
           `RPE1_ip53_+D` = "ip53_on", `RM_10_ip53_+D` = "ip53_on", `RM_13_ip53_+D` = "ip53_on")
stopifnot(all(c150$base %in% names(KARYO)), all(c150$base %in% names(P53BG)))
c150[, karyotype := KARYO[base]]
c150[, p53_bg := P53BG[base]]
c150[, condition := fifelse(karyotype == "diploid", "WT", "loss")]
c150[, arm_lost := fifelse(condition == "loss", karyotype, NA_character_)]
c150[, cell_line := "RPE1"]
c150[, clone := fifelse(condition == "loss", base, NA_character_)]
c150[, unit := sample]                       # 3 independent reps per line; no collapsing
c150[, dataset := "GSE150686"]


# ================= PRJNA930361 (ARPE-19 chrY elimination) ========================
# Our own realignment, not a GEO deposit: no count matrix was ever deposited for
# this study. STAR + featureCounts against the PAR-masked 1000G GRCh38 analysis
# set + GENCODE v49 (same annotation as the gene table used here), so the ids
# already carry a version suffix. Whole chrY is the lost segment.
# WEAK CONTRAST, label it as such: n=1 clone x3 replicates, and the comparator is
# the UNTREATED parental line -- no non-targeting or cutting-but-neutral guide,
# so Cas9 damage is not separated from Y dosage (same weakness as GSE336576).
f930 <- "/global/home/hpc5173/projects/LOY/RNA/updated/PRJNA930361/matrix/counts_matrix.tsv"
m930 <- fread(f930)
g930 <- strip_ver(m930[[1]]); m930 <- as.matrix(m930[, -1]); rownames(m930) <- g930
stopifnot(!anyDuplicated(rownames(m930)))

c930 <- data.table(sample = colnames(m930))
c930[, condition := fifelse(grepl("Ynull", sample), "loss", "WT")]
c930[, arm_lost  := fifelse(condition == "loss", "chrY", NA_character_)]
c930[, cell_line := "ARPE19"]
c930[, clone     := fifelse(condition == "loss", "A7123", NA_character_)]
c930[, unit      := sample]          # 3 independent replicates of one clone
c930[, dataset   := "PRJNA930361"]

# ================= assemble ======================================================
common <- c("dataset", "sample", "unit", "cell_line", "condition", "arm_lost", "clone")
fill <- function(d) { for (k in setdiff(common, names(d))) d[[k]] <- NA_character_; d }

ds <- list()
for (b in names(m222)) ds[[b]] <- list(counts = m222[[b]],
                                       coldata = c222[match(colnames(m222[[b]]), sample)])
ds[["GSE336576"]] <- list(counts = m336, coldata = c336)
ds[["GSE150686"]] <- list(counts = m150, coldata = c150)
ds[["PRJNA930361"]] <- list(counts = m930, coldata = c930)
for (b in names(ds)) ds[[b]]$coldata[, block := b]
cold <- rbindlist(lapply(ds, function(x) fill(copy(x$coldata))), fill = TRUE)

saveRDS(ds, file.path(OUT, "datasets.rds"))
fwrite(cold, file.path(OUT, "coldata_all.tsv"), sep = "\t")

for (k in names(ds)) {
  cat(sprintf("\n=== %s: %d genes x %d libs, %d replicate units\n", k,
              nrow(ds[[k]]$counts), ncol(ds[[k]]$counts), uniqueN(ds[[k]]$coldata$unit)))
  print(ds[[k]]$coldata[, .N, by = .(cell_line, condition, arm_lost)][order(cell_line, condition, arm_lost)])
}
cat("\nwrote", file.path(OUT, "datasets.rds"), "\n")
