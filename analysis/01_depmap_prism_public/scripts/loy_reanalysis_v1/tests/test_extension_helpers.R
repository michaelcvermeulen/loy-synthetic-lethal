# Unit checks for extension_helpers.R (scRNA environment; needs limma).
# Run: Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_extension_helpers.R

script_dir <- dirname(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
))))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
if (file.exists(file.path(script_dir, "extension_helpers.R"))) source(file.path(script_dir, "extension_helpers.R"))

n_pass <- 0L
n_fail <- 0L
check <- function(name, ok) {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
  } else {
    n_fail <<- n_fail + 1L
    cat("FAIL:", name, "\n")
  }
}
# A check whose expression stops counts as a failure, with the reason shown.
safe <- function(expr) tryCatch(expr, error = function(e) {
  cat("  error:", conditionMessage(e), "\n")
  FALSE
})
fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
tmp <- tempfile("ext_test_")
dir.create(tmp)

# --- guarded output ---------------------------------------------------------------
check("write_ext writes a prefixed table to the given folder", safe({
  write_ext(data.table(a = 1:2), "organoid_test", dir = tmp)
  identical(read_ext("organoid_test", dir = tmp)$a, 1:2)
}))
check("write_ext refuses a frozen table name", safe(fails(write_ext(data.table(a = 1), "cohort", dir = tmp))))
check("write_ext refuses any unprefixed name", safe(fails(write_ext(data.table(a = 1), "gene_annotation", dir = tmp))))
check("write_ext refuses the frozen folder even with a valid name", safe(fails(write_ext(data.table(a = 1), "organoid_test", dir = derived_dir))))

# --- second freeze ----------------------------------------------------------------
frozen_file <- file.path(tmp, "organoid_cohort.csv")
writeLines(c("a,b", "1,2"), frozen_file)
check("an extension freeze can be written and then holds", safe({
  write_extension_freeze("organoid", paths = c(cohort = frozen_file), sections = c("shared", "organoid"), dir = tmp)
  isTRUE(assert_extension_freeze("organoid", dir = tmp))
}))
check("a missing extension freeze stops the stage", safe(fails(assert_extension_freeze("beataml", dir = tmp))))
check("the freeze fails when a frozen table changes", safe({
  writeLines(c("a,b", "1,3"), frozen_file)
  changed <- fails(assert_extension_freeze("organoid", dir = tmp))
  writeLines(c("a,b", "1,2"), frozen_file)
  changed && isTRUE(assert_extension_freeze("organoid", dir = tmp))
}))
check("the freeze fails when its parameter section changes", safe({
  keep <- ext_params$organoid$selective_min
  ext_params$organoid$selective_min <- keep + 1
  changed <- fails(assert_extension_freeze("organoid", dir = tmp))
  ext_params$organoid$selective_min <- keep
  changed && isTRUE(assert_extension_freeze("organoid", dir = tmp))
}))
check("the freeze ignores a change in another parameter section", safe({
  keep <- ext_params$beataml$floor_sd
  ext_params$beataml$floor_sd <- keep + 1
  held <- isTRUE(assert_extension_freeze("organoid", dir = tmp))
  ext_params$beataml$floor_sd <- keep
  held
}))
check("a different freeze is not overwritten silently", safe({
  writeLines(c("a,b", "9,9"), frozen_file)
  refused <- fails(write_extension_freeze("organoid", paths = c(cohort = frozen_file), sections = c("shared", "organoid"), dir = tmp))
  write_extension_freeze("organoid", paths = c(cohort = frozen_file), sections = c("shared", "organoid"), dir = tmp, replace = TRUE)
  refused && isTRUE(assert_extension_freeze("organoid", dir = tmp))
}))

# --- downloads --------------------------------------------------------------------
download <- file.path(tmp, "download.txt")
writeLines("abc", download)
check("verify_download accepts the right checksum", safe(verify_download(download, digest::digest(file = download, algo = "md5"), "md5")$bytes == file.size(download)))
check("verify_download stops on a wrong checksum", safe(fails(verify_download(download, "00", "md5"))))
check("verify_download stops on a missing file", safe(fails(verify_download(file.path(tmp, "absent"), "00", "md5"))))

# --- RNA call on another dataset's scale ---------------------------------------------
calls <- read_result("y_status_calls")
thresholds <- jsonlite::fromJSON(file.path(derived_dir, "call_thresholds.json"))
core <- params$y_genes_core
with_rna <- calls[has_rna == TRUE]
expr <- as.matrix(with_rna[, paste0("rna_", core), with = FALSE])
dimnames(expr) <- list(with_rna$ModelID, core)
# The female maxima are taken from the data, as stage 03 does; the JSON copy is rounded.
cell_line_rule <- list(gene_max = apply(expr[with_rna$Sex == "Female", ], 2, max), depleted_mean_max = params$rna_call$depleted_mean_max, expressed_mean_min = params$rna_call$expressed_mean_min)
check("call_rna_y with the cell-line thresholds reproduces every stored RNA call", safe({
  out <- call_rna_y(expr, thresholds = cell_line_rule)
  identical(unname(out$call), with_rna$rna_call) && identical(unname(out$n_on), with_rna$rna_n_genes_on)
}))
expressed_fraction <- (2^params$rna_call$expressed_mean_min - 1) / (2^thresholds$expressed_male_median - 1)
derived <- safe(call_rna_y(expr, female_ids = with_rna[Sex == "Female", ModelID], male_ids = with_rna[Sex == "Male", ModelID], expressed_fraction = expressed_fraction))
check("the derived depleted threshold is the female maximum", safe(abs(derived$thresholds$depleted_mean_max - thresholds$female_mean_max) < 1e-6))
check("the derived expressed threshold lands on the cell-line threshold", safe(abs(derived$thresholds$expressed_mean_min - params$rna_call$expressed_mean_min) < 0.05))
check("derived thresholds give the stored calls for nearly every male line", safe({
  males <- with_rna$Sex == "Male"
  mean(derived$call[males] == with_rna$rna_call[males]) > 0.995
}))
check("call_rna_y stops when a reference female has male-level Y", safe({
  mislabelled <- with_rna[Sex == "Male" & rna_call == "expressed", ModelID][1]
  fails(call_rna_y(expr, female_ids = c(with_rna[Sex == "Female", ModelID], mislabelled), male_ids = with_rna[Sex == "Male", ModelID], expressed_fraction = expressed_fraction))
}))

# --- integration of RNA and DNA --------------------------------------------------------
check("integrate_y_calls reproduces every stored call and tier", safe({
  out <- integrate_y_calls(calls$rna_call, calls$dna_call, calls$Sex == "Male")
  same <- function(a, b) all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
  same(out$call, calls$y_call) && same(out$tier, calls$y_tier)
}))

# --- DNA call from copy-number segments --------------------------------------------------
segments <- rbind(
  data.table(sample = "lost", start = c(0, 500), end = c(500, 1000), ratio = c(0, 0.01)),
  data.table(sample = "kept", start = 0, end = 1000, ratio = 0.5),
  data.table(sample = "arm", start = c(0, 500), end = c(500, 1000), ratio = c(0, 0.5)),
  data.table(sample = "mosaic", start = 0, end = 1000, ratio = 0.1),
  data.table(sample = "sparse", start = 0, end = 100, ratio = 0),
  data.table(sample = "outside", start = c(-500, 900), end = c(100, 2000), ratio = c(0.5, 0.5))
)
dna <- safe(dna_segment_summary(segments, start = 0, end = 1000, loss_max = 0.03, retained_min = 0.2))
dna_call <- safe(setNames(call_dna_segments(dna, floor_min = 0.9, covered_min = 0.5), dna$sample))
check("segments at zero are a loss", safe(dna_call[["lost"]] == "loss"))
check("segments at one copy are retained", safe(dna_call[["kept"]] == "retained"))
check("half lost and half kept is partial", safe(dna_call[["arm"]] == "partial"))
check("an intermediate ratio throughout is mosaic", safe(dna_call[["mosaic"]] == "mosaic"))
check("a thinly covered sample gets no DNA call", safe(is.na(dna_call[["sparse"]])))
check("segments are clipped to the region", safe(abs(dna[sample == "outside", covered] - 0.2) < 1e-12 && abs(dna[sample == "arm", floor] - 0.5) < 1e-12))

# --- Y-loss fraction ----------------------------------------------------------------------
set.seed(11)
male_level <- matrix(rnorm(60 * 7, 6, 0.1), 60, 7, dimnames = list(paste0("m", 1:60), core))
typical <- apply(male_level, 2, median)
mixture <- function(alpha) typical + log2(alpha)
level <- rbind(male_level, mix75 = mixture(0.75), mix40 = mixture(0.4), mix10 = mixture(0.1), none = typical - 12, at_median = typical, above = typical + 1)
fraction <- safe(y_fraction(level, male_ids = rownames(male_level), floor_sd = 2))
value <- function(id, column = "f_raw") fraction[[column]][match(id, fraction$sample)]
check("y_fraction is linear in a mixture", safe(max(abs(value(c("mix75", "mix40", "mix10")) - c(0.25, 0.6, 0.9))) < 1e-8))
check("y_fraction is 0 at the male median and 1 without expression", safe(abs(value("at_median")) < 1e-8 && value("none") > 0.999))
check("y_fraction never leaves 0 to 1", safe(all(fraction$f_raw >= 0 & fraction$f_raw <= 1) && value("above") == 0))
check("values inside the male noise band are set to 0", safe({
  floor_value <- fraction$floor[1]
  floor_value > 0 && floor_value < 0.3 && all(fraction[f_raw < floor_value, f] == 0) && all(fraction[f_raw >= floor_value, f == f_raw])
}))
check("the floor is two robust SDs of the male bulk", safe({
  ratio <- log2(fraction$ratio[match(rownames(male_level), fraction$sample)])
  abs(fraction$floor[1] - (1 - 2^(-2 * robust_sd(ratio)))) < 1e-10
}))
check("one aberrant gene does not move the fraction", safe({
  odd <- rbind(level, odd = typical + c(-8, rep(0, 6)))
  abs(y_fraction(odd, male_ids = rownames(male_level), floor_sd = 2)[sample == "odd", f_raw]) < 1e-8
}))
check("females define a fraction of 1", safe({
  with_females <- rbind(level, f1 = typical - 9, f2 = typical - 9.2, f3 = typical - 8.8)
  out <- y_fraction(with_females, male_ids = rownames(male_level), female_ids = c("f1", "f2", "f3"), floor_sd = 2)
  abs(out[sample == "f1", f_raw] - 1) < 1e-8 && abs(out[sample == "mix40", f_raw] - 0.6) < 2e-3
}))

# --- karyotypes --------------------------------------------------------------------------
karyotypes <- c(
  "48,X,-Y,+4,t(6;9)(p22;q34),+8,+10,del(15)(q11.2q15)[20]" = 1,
  "45,X,-Y[19]/46,XY[1]" = 0.95,
  "45,X,-Y,t(8;21)(q22;q22)[4]/45,idem,del(9)(q12q22)[13]/46,XY[3]" = 0.85,
  "89<4n>,XXY,-Y,-7,-12,+13,+13,-17,-21[2]/46,XY[18]" = 0,
  "46,XY,t(10;11)(p13;q14)[11]/45,sl,-Y[7]/46,sdl1,t(1;8)(q12;p21)[2]" = 0.45,
  "43~44,X,-Y,t(3;11)(p21.3;p11.2),del(5)(q11.1;q31.1),add(12)(p11.2),-21[cp17]/44~45,idem,-t(3;11),+9[cp4]" = 1,
  "45,X,-Y,del(7)(p15)[4]/46,XY[16]" = 0.2,
  "45,X,-Y[3]/46,XY[17]" = 0.15,
  "45,X,-Y,t(8;21)(q22;q22)[20]" = 1,
  "45~49,X,-Y,add(5)(q11.2),+der(9;8)(q10;q10),+10,der(11)t(11;17)(p15;q11.2)del(11)(q21q22),add(12)(q12),del(13)(q12q22),-17,add(18)(p11.2),add(19)(q13.3),-20,+1~4mar[cp4]/46,XX[16]" = 0.2,
  "45,X,-Y,t(8;2;21)(q22;p21;q22)[20]" = 1,
  "45,X,-Y[20]" = 1,
  "43-45,X,-Y,del(5)(q22q35),del(7)(q22),+8,-9,add(9)(p13),+11,-13,-14,add(16)(q24),-17,-18,del(20)(q11.2q13.2),+2-3mar[cp21]" = 1,
  "45,X,-Y,der(2)t(2;5)(q37;p13),-5,-18,add(19)(p13.3),+21,+21[cp13]/45,sl,-der(2)t(2;5)(q37;p13),add(2)(q11.2),add(12)(q15)[7]" = 1,
  "40~42,X,-Y,dic(5;12)(q13;p11.2),del(7)(q22q36),add(9)(q31),add(13)(q22),ins(?;14)(?;q11.2),-15,add(15)(p11.2),del(16)(q11.1),-17,-18,add(19)(q13.4),-20,der(21)t(3;21)(q22.3;q13),+mar[cp12]/40~41,sl,add(2)(q31)[cp7]/46,XY[1] " = 0.95,
  "39~41,X,-Y,dic(5;12)(q13;p11.2),del(7)(q22q36),add(9)(q31),add(13)(q22),ins(?;14)(?;q11.2),-15,add(15)(p11.2),del(16)(q11.1),-17,-18.add(19)(q13.4),-20,der(21)t(3;21)(q22.3;q13)[cp16]/46,XY[5] " = 16 / 21,
  "45,X,-Y[4]/46,XY[16]" = 0.2,
  "45,X,-Y,del(11)(q12)[20]" = 1,
  "45-46,X,-Y,t(3;3)(q21;q26 .2),del(4)(q12),-5,del(6) (q21),del(9)(q22),del(11) (p11.2),ins(15;5)(q12;p15 .2p15.2),add(16)(q22),-17,-20,-21,-22,+mar1-5[c p18]/ 46,XY[2]" = 0.9,
  "46,XY,t(10;11)(p12;q23)[8]/45,sl,-Y[12]" = 0.6
)
parsed <- safe(parse_karyotype_y(names(karyotypes)))
check("every hand-checked karyotype with -Y is read", safe(all(parsed$parsed) && all(parsed$mentions_y_loss)))
check("the fraction of cells without a Y matches the hand count", safe(max(abs(parsed$fraction - unname(karyotypes))) < 1e-10))
check("cell totals are read from plain and composite counts", safe(identical(parsed$cells[c(1, 2, 6, 13)], c(20L, 20L, 21L, 21L))))
check("a clone that loses one of two Y chromosomes is still Y-positive", safe(parsed$fraction[4] == 0 && parsed$cells_without_y[4] == 0L))
check("a female clone in a male karyotype is flagged", safe(identical(which(parsed$female_clone), 10L)))
other <- safe(parse_karyotype_y(c("46,XY[20]", NA, "", "46,XX[20]", "46,XY[3]//46,XX[17]", "45,X,-Y", "92,XXYY[9]/46,XY[11]", "48,XYY,+8[21]", "normal", "46,XY[18]/45,X,-Y")))
check("a normal male karyotype has fraction 0", safe(other$parsed[1] && other$fraction[1] == 0 && other$cells[1] == 20L && !other$mentions_y_loss[1]))
check("missing and unreadable karyotypes are flagged, not guessed", safe(!any(other$parsed[c(2, 3, 9)]) && all(is.na(other$fraction[c(2, 3, 9)]))))
check("a wholly female karyotype is flagged", safe(other$female_clone[4] && other$fraction[4] == 0))
check("donor chimerism is flagged", safe(other$chimerism[5] && !any(other$chimerism[-5], na.rm = TRUE)))
check("a single clone without a count still gives its fraction", safe(other$parsed[6] && other$fraction[6] == 1 && is.na(other$cells[6]) && other$counts_missing[6]))
check("gained Y chromosomes are not a loss", safe(all(other$fraction[7:8] == 0)))
check("several clones with a count missing give no fraction", safe(other$parsed[10] && is.na(other$fraction[10]) && other$counts_missing[10] && other$mentions_y_loss[10]))
# Formats met in BeatAML2: a no-break space, a note after the last count, FISH results after the
# karyotype, a doubled clone, and a mistyped bracket.
messy <- safe(parse_karyotype_y(c(
  "46,XY[19]\u00a0", "47,XY,+8[6]/46,XY[13] (Specimen Source: Blood)", "46,XY,del(11)(q14q23)[13]/46,XY[7].ishdel(11)(5'MLL-,3'MLL+)[2].",
  "47,XY,+der(1;9)(q10;p10)[17]/94,idemx2[2]/46,XY[1]", "45,X,-Y[10]/90,idemx2[5]/46,XY[5]", "46,XY[20)"
)))
check("a no-break space after the count is ignored", safe(messy$parsed[1] && messy$cells[1] == 19L && messy$fraction[1] == 0))
check("a note after the last count is ignored", safe(messy$parsed[2] && messy$cells[2] == 19L && messy$fraction[2] == 0))
check("FISH results after the karyotype are ignored", safe(messy$parsed[3] && messy$cells[3] == 20L && messy$fraction[3] == 0))
check("a doubled clone inherits the sex chromosomes of its parent", safe(messy$parsed[4] && messy$fraction[4] == 0 && messy$parsed[5] && messy$fraction[5] == 0.75))
check("a mistyped count is treated as missing, not as unreadable", safe(messy$parsed[6] && messy$fraction[6] == 0 && messy$counts_missing[6]))

# --- context fits by line set ----------------------------------------------------------------
set.seed(5)
n_lines <- 150
ctx_frame <- data.table(
  ModelID = sprintf("L%03d", seq_len(n_lines)), OncotreeLineage = sample(c("A", "B", "C"), n_lines, replace = TRUE),
  exposed = rbinom(n_lines, 1, 0.45), ctx = rbinom(n_lines, 1, 0.3), LibraryKY = rbinom(n_lines, 1, 0.3)
)
ctx_frame[, interaction := exposed * ctx]
mat_raw <- matrix(rnorm(80 * n_lines), 80, n_lines, dimnames = list(sprintf("g%02d", 1:80), ctx_frame$ModelID))
mat_raw[1, ] <- mat_raw[1, ] - 1.5 * ctx_frame$interaction
rank_rows <- function(m) {
  out <- t(apply(m, 1, rank_normalise))
  dimnames(out) <- dimnames(m)
  out
}
mat_rank <- rank_rows(mat_raw)
reference <- suppressWarnings(fit_context_tests(ctx_frame, mat_rank, "LibraryKY"))
by_pattern <- safe(fit_context_patterns(ctx_frame, mat_rank, "LibraryKY", min_cell = 5))
check("fit_context_patterns equals fit_context_tests on a complete matrix", safe({
  max(abs(log10(by_pattern$interaction_p) - log10(reference$interaction_p)), abs(log10(by_pattern$inside_p) - log10(reference$inside_p))) < 1e-8 &&
    max(abs(by_pattern$interaction_estimate - reference$interaction_estimate), abs(by_pattern$inside_estimate - reference$inside_estimate)) < 1e-10
}))
check("fit_context_patterns reports the cells of the context", safe(all(by_pattern$n_neg_in == sum(ctx_frame$interaction)) && all(by_pattern$n == n_lines)))
holes <- mat_raw
set_a <- ctx_frame$ModelID[1:40]
set_b <- ctx_frame$ModelID[seq(3, 150, by = 3)]
holes[11:30, set_a] <- NA
holes[31:45, set_b] <- NA
# Gene 46 keeps only four Y-negative lines inside the context: below the gate.
inside_negative <- ctx_frame[interaction == 1, ModelID]
holes[46, inside_negative[-(1:4)]] <- NA
holes_rank <- rank_rows(holes)
groups <- safe(pattern_groups(holes_rank))
check("pattern_groups puts genes with the same lines together", safe(length(groups) == 4 && identical(sort(vapply(groups, function(g) length(g$features), integer(1))), c(1L, 15L, 20L, 44L))))
gapped <- safe(fit_context_patterns(ctx_frame, holes_rank, "LibraryKY", min_cell = 5))
lm_interaction <- function(gene) {
  data <- cbind(ctx_frame, y = holes_rank[gene, ])
  coef(summary(lm(y ~ interaction + exposed + ctx + OncotreeLineage + LibraryKY, data = data)))["interaction", ]
}
check("with missing lines each gene's estimate is the least-squares one on its own lines", safe({
  genes <- c("g05", "g12", "g40")
  max(abs(gapped[match(genes, feature), interaction_estimate] - vapply(genes, function(g) lm_interaction(g)[["Estimate"]], numeric(1)))) < 1e-10
}))
check("one variance prior is pooled over all genes of the context", safe({
  fitted <- gapped[!is.na(interaction_p)]
  stats <- t(vapply(fitted$feature, function(g) {
    data <- cbind(ctx_frame, y = holes_rank[g, ])
    fit <- lm(y ~ interaction + exposed + ctx + OncotreeLineage + LibraryKY, data = data)
    c(sigma(fit)^2, fit$df.residual, coef(summary(fit))["interaction", "Estimate"], coef(summary(fit))["interaction", "Std. Error"] / sigma(fit))
  }, numeric(4)))
  squeezed <- squeezeVar(stats[, 1], stats[, 2], robust = TRUE)
  df_total <- pmin(stats[, 2] + squeezed$df.prior, sum(stats[, 2]))
  p <- 2 * pt(-abs(stats[, 3] / stats[, 4] / sqrt(squeezed$var.post)), df_total)
  max(abs(log10(p) - log10(fitted$interaction_p))) < 1e-8
}))
check("a gene whose own lines leave a short cell is not tested", safe(is.na(gapped[feature == "g46", interaction_p]) && is.na(gapped[feature == "g46", inside_p]) && gapped[feature == "g46", n_neg_in] == 4))
check("the other genes keep their tests and their own counts", safe(sum(!is.na(gapped$interaction_p)) == 79 && gapped[feature == "g12", n] == n_lines - 40))
check("moderated_test equals limma's robust empirical Bayes without trend", safe({
  design <- build_design(ctx_frame, tested = "exposed", covariates = "LibraryKY")$design
  fit <- lmFit(mat_rank, design)
  eb <- eBayes(fit, trend = FALSE, robust = TRUE)
  out <- moderated_test(fit$coefficients[, "exposed"], fit$stdev.unscaled[, "exposed"], fit$sigma, fit$df.residual)
  max(abs(out$t - eb$t[, "exposed"])) < 1e-10 && max(abs(log10(out$p) - log10(eb$p.value[, "exposed"]))) < 1e-10
}))

# --- weighted scan ------------------------------------------------------------------------------
design <- build_design(ctx_frame, tested = "exposed", covariates = "LibraryKY")$design
plain <- limma_scan(mat_raw, design, "exposed")
equal <- safe(weighted_limma_scan(mat_raw, design, "exposed", weights = rep(1, n_lines)))
check("weighted_limma_scan with equal weights equals limma_scan", safe(max(abs(equal$estimate - plain$estimate), abs(equal$t - plain$t)) < 1e-10 && max(abs(log10(equal$p) - log10(plain$p))) < 1e-10))
set.seed(8)
sample_weights <- runif(n_lines, 0.3, 2)
weighted <- safe(weighted_limma_scan(mat_raw, design, "exposed", weights = sample_weights))
check("weighted_limma_scan estimates equal weighted least squares", safe({
  genes <- c(1, 17, 63)
  by_lm <- vapply(genes, function(i) coef(lm(mat_raw[i, ] ~ 0 + design, weights = sample_weights))[["designexposed"]], numeric(1))
  max(abs(weighted$estimate[genes] - by_lm)) < 1e-10
}))
check("the design unit and the standardised effect are returned", safe({
  unit <- sqrt(solve(crossprod(design * sqrt(sample_weights)))["exposed", "exposed"])
  max(abs(weighted$unit - unit)) < 1e-10 && max(abs(weighted$smd - weighted$t * weighted$unit)) < 1e-10
}))
check("weighted_limma_scan refuses missing values", safe(fails(weighted_limma_scan(holes, design, "exposed", weights = sample_weights))))

# --- combination across sources -------------------------------------------------------------------
check("signed_z carries the sign of the estimate", safe(all.equal(signed_z(c(0.05, 0.05, 1), c(2, -2, 3)), c(qnorm(0.975), -qnorm(0.975), 0))))
check("signed_z keeps a finite value at p = 0", safe(is.finite(signed_z(0, -1)) && signed_z(0, -1) < -30))
z_two <- rbind(c(2, 2), c(2, -2), c(3, NA), c(-1.5, -2.5))
units <- c(0.1, 0.25)
pooled <- safe(stouffer(z_two, 1 / units))
check("stouffer of two equal sources is z times root two", safe(abs(stouffer(matrix(c(2, 2), 1), c(1, 1))$z - 2 * sqrt(2)) < 1e-12))
check("stouffer equals fixed-effect pooling of standardised effects", safe({
  d <- sweep(z_two, 2, units, "*")
  weight <- 1 / units^2
  by_hand <- vapply(c(1, 2, 4), function(i) sum(weight * d[i, ]) / sum(weight) * sqrt(sum(weight)), numeric(1))
  max(abs(pooled$z[c(1, 2, 4)] - by_hand)) < 1e-12
}))
check("stouffer reverses with the sign and leaves a single source as it is", safe(abs(stouffer(-z_two, 1 / units)$z[1] + pooled$z[1]) < 1e-12 && pooled$z[3] == 3 && pooled$n_sources[3] == 1L))
check("stouffer p-values are two-sided", safe(abs(pooled$p[1] - 2 * pnorm(-abs(pooled$z[1]))) < 1e-15))
difference <- safe(source_difference(c(0.5, 0.1), c(0.1, 0.2)))
check("source_difference is the difference with its pooled standard error", safe(abs(difference$difference - 0.4) < 1e-12 && abs(difference$se - sqrt(0.05)) < 1e-12 && abs(difference$p - 2 * pnorm(-0.4 / sqrt(0.05))) < 1e-12))
classes <- safe(evidence_class(rbind(c(4, 2.5), c(5, 0.3), c(5, -2.5), c(-4, -2.2), c(4, NA)), c(10, 4), alpha = 0.05))
check("evidence_class separates support, one source and disagreement", safe(identical(classes, c("supported by both", "carried by one source", "sources disagree", "supported by both", "carried by one source"))))

# --- least squares and permutation for the drug screen ------------------------------------------------
set.seed(21)
n_spec <- 140
covariate <- rnorm(n_spec)
exposure <- pmax(0, 0.6 * covariate + rnorm(n_spec)) * rbinom(n_spec, 1, 0.15)
nuisance <- cbind(1, covariate)
response <- matrix(rnorm(n_spec * 6), n_spec, 6, dimnames = list(NULL, paste0("d", 1:6))) + covariate
response[, 1] <- response[, 1] - 1.2 * exposure
response[c(3, 50, 51), 2] <- NA
ols <- safe(scan_ols(response, exposure, nuisance))
check("scan_ols equals lm per column, with missing responses dropped", safe({
  by_lm <- vapply(1:6, function(j) coef(summary(lm(response[, j] ~ exposure + covariate)))["exposure", c("Estimate", "Pr(>|t|)")], numeric(2))
  max(abs(ols$estimate - by_lm[1, ]), abs(ols$p - by_lm[2, ])) < 1e-10 && ols$n[2] == n_spec - 3
}))
permutations <- safe(permutation_matrix(n_spec, 4000, seed = 4))
check("permutation_matrix holds permutations and is reproducible", safe(all(apply(permutations[, 1:50], 2, function(v) identical(sort(v), seq_len(n_spec)))) && identical(permutations, permutation_matrix(n_spec, 4000, seed = 4))))
keep <- setdiff(seq_len(n_spec), c(3, 50, 51))
restricted <- safe(restrict_permutations(permutations, keep))
check("restrict_permutations gives permutations of the kept samples", safe(nrow(restricted) == length(keep) && all(apply(restricted[, 1:50], 2, function(v) identical(sort(v), seq_along(keep))))))
check("a restricted permutation keeps the order of the full one", safe(identical(order(permutations[keep, 7]), order(restricted[, 7]))))
real <- safe(freedman_lane(response[, 1], exposure, nuisance, permutations))
check("freedman_lane finds a real effect and reports the least-squares statistic", safe(real$p < 0.01 && abs(real$t - ols$t[1]) < 1e-10 && length(real$t_perm) == 4000))
check("freedman_lane agrees with the classical test when errors are normal", safe({
  out <- freedman_lane(response[, 3], exposure, nuisance, permutations)
  abs(out$p - ols$p[3]) < 0.05
}))
check("freedman_lane is uniform under a null with a correlated covariate", safe({
  set.seed(33)
  small <- permutation_matrix(n_spec, 999, seed = 9)
  p_null <- vapply(seq_len(300), function(i) freedman_lane(covariate + rt(n_spec, 3), exposure, nuisance, small)$p, numeric(1))
  ks.test(p_null, "punif")$p.value > 0.01 && abs(mean(p_null < 0.05) - 0.05) < 0.035
}))

# --- follow-up checks -----------------------------------------------------------------------------------
set.seed(40)
n_fu <- 90
strata <- rep(c("a", "b"), each = n_fu / 2)
group <- rep(c(0, 1), n_fu / 2)
fu_nuisance <- cbind(1, as.numeric(strata == "b"))
stable_y <- rnorm(n_fu) - 1.5 * group
stable <- safe(followup_checks(stable_y, group, fu_nuisance, strata, n_perm = 2000, seed = 2))
check("a stable effect passes the follow-up checks", safe(stable$loo_sign_stable && stable$loo_max_p < 0.01 && stable$permutation_p < 0.01 && abs(stable$estimate + 1.5) < 0.5))
outlier_y <- rnorm(n_fu, sd = 0.3)
outlier_y[2] <- -14
driven <- safe(followup_checks(outlier_y, group, fu_nuisance, strata, n_perm = 2000, seed = 2))
check("an effect carried by one sample fails leave-one-out", safe(driven$loo_max_p > 0.05 && driven$loo_min_abs < abs(driven$estimate) / 2))
check("follow-up estimates equal least squares, also with weights", safe({
  w <- runif(n_fu, 0.5, 2)
  out <- followup_checks(stable_y, group, fu_nuisance, strata, weights = w, n_perm = 200, seed = 2)
  abs(out$estimate - coef(lm(stable_y ~ group + strata, weights = w))[["group"]]) < 1e-10 && abs(stable$estimate - coef(lm(stable_y ~ group + strata))[["group"]]) < 1e-10
}))
check("the follow-up permutation is reproducible under its seed", safe(identical(stable$permutation_p, followup_checks(stable_y, group, fu_nuisance, strata, n_perm = 2000, seed = 2)$permutation_p)))

# --- report sections ---------------------------------------------------------------------------------------
html <- file.path(tmp, "report.html")
writeLines(c("<html><body><p>front</p>", "<h1 id='a'>Summary</h1><p>One <b>two</b>.</p>", "<h1>Limits</h1><ul><li>three</li></ul>", "</body></html>"), html)
check("report_sections splits the report at its top-level headings", safe({
  out <- report_sections(html)
  identical(names(out), c("Summary", "Limits")) && out[["Summary"]] == "One two ." && out[["Limits"]] == "three"
}))
check("robust_sd is the scaled median absolute deviation", safe(abs(robust_sd(c(1, 2, 3, 4, 100)) - 1.4826 * 1) < 1e-12))

# --- report checks used by stage 23 ---------------------------------------------------------
blocks_html <- file.path(tmp, "blocks.html")
writeLines(c("<html><head><style>p { color: red; }</style></head><body>", "<h1>Summary</h1><p>One   <b>two</b>.</p>",
  "<ul><li><p>nested item</p></li><li>plain item</li></ul>", "<table><caption>A caption</caption><tr><td>a</td><td>1</td></tr></table>",
  "<p><img src='data:image/png;base64,AAAA'/></p>", "</body></html>"), blocks_html)
check("report_blocks returns each paragraph, list item, table row and heading once, in order", safe({
  identical(report_blocks(blocks_html), c("Summary", "One two.", "nested item", "plain item", "A caption", "a 1"))
}))
rmd <- c("# First", "text 1", "```{r}", "# a comment, not a heading", "x <- 2", "```", "## Sub", "more", "# Second", "other")
check("rmd_section returns one top-level section and is not fooled by comments in code", safe({
  identical(rmd_section(rmd, "First"), rmd[2:8]) && identical(rmd_section(rmd, "Second"), "other")
}))
check("rmd_section stops when the heading is missing", safe(inherits(try(rmd_section(rmd, "Third"), silent = TRUE), "try-error")))
prose <- c("Text with 3 genes and `r 5 + 2` inline, TP53 and Figure 11a at 0.05.", "```{r}", "x <- 42", "```", "On 2026-10-04, 20,000 shuffles; BeatAML2 at 80%.")
check("prose_numbers lists the numbers typed into the prose, not those in code or names", safe({
  identical(prose_numbers(prose), c("3", "0.05", "2026", "10", "04", "20,000", "80"))
}))
check("prose_numbers does not read a translocation name as numbers", safe(identical(prose_numbers("Without the t(8;21) cases and inv(16), 3 remain."), "3")))

# --- a table that only repeats an earlier result keeps its bytes ------------------------------
stable_dir <- file.path(tmp, "stable")
dir.create(stable_dir)
first <- data.table(symbol = c("A", "B"), estimate = c(-0.5, 1.25), p = c(3.7e-64, 0.02))
bytes_of <- function(name) unname(tools::md5sum(file.path(stable_dir, paste0(name, c(".csv", ".rds")))))
check("write_result_stable writes a table that is not stored yet", safe({
  written <- write_result_stable(first, "extension_refit", dir = stable_dir)
  isTRUE(written) && all(file.exists(file.path(stable_dir, c("extension_refit.csv", "extension_refit.rds"))))
}))
stored_bytes <- bytes_of("extension_refit")
check("write_result_stable keeps the stored bytes when values differ only in their last digits", safe({
  noisy <- copy(first)[, estimate := estimate * (1 + 1e-13)][, p := p * (1 - 1e-12)]
  written <- write_result_stable(noisy, "extension_refit", dir = stable_dir)
  identical(written, FALSE) && identical(bytes_of("extension_refit"), stored_bytes)
}))
check("write_result_stable rewrites when a small p-value changes by more than the tolerance", safe({
  moved <- copy(first)[1, p := 3.8e-64]
  written <- write_result_stable(moved, "extension_refit", dir = stable_dir)
  isTRUE(written) && read_result("extension_refit", dir = stable_dir)$p[1] == 3.8e-64
}))
check("write_result_stable rewrites when a label or the number of rows changes", safe({
  renamed <- copy(first)[1, p := 3.8e-64][2, symbol := "C"]
  isTRUE(write_result_stable(renamed, "extension_refit", dir = stable_dir)) && isTRUE(write_result_stable(renamed[1], "extension_refit", dir = stable_dir)) &&
    nrow(read_result("extension_refit", dir = stable_dir)) == 1
}))
check("write_result_stable treats a value that became missing as a change", safe({
  one <- read_result("extension_refit", dir = stable_dir)
  isTRUE(write_result_stable(copy(one)[, estimate := NA_real_], "extension_refit", dir = stable_dir))
}))

# --- NextGen (stages 31 to 33) ------------------------------------------------------------------
check("write_ext accepts the nextgen prefix", safe({
  write_ext(data.table(a = 1), "nextgen_test", dir = tmp)
  identical(read_ext("nextgen_test", dir = tmp)$a, 1)
}))

# Donors: HCMI donor ID, name root, genotype.
check("hcmi_donor reads the donor part of an HCMI name and nothing from other names", safe(
  identical(hcmi_donor(c("HCM-SANG-0299-C15-B", "HCM-BROD-0115-C16", "PANFR0440", NA, "WTSI-COLIVM_080")), c("HCM-SANG-0299", "HCM-BROD-0115", NA, NA, NA))
))
mutated <- rbind(a = c(1, 1, 1, 0, 0, 0, 1, 1), b = c(1, 1, 1, 0, 0, 0, 1, 0), c = c(0, 0, 0, 1, 1, 0, 1, 0), d = c(0, 0, 0, 0, 0, 0, 1, 0), e = c(0, 0, 0, 0, 0, 1, 1, 2))
similarity <- safe(genotype_similarity(mutated, max_fraction = 0.6))
check("genotype_similarity is the Jaccard index on calls that are rare across models", safe(
  abs(similarity["a", "b"] - 0.75) < 1e-12 && similarity["a", "c"] == 0 && abs(similarity["a", "e"] - 0.2) < 1e-12 && isSymmetric(similarity)
))
check("a model without any rare call is similar to nobody", safe(all(similarity["d", c("a", "b", "c", "e")] == 0) && !anyNA(similarity)))
donor_ids <- paste0("m", 1:6)
donor_keys <- list(hcmi = c("HCM-X-1", "HCM-X-1", NA, NA, NA, NA), root = c("A", "B", "P1", "P1", "Q", "R"))
donor_similarity <- matrix(0, 6, 6, dimnames = list(donor_ids, donor_ids))
donor_similarity["m5", "m6"] <- donor_similarity["m6", "m5"] <- 0.5
donor_similarity["m3", "m6"] <- donor_similarity["m6", "m3"] <- 0.1
grouped <- safe(donor_groups(donor_ids, donor_keys, donor_similarity, cutoff = 0.4))
check("donor_groups joins models by a shared ID, a shared name root or a genotype above the cut", safe(
  identical(unname(grouped), c("m1", "m1", "m3", "m3", "m5", "m5"))
))
check("two missing IDs are never the same donor", safe(length(unique(donor_groups(donor_ids, list(hcmi = rep(NA_character_, 6)), NULL, cutoff = 0.4))) == 6))
check("donor_groups is transitive and does not depend on row order", safe({
  chained <- donor_similarity
  chained["m4", "m5"] <- chained["m5", "m4"] <- 0.45
  forward <- donor_groups(donor_ids, donor_keys, chained, cutoff = 0.4)
  shuffled <- c(4, 6, 1, 5, 2, 3)
  backward <- donor_groups(donor_ids[shuffled], lapply(donor_keys, `[`, shuffled), chained, cutoff = 0.4)
  identical(unname(forward), c("m1", "m1", "m3", "m3", "m3", "m3")) && identical(unname(backward[order(shuffled)]), unname(forward))
}))
check("a similarity matrix that covers only some models leaves the others to the keys", safe(
  identical(unname(donor_groups(donor_ids, donor_keys, donor_similarity[5:6, 5:6], cutoff = 0.4)), c("m1", "m1", "m3", "m3", "m5", "m5"))
))
check("best_of_group keeps the better screen of each donor whatever the row order", safe({
  group <- c("m1", "m1", "m3", "m3", "m5")
  quality <- c(-6.1, -8.6, -4.7, -7.3, -5)
  keep <- best_of_group(donor_ids[1:5], group, quality)
  flipped <- 5:1
  identical(donor_ids[1:5][keep], c("m2", "m4", "m5")) && identical(sort(donor_ids[1:5][flipped][best_of_group(donor_ids[1:5][flipped], group[flipped], quality[flipped])]), c("m2", "m4", "m5"))
}))

# The RNA rule where recorded females show no Y signal at all.
set.seed(31)
trace_expr <- rbind(
  matrix(0, 12, 7), matrix(rnorm(70, 4, 0.3), 10, 7),
  trace = c(0.03, 0, 0.05, 0, 0.02, 0, 0), low = rep(1, 7), none = rep(0, 7)
)
trace_expr[3, 2] <- 0.05
dimnames(trace_expr) <- list(c(paste0("f", 1:12), paste0("m", 1:10), "trace", "low", "none"), paste0("Y", 1:7))
trace_females <- paste0("f", 1:12)
trace_males <- c(paste0("m", 1:10), "trace", "low", "none")
plain_rule <- safe(call_rna_y(trace_expr, female_ids = trace_females, male_ids = trace_males, expressed_fraction = 0.174))
ceiling_rule <- safe(call_rna_y(trace_expr, female_ids = trace_females, male_ids = trace_males, expressed_fraction = 0.174, ceiling_fraction = 0.02))
check("without a ceiling, a trace of Y signal in three genes is not called depleted", safe(plain_rule$call[["trace"]] == "intermediate" && plain_rule$call[["none"]] == "depleted"))
check("with the ceiling, a trace of Y signal is depleted and a tenth of the male level is not", safe(
  ceiling_rule$call[["trace"]] == "depleted" && ceiling_rule$call[["low"]] == "intermediate" && ceiling_rule$call[["none"]] == "depleted" && all(ceiling_rule$call[paste0("m", 1:10)] == "expressed")
))
check("the ceiling is the stated share of the typical male level", safe({
  expected <- log2(0.02 * (2^ceiling_rule$thresholds$male_reference - 1) + 1)
  abs(ceiling_rule$thresholds$depleted_mean_max - expected) < 1e-12 && ceiling_rule$thresholds$depleted_mean_max > plain_rule$thresholds$depleted_mean_max
}))
check("the ceiling never lowers a threshold the females set", safe({
  noisy <- trace_expr
  noisy["f1", ] <- 0.9
  high <- call_rna_y(noisy, female_ids = trace_females, male_ids = trace_males, expressed_fraction = 0.174, ceiling_fraction = 0.02)
  abs(high$thresholds$depleted_mean_max - 0.9) < 1e-12 && all(high$thresholds$gene_max >= 0.9)
}))

# chrY copy number from a gene-level table.
gene_cn <- rbind(lost = rep(0, 7), nearly = c(rep(0, 6), 0.5), kept = rep(0.5, 7), arm = c(0, 0, 0, 0.5, 0.5, 0.5, 0.5), mosaic = rep(0.1, 7), sparse = c(0, 0, 0, 0, 0, NA, NA))
gene_call <- safe(dna_gene_call(gene_cn, loss_max = 0.03, retained_min = 0.20, min_genes = 6))
check("dna_gene_call needs six of seven genes for loss or retained", safe(
  identical(gene_call[match(c("lost", "nearly", "kept"), sample), call], c("loss", "loss", "retained"))
))
check("dna_gene_call marks a split or an in-between profile as unclear, and too few genes as missing", safe(
  identical(gene_call[match(c("arm", "mosaic"), sample), call], c("partial", "mosaic")) && is.na(gene_call[sample == "sparse", call]) && abs(gene_call[sample == "nearly", ratio] - 0.5 / 7) < 1e-12
))

# Ranks within a source, the gate for pooled cancer types, and the stratified test.
set.seed(32)
pooled_mat <- matrix(rnorm(5 * 10), 5, 10, dimnames = list(paste0("g", 1:5), paste0("s", 1:10)))
pooled_source <- rep(c("cell", "nextgen"), c(6, 4))
check("rank_within_source equals rank_normalise applied source by source", safe({
  out <- rank_within_source(pooled_mat, pooled_source)
  by_hand <- cbind(t(apply(pooled_mat[, 1:6], 1, rank_normalise)), t(apply(pooled_mat[, 7:10], 1, rank_normalise)))
  max(abs(out - by_hand)) < 1e-12 && identical(dimnames(out), dimnames(pooled_mat))
}))
gate_counts <- data.table(
  type = rep(c("Pancreas", "Bowel", "Prostate"), c(3, 3, 2)), source = c("cell", "nextgen", "sanger", "cell", "nextgen", "sanger", "cell", "nextgen"),
  n_negative = c(21, 7, 0, 18, 4, 20, 1, 2), n_positive = c(3, 3, 1, 11, 4, 25, 4, 6)
)
check("pooled_gate keeps a cancer type only when no source can test it and the sources together can", safe(identical(pooled_gate(gate_counts, min_group = 5), "Pancreas")))

set.seed(33)
strat_mat <- matrix(rnorm(5 * 7), 5, 7, dimnames = list(paste0("g", 1:5), paste0("s", 1:7)))
strat_exposed <- c(1, 1, 0, 0, 1, 0, 0)
strat_strata <- c("a", "a", "a", "a", "b", "b", "b")
strat <- safe(stratified_permutation_test(strat_mat, strat_exposed, strat_strata))
arrangements <- as.matrix(expand.grid(first = seq_len(choose(4, 2)), second = seq_len(choose(3, 1))))
first_sets <- combn(4, 2)
brute <- sapply(seq_len(nrow(arrangements)), function(k) {
  x <- numeric(7)
  x[first_sets[, arrangements[k, 1]]] <- 1
  x[4 + arrangements[k, 2]] <- 1
  apply(strat_mat, 1, function(y) coef(lm(y ~ factor(strat_strata) + x))[["x"]])
})
brute_observed <- apply(strat_mat, 1, function(y) coef(lm(y ~ factor(strat_strata) + strat_exposed))[["strat_exposed"]])
brute_p <- pmin(1, 2 * pmin(rowMeans(brute >= brute_observed - 1e-12), rowMeans(brute <= brute_observed + 1e-12)))
check("the stratified estimate equals lm with the stratum as a factor", safe(max(abs(strat$table$estimate - brute_observed)) < 1e-10))
check("the exact p-value equals a brute-force enumeration of every arrangement", safe(
  strat$exact && strat$arrangements == 18 && max(abs(strat$table$p - brute_p)) < 1e-12
))
check("no exact p-value falls below two over the number of arrangements", safe(abs(strat$floor - 2 / 18) < 1e-12 && all(strat$table$p >= strat$floor - 1e-12)))
check("the family-wise p-value is the share of arrangements whose largest statistic is at least as large", safe({
  spread <- sqrt(rowMeans(brute^2))
  largest <- apply(abs(brute / spread), 2, max)
  expected <- vapply(abs(brute_observed / spread), function(z) mean(largest >= z - 1e-12), numeric(1))
  max(abs(strat$table$p_family - expected)) < 1e-12
}))
check("a stratum with one Y group changes nothing", safe({
  wider <- cbind(strat_mat, matrix(rnorm(5 * 3), 5, 3, dimnames = list(NULL, paste0("t", 1:3))))
  with_empty <- stratified_permutation_test(wider, c(strat_exposed, 1, 1, 1), c(strat_strata, "c", "c", "c"))
  with_empty$arrangements == 18 && max(abs(with_empty$table$estimate - strat$table$estimate)) < 1e-10 && max(abs(with_empty$table$p - strat$table$p)) < 1e-12
}))
check("too many arrangements are sampled, and the sampled p-value is close to the exact one", safe({
  bigger_x <- c(1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0)
  bigger_strata <- rep(c("a", "b"), each = 6)
  bigger <- matrix(rnorm(3 * 12), 3, 12, dimnames = list(paste0("g", 1:3), NULL))
  bigger[1, bigger_x == 1] <- bigger[1, bigger_x == 1] + 1.2
  exact <- stratified_permutation_test(bigger, bigger_x, bigger_strata)
  sampled <- stratified_permutation_test(bigger, bigger_x, bigger_strata, max_arrangements = 10, n_random = 20000, seed = 5)
  exact$exact && exact$arrangements == 400 && !sampled$exact && max(abs(sampled$table$p - exact$table$p)) < 0.03 && all(sampled$table$p >= 2 / 20001 - 1e-12)
}))
check("one feature alone gets the estimate and p-value it has among others, enumerated or sampled", safe({
  alone <- stratified_permutation_test(strat_mat[2, , drop = FALSE], strat_exposed, strat_strata)
  one_stratum <- stratified_permutation_test(strat_mat[2, 1:4, drop = FALSE], strat_exposed[1:4], strat_strata[1:4])
  sampled_alone <- stratified_permutation_test(strat_mat[2, , drop = FALSE], strat_exposed, strat_strata, max_arrangements = 5, n_random = 4000, seed = 3)
  abs(alone$table$estimate - strat$table$estimate[2]) < 1e-12 && abs(alone$table$p - strat$table$p[2]) < 1e-12 && alone$table$feature == "g2" &&
    one_stratum$arrangements == 6 && abs(one_stratum$table$estimate - (mean(strat_mat[2, 1:2]) - mean(strat_mat[2, 3:4]))) < 1e-12 &&
    !sampled_alone$exact && abs(sampled_alone$table$p - alone$table$p) < 0.05
}))
check("the stratified test stops on a missing value", safe(fails(stratified_permutation_test(replace(strat_mat, 3, NA), strat_exposed, strat_strata))))

# Calibration bar, replication outcome and the label of a combined hit.
check("family_wise_bar is the smallest p-value that the stated share of null runs reaches", safe(
  family_wise_bar((1:100) / 100, 0.05) == 0.05 && family_wise_bar(c(0.5, 1e-4, 0.2, 1e-6, 0.03, rep(0.9, 295)), 0.05) == 0.9 && family_wise_bar(c(1e-6, 1e-5, rep(0.9, 38)), 0.05) == 1e-5
))
check("replication_outcome gives each of its three answers", safe(identical(
  replication_outcome(estimate = c(-0.5, -0.05, -0.3, 0.5), adjusted_p = c(0.01, 0.6, 0.2, 0.01), reference = c(-0.7, -0.7, -0.6, -0.7), difference_p = c(0.3, 0.001, 0.2, 1e-4), alpha = 0.05),
  c("replicated", "not replicated", "inconclusive", "not replicated")
)))
check("an effect that is significant but smaller than the reference still replicates, and a larger one is never 'not replicated'", safe(identical(
  replication_outcome(estimate = c(-0.2, -1.6), adjusted_p = c(0.01, 0.2), reference = c(-0.7, -0.7), difference_p = c(0.001, 0.01), alpha = 0.05),
  c("replicated", "inconclusive")
)))
two_sources <- rbind(c(-8, -4), c(-8, 0.2), c(-8, 3.5))
check("evidence_summary agrees with evidence_class when two sources are present", safe({
  old <- evidence_class(two_sources, c(10, 5), alpha = 0.05)
  new <- evidence_summary(two_sources, c(10, 5), alpha = 0.05)
  identical(old, c("supported by both", "carried by one source", "sources disagree")) && identical(new$label, c("supported by every source", "carried by one source", "sources disagree"))
}))
three_sources <- rbind(c(-8, -3, -3), c(-8, -3, 0.2), c(-8, 0.1, 0.2), c(-8, NA, -3), c(-8, 3.5, -3), c(NA, -6, NA))
check("evidence_summary counts the supporting sources with three sources and with one missing", safe({
  out <- evidence_summary(three_sources, c(10, 5, 3), alpha = 0.05)
  identical(out$n_sources, c(3L, 3L, 3L, 2L, 3L, 1L)) && identical(out$n_supporting, c(3L, 2L, 1L, 2L, 2L, 1L)) &&
    identical(out$label, c("supported by every source", "supported by 2 of 3 sources", "carried by one source", "supported by every source", "sources disagree", "one source only"))
}))

# Follow-up of a hit combined across sources.
set.seed(34)
fu_n <- 40
fu_strata <- rep(c("a", "b", "c", "d"), each = 10)
fu_x <- rep(c(1, 0), 20)
fu_z <- cbind(1, as.numeric(fu_strata == "b"), as.numeric(fu_strata == "c"), as.numeric(fu_strata == "d"))
fu_y <- rnorm(fu_n) - 1.1 * fu_x
fu_w <- runif(fu_n, 0.5, 2)
fu_checks <- safe(followup_checks(fu_y, fu_x, fu_z, fu_strata, weights = fu_w, n_perm = 2000, seed = 9))
fu_parts <- safe(followup_statistics(fu_y, fu_x, fu_z, fu_strata, weights = fu_w, n_perm = 2000, seed = 9))
check("followup_statistics holds the pieces that followup_checks summarises", safe({
  root <- sqrt(fu_w)
  direct <- fit_ols(fu_y * root, fu_x * root, fu_z * root)
  abs(fu_parts$t - direct$t) < 1e-12 && fu_parts$df == direct$df && abs(fu_parts$estimate - fu_checks$estimate) < 1e-12 &&
    abs((1 + sum(abs(fu_parts$perm_t) >= abs(fu_parts$t))) / 2001 - fu_checks$permutation_p) < 1e-12 &&
    abs(max(2 * pt(-abs(fu_parts$loo_t), fu_parts$loo_df)) - fu_checks$loo_max_p) < 1e-12 && length(fu_parts$loo_t) == fu_n && length(fu_parts$perm_t) == 2000
}))
check("combine_followup with one source gives that source's own answers", safe({
  one <- combine_followup(list(only = fu_parts), c(only = 3))
  abs(one$z - signed_z(2 * pt(-abs(fu_parts$t), fu_parts$df), fu_parts$t)) < 1e-10 && abs(one$permutation_p - fu_checks$permutation_p) < 1e-12 &&
    one$loo_sign_stable == fu_checks$loo_sign_stable && abs(one$loo_max_p - 2 * pnorm(-min(abs(signed_z(2 * pt(-abs(fu_parts$loo_t), fu_parts$loo_df), fu_parts$loo_t))))) < 1e-12
}))
fu_noise <- safe(followup_statistics(rnorm(30), rep(c(1, 0), 15), cbind(rep(1, 30)), rep("s", 30), n_perm = 2000, seed = 10))
check("combine_followup pools sources as stouffer does, and one model can only move its own source", safe({
  both <- combine_followup(list(signal = fu_parts, noise = fu_noise), c(signal = 2, noise = 1))
  z_signal <- signed_z(2 * pt(-abs(fu_parts$t), fu_parts$df), fu_parts$t)
  z_noise <- signed_z(2 * pt(-abs(fu_noise$t), fu_noise$df), fu_noise$t)
  expected <- stouffer(cbind(z_signal, z_noise), c(2, 1))$z
  loo_signal <- signed_z(2 * pt(-abs(fu_parts$loo_t), fu_parts$loo_df), fu_parts$loo_t)
  loo_noise <- signed_z(2 * pt(-abs(fu_noise$loo_t), fu_noise$loo_df), fu_noise$loo_t)
  every <- c((2 * loo_signal + z_noise) / sqrt(5), (2 * z_signal + loo_noise) / sqrt(5))
  abs(both$z - expected) < 1e-12 && both$n_loo == 70 && abs(both$loo_max_p - max(2 * pnorm(-abs(every)))) < 1e-12 && abs(both$loo_min_abs - min(abs(every))) < 1e-12 &&
    both$permutation_p > 0 && both$permutation_p <= 1 && both$permutation_p < 0.05
}))
# Shuffles that reproduce the observed arrangement tie with it; a tie must count whatever the rounding.
tie_y <- c(0.1, 0.7, 0.2, 0.9, 0.3, 0.6, 0.35, 0.05)
tie_x <- c(1, 0, 1, 0, 1, 0, 0, 1)
tie_strata <- rep(c("a", "b"), each = 4)
tie_z <- cbind(1, as.numeric(tie_strata == "b"))
tie_expected <- safe({
  set.seed(21)
  shuffles <- vapply(seq_len(600), function(b) permute_within(tie_x, tie_strata), numeric(8))
  one_t <- function(x) abs(coef(summary(lm(tie_y ~ factor(tie_strata) + x)))["x", "t value"])
  observed <- one_t(tie_x)
  keys <- apply(shuffles, 2, paste, collapse = "")
  by_key <- vapply(unique(keys), function(key) one_t(as.numeric(strsplit(key, "")[[1]])), numeric(1))
  (1 + sum(by_key[keys] >= observed * (1 - 1e-9))) / 601
})
check("followup_checks counts every shuffle that ties with the observed arrangement when given a tie tolerance", safe({
  with_tolerance <- followup_checks(tie_y, tie_x, tie_z, tie_strata, n_perm = 600, seed = 21, tie_tolerance = 1e-9)
  set.seed(21)
  same <- sum(apply(vapply(seq_len(600), function(b) permute_within(tie_x, tie_strata), numeric(8)), 2, function(x) all(x == tie_x)))
  same > 5 && abs(with_tolerance$permutation_p - tie_expected) < 1e-12 && with_tolerance$permutation_p >= (1 + same) / 601
}))
check("without a tie tolerance followup_checks is what it was", safe({
  plain <- followup_checks(fu_y, fu_x, fu_z, fu_strata, weights = fu_w, n_perm = 2000, seed = 9)
  zero <- followup_checks(fu_y, fu_x, fu_z, fu_strata, weights = fu_w, n_perm = 2000, seed = 9, tie_tolerance = 0)
  identical(plain, zero) && identical(plain, fu_checks)
}))
check("combine_followup counts ties as followup_checks does with the tolerance", safe({
  parts <- followup_statistics(tie_y, tie_x, tie_z, tie_strata, n_perm = 600, seed = 21)
  abs(combine_followup(list(only = parts), c(only = 1))$permutation_p - tie_expected) < 1e-12
}))
# Findings of the review of 2026-10-05.
check("followup_statistics gives the same shuffles whatever the chunk size, and can skip the leave-one-out refits", safe({
  small <- followup_statistics(fu_y, fu_x, fu_z, fu_strata, weights = fu_w, n_perm = 2000, seed = 9, chunk = 300)
  quick <- followup_statistics(fu_y, fu_x, fu_z, fu_strata, weights = fu_w, n_perm = 2000, seed = 9, loo = FALSE)
  identical(small$perm_t, fu_parts$perm_t) && identical(small$loo_t, fu_parts$loo_t) && identical(quick$perm_t, fu_parts$perm_t) && identical(quick$t, fu_parts$t) && all(is.na(quick$loo_t))
}))
check("combine_followup works without leave-one-out pieces", safe({
  quick <- followup_statistics(fu_y, fu_x, fu_z, fu_strata, weights = fu_w, n_perm = 2000, seed = 9, loo = FALSE)
  out <- combine_followup(list(only = quick), c(only = 1))
  abs(out$permutation_p - fu_checks$permutation_p) < 1e-12 && is.na(out$loo_max_p) && is.na(out$loo_sign_stable)
}))
check("partner_model prefers the same name, then the donor's preferred model, then its first, and never joins a missing donor", safe(identical(
  partner_model(
    name = c("HCM-X-9-C15-B", "HCM-X-7-C18", "HCM-X-5-C20", "OTHER", "HCM-X-3-C15-C"), donor = c("HCM-X-9", "HCM-X-7", "HCM-X-5", NA, "HCM-X-3"),
    candidates = c("HCM-X-9-C15-A", "HCM-X-9-C15-B", "HCM-X-7-C18-A", "HCM-X-7-C18-B", "WTSI-1", "HCM-X-3-C15-A", "HCM-X-3-C15-B"),
    candidate_donor = c("HCM-X-9", "HCM-X-9", "HCM-X-7", "HCM-X-7", NA, "HCM-X-3", "HCM-X-3"), preferred = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
  ),
  c("HCM-X-9-C15-B", "HCM-X-7-C18-B", NA, NA, "HCM-X-3-C15-A")
)))
check("a missing source is left out of the combined follow-up", safe({
  with_gap <- combine_followup(list(signal = fu_parts, absent = NULL), c(signal = 2, absent = 1))
  alone <- combine_followup(list(signal = fu_parts), c(signal = 2))
  identical(with_gap$z, alone$z) && identical(with_gap$permutation_p, alone$permutation_p) && with_gap$n_sources == 1
}))

# --- driver status inside cancer types (stages 34 to 36) ---------------------------------
dt_models <- paste0("M", 1:6)
dt_status <- safe(gene_status(
  dt_models, c("KRAS", "MYC", "TP53"),
  events = data.table(ModelID = c("M1", "M2", "M3", "M9", "M1"), gene = c("KRAS", "KRAS", "MYC", "KRAS", "OTHER")),
  mutation_profiled = paste0("M", 1:5), copy_missing = data.table(ModelID = c("M4", "M3"), gene = c("MYC", "MYC"))
))
check("gene_status marks a gene altered with any listed event", safe(identical(unname(dt_status[1:3, "KRAS"]), c(1, 1, 0)) && dt_status["M3", "MYC"] == 1))
check("gene_status leaves a line without mutation calls missing, never unaltered", safe(all(is.na(dt_status["M6", ]))))
check("gene_status leaves a copy-number gene missing where its ratio is missing", safe(is.na(dt_status["M4", "MYC"]) && dt_status["M4", "KRAS"] == 0 && dt_status["M5", "MYC"] == 0))
check("gene_status counts an event even where another kind of data is missing", safe(dt_status["M3", "MYC"] == 1))
check("gene_status returns one row per model and one column per gene asked for", safe(identical(dimnames(dt_status), list(dt_models, c("KRAS", "MYC", "TP53")))))

# 20 lines of type A and 4 of type B; four genes with four different outcomes at a gate of five.
cs_type <- c(rep("A", 20), rep("B", 4))
cs_exposed <- c(rep(c(1, 0), each = 10), 1, 1, 0, 0)
cs_status <- cbind(
  G1 = c(rep(c(1, 0), each = 5), rep(c(1, 0), each = 5), 1, 0, 1, 0),
  G2 = c(rep(1, 6), rep(0, 2), NA, NA, rep(1, 5), rep(0, 5), 0, 0, 0, 0),
  G3 = c(rep(1, 3), rep(0, 7), rep(1, 5), rep(0, 5), 0, 0, 0, 0),
  G4 = c(rep(0, 10), rep(1, 4), rep(0, 6), 0, 0, 0, 0),
  G5 = rep(0, 24)
)
census <- safe(split_census(cs_status, cs_type, cs_exposed, min_cell = 5))
cs_row <- function(g, t = "A") census[gene == g & type == t]
check("split_census counts the four cells of a gene in a cancer type", safe(identical(unlist(cs_row("G1")[, .(neg_alt, pos_alt, neg_wt, pos_wt)], use.names = FALSE), c(5, 5, 5, 5))))
check("split_census leaves lines with a missing status out of the unaltered cells", safe(identical(unlist(cs_row("G2")[, .(neg_alt, pos_alt, neg_wt, pos_wt)], use.names = FALSE), c(6, 5, 2, 5))))
check("split_census labels each pair by the gate it reaches", safe(identical(
  c(cs_row("G1")$label, cs_row("G2")$label, cs_row("G3")$label, cs_row("G4")$label),
  c("specificity testable", "Y effect among altered only", "below the gate", "one Y group missing")
)))
check("split_census lists only pairs with an altered line", safe(nrow(census[gene == "G5"]) == 0 && nrow(cs_row("G1", "B")) == 1 && nrow(cs_row("G3", "B")) == 0))
check("same_lines_as names the first split with the same gene and the same lines", safe(identical(
  same_lines_as(c("a", "b", "c", "d"), c("K", "K", "K", "T"), list(c("M1", "M2"), c("M2", "M1"), c("M1", "M3"), c("M1", "M2"))),
  c(NA, "a", NA, NA)
)))
check("overlap_groups joins sets that share most of their members and no others", safe({
  groups <- overlap_groups(list(1:10, c(1:9, 11), 20:25, c(1:9, 12)), cut = 0.8)
  groups[1] == groups[2] && groups[1] == groups[4] && groups[3] != groups[1]
}))

set.seed(34)
sp_n <- 120
sp_source <- data.table(
  ModelID = sprintf("ACH-%06d", seq_len(sp_n)), OncotreeLineage = rep(c("Lung", "Skin"), c(80, 40)),
  exposed = rbinom(sp_n, 1, 0.45), LibraryKY = rbinom(sp_n, 1, 0.3)
)
sp_value <- rbinom(sp_n, 1, 0.4)
sp_value[c(3, 50)] <- NA
sp_frame <- safe(split_frame(sp_source, sp_source$OncotreeLineage == "Lung", sp_value))
check("split_frame keeps the lines of the type with a known status", safe(nrow(sp_frame) == 78 && all(sp_frame$OncotreeLineage == "Lung") && !anyNA(sp_frame$ctx)))
check("split_frame adds the Y columns of the model", safe(identical(sp_frame$exposed_in, sp_frame$exposed * sp_frame$ctx) && identical(sp_frame$exposed_out, sp_frame$exposed * (1 - sp_frame$ctx)) && identical(sp_frame$interaction, sp_frame$exposed_in)))
sp_genes <- paste0("G", 1:200, " (", 1:200, ")")
sp_raw <- matrix(rnorm(200 * sp_n), 200, sp_n, dimnames = list(sp_genes, sp_source$ModelID))
sp_rank <- t(apply(sp_raw, 1, rank_normalise))
dimnames(sp_rank) <- dimnames(sp_raw)
sp_fit <- safe(fit_split(sp_frame, sp_rank, "LibraryKY", min_cell = 5))
check("fit_split gives the interaction of the screen's own fit on a split with four cells", safe(isTRUE(all.equal(sp_fit$interaction_p, fit_context_tests(sp_frame, sp_rank, "LibraryKY")$interaction_p))))
sp_lm <- safe(coef(lm(sp_rank[7, sp_frame$ModelID] ~ LibraryKY + ctx + exposed_in + exposed_out, data = sp_frame)))
check("fit_split estimates equal plain lm", safe(
  abs(sp_fit$inside_estimate[7] - sp_lm[["exposed_in"]]) < 1e-10 && abs(sp_fit$outside_estimate[7] - sp_lm[["exposed_out"]]) < 1e-10 &&
    abs(sp_fit$interaction_estimate[7] - (sp_lm[["exposed_in"]] - sp_lm[["exposed_out"]])) < 1e-10
))
check("fit_split returns one row per feature with the cell counts", safe(identical(sp_fit$feature, sp_genes) && all(sp_fit$n_neg_in == sum(sp_frame$exposed == 1 & sp_frame$ctx == 1))))
sp_inside <- sp_frame[ctx == 1]
sp_alone <- safe(fit_split(sp_inside, sp_rank, "LibraryKY", min_cell = 5))
check("fit_split fits a split without unaltered lines and reports no specificity", safe({
  design <- build_design(sp_inside, tested = "exposed", covariates = "LibraryKY")
  direct <- limma_scan(sp_rank[, sp_inside$ModelID], design$design, "exposed", trend = FALSE)
  all(c("interaction_p", "outside_estimate") %in% names(sp_alone)) && isTRUE(all.equal(sp_alone$inside_p, direct$p)) && all(is.na(sp_alone$interaction_p)) && all(is.na(sp_alone$outside_estimate))
}))
# Three unaltered Y-positive lines: the Y effect among altered lines is fitted, specificity is not testable.
sp_thin <- sp_frame[!(ctx == 0 & exposed == 0) | ModelID %in% sp_frame[ctx == 0 & exposed == 0, ModelID][1:3]]
sp_thin_fit <- safe(fit_split(sp_thin, sp_rank, "LibraryKY", min_cell = 5))
check("fit_split reports specificity only when all four cells reach the gate", safe(all(is.finite(sp_thin_fit$inside_p)) && all(is.na(sp_thin_fit$interaction_p)) && all(is.finite(sp_thin_fit$outside_estimate))))
check("information_kept is 1 for a plain two-group comparison", safe(abs(information_kept(cbind(intercept = 1, exposed = rep(c(1, 0), c(4, 6))), "exposed", 4, 6) - 1) < 1e-12))
check("information_kept falls when a covariate follows Y status", safe({
  follows <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1)
  kept <- information_kept(cbind(intercept = 1, exposed = rep(c(1, 0), c(4, 6)), library = follows), "exposed", 4, 6)
  abs(kept - (1 - cor(rep(c(1, 0), c(4, 6)), follows)^2)) < 1e-12 && kept < 0.7
}))

# One more gene, tested as if it alone were added to the split's fit.
sp_target <- sp_source$ModelID %in% sp_frame[exposed == 1 & ctx == 1, ModelID]
sp_spiked <- spike_cell(sp_raw, sp_genes[1:6], sp_target, -3)
split_reference <- function(frame, rows) {
  rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    one <- matrix(rank_normalise(rows[i, ]), 1, dimnames = list("sim__1", colnames(rows)))
    suppressWarnings(fit_split(frame, rbind(sp_rank, one), "LibraryKY", min_cell = 5))[.N, .(inside_p, interaction_p)]
  }))
}
sp_parts <- safe(split_fit_parts(sp_frame, sp_rank, "LibraryKY", min_cell = 5))
sp_fast <- safe(one_more_split_tests(sp_parts, sp_spiked))
check("one_more_split_tests equals refitting the split with the gene added", safe({
  slow <- split_reference(sp_frame, sp_spiked)
  isTRUE(all.equal(sp_fast$inside_p, slow$inside_p, tolerance = 1e-8)) && isTRUE(all.equal(sp_fast$interaction_p, slow$interaction_p, tolerance = 1e-8))
}))
check("one_more_split_tests works on a split without a specificity test", safe({
  fast <- one_more_split_tests(split_fit_parts(sp_inside, sp_rank, "LibraryKY", min_cell = 5), sp_spiked)
  slow <- split_reference(sp_inside, sp_spiked)
  isTRUE(all.equal(fast$inside_p, slow$inside_p, tolerance = 1e-8)) && all(is.na(fast$interaction_p))
}))
check("a large shift in the altered Y-negative lines gives small p-values", safe(median(sp_fast$inside_p) < 1e-3 && median(sp_fast$interaction_p) < 0.05))
check("one_more_split_tests refuses rows that do not cover the screen's lines", safe(fails(one_more_split_tests(sp_parts, sp_spiked[, 1:20]))))

# Enlarged copies: each row is one copy of the design; the rows of a group are tested together.
check("enlarged_split_tests with one copy per gene equals one_more_split_tests", safe({
  single <- enlarged_split_tests(sp_parts, sp_spiked, group = seq_len(nrow(sp_spiked)))
  identical(single$inside_p, sp_fast$inside_p) && identical(single$interaction_p, sp_fast$interaction_p) && all(single$copies == 1)
}))
sp_four <- safe(enlarged_split_tests(sp_parts, sp_spiked[1:4, ], group = rep(1, 4)))
check("enlarged_split_tests equals a least-squares fit of the copies with their own nuisance terms", safe({
  ranked <- t(apply(sp_spiked[1:4, ], 1, rank_normalise))[, sp_frame$ModelID]
  stacked <- rbindlist(lapply(1:4, function(k) cbind(sp_frame[, .(LibraryKY, ctx, exposed_in, exposed_out)], y = ranked[k, ], copy = factor(k))))
  fit <- lm(y ~ 0 + copy + copy:LibraryKY + copy:ctx + copy:exposed_out + exposed_in, data = stacked)
  abs(sp_four$inside_estimate - coef(fit)[["exposed_in"]]) < 1e-10 && abs(sp_four$inside_sigma2 - summary(fit)$sigma^2) < 1e-10 && sp_four$inside_df == fit$df.residual
}))
check("the standard error of an enlarged split falls as the square root of the copies", safe(abs(sp_four$inside_unscaled - sp_parts$inside$unscaled / 2) < 1e-12 && sp_four$copies == 4))
check("enlarged_split_tests needs the same number of copies in every group", safe(fails(enlarged_split_tests(sp_parts, sp_spiked[1:5, ], group = c(1, 1, 1, 2, 2)))))

check("lines_needed is 1 when the effect is found today", safe(lines_needed(c(1, 2, 4), c(1, 0.7, 0.5), 1.5) == 1))
check("lines_needed interpolates between two multiples", safe(abs(lines_needed(c(1, 2, 4), c(4, 2, 1), 1.5) - 4 / 1.5) < 1e-10))
check("lines_needed is infinite when no multiple reaches the effect", safe(is.infinite(lines_needed(c(1, 2, 4), c(8, 6, 5), 1.5)) && is.infinite(lines_needed(c(1, 2, 4), c(NA, NA, NA), 1.5))))
check("lines_needed takes the first multiple that reaches the effect when smaller ones never do", safe(lines_needed(c(1, 2, 4), c(NA, NA, 1), 1.5) == 4))
check("lines_needed works for several effects at once", safe(isTRUE(all.equal(lines_needed(c(1, 2, 4), c(4, 2, 1), c(5, 3, 0.5)), c(1, 4 / 3, Inf)))))

# --- added after the review of the driver-by-type stages -----------------------------------
check("permutation_floor is the smallest p-value the follow-up's shuffles can give", safe(
  abs(permutation_floor(7, 6, 20000) - (1 / choose(13, 6) + 1 / 20001)) < 1e-15 && abs(permutation_floor(5, 5, 20000) - (2 / choose(10, 5) + 1 / 20001)) < 1e-15 &&
    permutation_floor(7, 6, 20000) < 0.001 && permutation_floor(5, 7, 20000) > 0.001 && identical(length(permutation_floor(c(5, 6), c(5, 8), 20000)), 2L)
))
sp_ols <- safe(fit_split_ols(sp_frame, sp_rank, "LibraryKY", min_cell = 5))
check("fit_split_ols gives the p-values of plain lm for the Y effect among altered lines and for the difference", safe({
  fit <- summary(lm(sp_rank[7, sp_frame$ModelID] ~ LibraryKY + ctx + exposed_in + exposed_out, data = sp_frame))$coefficients
  difference <- summary(lm(sp_rank[7, sp_frame$ModelID] ~ LibraryKY + ctx + exposed + interaction, data = sp_frame))$coefficients
  identical(sp_ols$feature, sp_genes) && abs(sp_ols$inside_p_ols[7] / fit["exposed_in", "Pr(>|t|)"] - 1) < 1e-8 && abs(sp_ols$interaction_p_ols[7] / difference["interaction", "Pr(>|t|)"] - 1) < 1e-8
}))
check("fit_split_ols fits a gene with missing values on its own lines and reports no difference without four cells", safe({
  holes <- sp_rank
  holes[3, sp_frame$ModelID[1:4]] <- NA
  with_holes <- fit_split_ols(sp_frame, holes, "LibraryKY", min_cell = 5)
  kept <- sp_frame[-(1:4)]
  fit <- summary(lm(holes[3, kept$ModelID] ~ LibraryKY + ctx + exposed_in + exposed_out, data = kept))$coefficients
  abs(with_holes$inside_p_ols[3] / fit["exposed_in", "Pr(>|t|)"] - 1) < 1e-8 && identical(with_holes$inside_p_ols[-3], sp_ols$inside_p_ols[-3]) && all(is.na(fit_split_ols(sp_thin, sp_rank, "LibraryKY", min_cell = 5)$interaction_p_ols))
}))
check("ols_split_tests gives plain lm p-values for values on any scale, without ranking them", safe({
  rows <- sp_raw[1:5, ] * 3 + 1
  fast <- ols_split_tests(sp_parts, rows)
  slow <- vapply(1:5, function(i) summary(lm(rows[i, sp_frame$ModelID] ~ LibraryKY + ctx + exposed_in + exposed_out, data = sp_frame))$coefficients["exposed_in", "Pr(>|t|)"], numeric(1))
  slow_difference <- vapply(1:5, function(i) summary(lm(rows[i, sp_frame$ModelID] ~ LibraryKY + ctx + exposed + interaction, data = sp_frame))$coefficients["interaction", "Pr(>|t|)"], numeric(1))
  isTRUE(all.equal(fast$inside_p, slow, tolerance = 1e-8)) && isTRUE(all.equal(fast$interaction_p, slow_difference, tolerance = 1e-8))
}))
check("the p-value of an enlarged split is limma's moderated test of the stacked gene beside today's real genes", safe({
  two <- enlarged_split_tests(sp_parts, sp_spiked[1:2, ], group = c(1, 1))
  part <- sp_parts$inside
  squeezed <- squeezeVar(c(part$sigma2, two$inside_sigma2), c(part$df_residual, two$inside_df), robust = TRUE)
  last <- length(part$sigma2) + 1
  df_prior <- if (length(squeezed$df.prior) > 1) squeezed$df.prior[last] else squeezed$df.prior
  by_hand <- 2 * pt(-abs(two$inside_estimate / two$inside_unscaled / sqrt(squeezed$var.post[last])), min(two$inside_df + df_prior, sum(part$df_residual) + two$inside_df))
  two$copies == 2 && abs(two$inside_p / by_hand - 1) < 1e-10 && two$inside_p < sp_fast$inside_p[1]
}))
check("lines_needed_bound marks a multiple that could not be interpolated", safe(identical(
  lines_needed_bound(c(1, 2, 4), c(NA, 2, 1), c(3, 1.5, 0.5)), c(TRUE, FALSE, FALSE)
) && identical(lines_needed_bound(c(1, 2, 4), c(4, 2, 1), c(5, 1.5)), c(FALSE, FALSE))))

# --- added after the rerun of the driver-by-type stages on the other processor type --------
check("tied_order keeps the given order for p-values that differ only in the last bits", safe({
  p <- c(2e-5, 1e-5 * (1 + 4e-15), 1e-5, 1e-7)
  identical(order(p), c(4L, 3L, 2L, 1L)) && identical(tied_order(p, c(1, 2, 3, 4)), c(4L, 2L, 3L, 1L)) &&
    identical(tied_order(p, c(1, 3, 2, 4)), c(4L, 3L, 2L, 1L)) && identical(tied_order(c(3e-4, 1e-4, NA), c(1, 2, 3)), c(2L, 1L, 3L)) &&
    identical(tied_order(p, c(1, 2, 3, 4), first = c("b", "a", "a", "b")), c(2L, 3L, 4L, 1L))
}))

unlink(tmp, recursive = TRUE)
cat(sprintf("%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
