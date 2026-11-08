# Stage 01: confirm that every source table is the file the collection documents.
# - SHA-256 of each input against docs/SHA256SUMS (a local relocation check, not publisher verification)
# - row and column counts against docs/table_summary.csv
# Stops on any mismatch and writes the manifest that later stages cite.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))

listed <- readLines(file.path(root_dir, "docs", "SHA256SUMS"))
listed <- data.table(listed_sha256 = substr(listed, 1, 64), path = substring(listed, 67))
summary_tab <- fread(file.path(root_dir, "docs", "table_summary.csv"), select = c("file", "rows", "columns"))
setnames(summary_tab, c("path", "expected_rows", "expected_columns"))

manifest <- data.table(input = names(input_paths), path = unname(input_paths))
manifest[, exists := file.exists(file.path(root_dir, path))]
if (!all(manifest$exists)) stop("missing inputs: ", paste(manifest[exists == FALSE, path], collapse = ", "))
manifest[, bytes := file.size(file.path(root_dir, path))]

log_msg("hashing %d inputs (%.1f GiB)", nrow(manifest), sum(manifest$bytes) / 2^30)
manifest[, sha256 := unlist(parallel::mclapply(
  file.path(root_dir, path),
  function(p) sub(" .*$", "", system2("sha256sum", shQuote(p), stdout = TRUE)),
  mc.cores = n_cores
))]
manifest <- merge(manifest, listed, by = "path", all.x = TRUE, sort = FALSE)
manifest[, checksum_status := fifelse(is.na(listed_sha256), "not listed", fifelse(sha256 == listed_sha256, "match", "MISMATCH"))]

# Dimensions: data rows (header excluded) and total columns, as table_summary.csv counts them.
count_dims <- function(p) {
  header <- read_header(p)
  rows <- nrow(fread(p, select = 1L, showProgress = FALSE))
  c(rows = rows, columns = length(header))
}
manifest <- merge(manifest, summary_tab, by = "path", all.x = TRUE, sort = FALSE)
dims <- parallel::mclapply(file.path(root_dir, manifest$path), function(p) {
  if (grepl("\\.csv$", p)) count_dims(p) else c(rows = NA_real_, columns = NA_real_)
}, mc.cores = n_cores)
manifest[, rows := vapply(dims, function(d) d[["rows"]], numeric(1))]
manifest[, columns := vapply(dims, function(d) d[["columns"]], numeric(1))]
manifest[, dims_status := fifelse(
  is.na(expected_rows), "not listed",
  fifelse(rows == expected_rows & columns == expected_columns, "match", "MISMATCH")
)]

write_result(manifest, "input_manifest")
print(manifest[, .(input, checksum_status, dims_status, rows, columns)])

if (any(manifest$checksum_status == "MISMATCH")) stop("checksum mismatch: ", paste(manifest[checksum_status == "MISMATCH", path], collapse = ", "))
if (any(manifest$dims_status == "MISMATCH")) stop("dimension mismatch: ", paste(manifest[dims_status == "MISMATCH", path], collapse = ", "))
log_msg(
  "inputs verified: %d checksums match, %d not listed; %d dimension checks match",
  sum(manifest$checksum_status == "match"), sum(manifest$checksum_status == "not listed"),
  sum(manifest$dims_status == "match")
)
