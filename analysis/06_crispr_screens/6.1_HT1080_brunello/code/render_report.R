# Render the HT1080 CRISPR screen analysis report
# Run with: Rscript render_report.R

setwd("A:/Dropbox/LOY_paralogs/projects/crispr_screen/ht1080_loy_screen_2/code")

rmarkdown::render(
  input = "ht1080_screen2_analysis.Rmd",
  output_dir = "../results",
  output_file = "ht1080_screen2_analysis.html"
)

cat("Report rendered successfully to: ../results/ht1080_screen2_analysis.html\n")
