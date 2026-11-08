# Re-render saved plot grobs under color-vision simulations; primary figure files are untouched.
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
setwd(dirname(normalizePath(script,winslash='/')))
source('R/05_figures.R')
ggplot2::theme_set(nature_theme())
dir.create('provenance/figure_QA',recursive=TRUE,showWarnings=FALSE)
transform_colors <- function(v, fun) {
  keep <- which(!is.na(v) & nzchar(v))
  if (length(keep)) {
    rgba <- grDevices::col2rgb(v[keep], alpha=TRUE) / 255
    hex <- grDevices::rgb(rgba[1,],rgba[2,],rgba[3,],rgba[4,])
    v[keep] <- fun(hex)
  }
  v
}
recolor <- function(g, fun) {
  # ggrepel creates label grobs during drawing, so transform its stored aesthetics too.
  if (!is.null(g$data) && is.data.frame(g$data)) {
    for (key in intersect(c('colour','color','fill','segment.colour'),names(g$data)))
      if (is.character(g$data[[key]])) g$data[[key]] <- transform_colors(g$data[[key]],fun)
  }
  for (key in c('col','fill')) if (!is.null(g$gp[[key]]) && is.character(g$gp[[key]])) {
    g$gp[[key]] <- transform_colors(g$gp[[key]],fun)
  }
  if (!is.null(g$raster) && is.character(g$raster)) {
    dims <- dim(g$raster); g$raster <- matrix(transform_colors(as.character(g$raster),fun), nrow=dims[1], ncol=dims[2])
  }
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs,recolor,fun=fun)
  if (!is.null(g$children)) g$children <- do.call(grid::gList,lapply(g$children,recolor,fun=fun))
  g
}
manifest <- data.table::fread('outputs/figures/figure_manifest.csv')
figs <- unique(manifest[,c('figure','height_mm'),with=FALSE])
requested <- commandArgs(trailingOnly = TRUE)
if (length(requested)) {
  stopifnot(all(requested %in% figs$figure))
  figs <- figs[figure %in% requested]
}
for (i in seq_len(nrow(figs))) {
  id <- figs$figure[i]; message('Color review: ',id)
  plot <- readRDS(file.path('outputs/models',paste0(id,'_plot.rds')))
  # A real device supplies the same Arial metrics as figure export.
  ragg::agg_png(file.path('provenance/figure_QA',paste0(id,'_simulations.png')),
                width=183*3,height=figs$height_mm[i]+8,units='mm',res=100,background='white')
  grob <- patchwork::patchworkGrob(plot)
  grid::grid.newpage()
  modes <- list(protanopia=colorspace::protan,deuteranopia=colorspace::deutan,grayscale=colorspace::desaturate)
  for (j in seq_along(modes)) {
    grid::pushViewport(grid::viewport(x=(j-0.5)/3,y=0.49,width=1/3,height=0.94))
    grid::grid.draw(recolor(grob,modes[[j]])); grid::popViewport()
    grid::grid.text(names(modes)[j],x=(j-0.5)/3,y=0.98,gp=grid::gpar(fontfamily='Arial',fontsize=9))
  }
  grDevices::dev.off()
}
