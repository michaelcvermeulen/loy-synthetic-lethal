library(tidyverse)


list.files(path = "A:/Dropbox/LOY_paralogs/data/ddPCR/control_test/attempt3/raw_files/", full.names = T) -> files

lapply(X = files, FUN = function(y){
  data.table::fread(file = y, skip = 4) %>% dplyr::select(-V5) %>% return()
}) -> l
names(l) <- basename(files)


lapply(X = 1:length(l), FUN = function(x){
  l[[x]] -> df
  ifelse(test = df$AMELX == 0 & df$AMELY == 1, yes = "AMELY", 
         no = ifelse(test = df$AMELX == 1 & df$AMELY == 0, yes = "AMELX", 
                     no = ifelse(test = df$AMELX == 1 & df$AMELY == 1, yes = "X+Y", no = "empty"))) -> df$call
  basename(files)[x] %>% stringr::str_split(pattern = "_", simplify = T) -> tmp
  tmp[,2] -> df$well
  return(df)
}) -> l
data.table::rbindlist(l) -> l
l %>% tidyr::drop_na(call) -> l

l$Ch1Amplitude %>% as.numeric() -> l$Ch1Amplitude
l$Ch2Amplitude %>% as.numeric() -> l$Ch2Amplitude
l$AMELX %>% as.numeric() -> l$AMELX
l$AMELY %>% as.numeric() -> l$AMELY

data.table::fread("A:/Dropbox/LOY_paralogs/data/ddPCR/control_test/attempt3/AMELX 07May2024.csv") -> plate

plate %>% dplyr::distinct(Well, .keep_all = T) -> plate

dplyr::left_join(x = l, y = plate %>% dplyr::select(Well, `Sample description 1`), by = c("well"="Well")) -> l

l %>% dplyr::filter(!(well %in% c("B07","A07","B08","E09","H07"))) -> l
l %>% dplyr::filter(`Sample description 1` != "NTC") -> l

table(l$well,l$`Sample description 1`)

paste0(l$well,"_",l$`Sample description 1`) -> l$name

l %>% dplyr::group_by(name, `Sample description 1`,call) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  reshape2::dcast(name ~ call, value.var = "n") -> tmp
tmp[is.na(tmp)] <- 0

stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,2] -> tmp$run
tmp %>% dplyr::arrange(run) -> tmp
tmp %>% dplyr::mutate(ratio = (AMELY + `X+Y`) / AMELX) -> tmp

data.table::fread(input = "A:/Dropbox/LOY_paralogs/data/ddPCR/control_test/attempt3/ratios.csv") -> df

df %>% reshape2::melt(id = c("Expected","Conc")) -> df
ggpubr::ggscatter(data = df, y = "value", x = "Expected", add.params = list(color = "red", size = 0.6), color = "variable",
                  facet.by = "Conc", add = "reg.line", ylab = "Observed", palette = "aaas") + ggpubr::stat_cor() -> p
ggpar(p, legend = "bottom", legend.title = "Technical replicate")

ggpubr::ggscatter(l, y = "Ch1Amplitude", x = "Ch2Amplitude", color = "call", palette = "aaas", 
                  size = 0.3, alpha = 0.5) -> p
facet(p = p, facet.by = "well", ncol = 6, nrow = 5) -> p
ggpubr::ggpar(p, ylim = c(0,15000), xlim = c(0,7500), xlab = "VIC Amplitude (AMELX)", ylab = "FAM Amplitude (AMELY)", legend = "none") 

ggpubr::ggscatter(l %>% dplyr::filter(well == "F05"), y = "Ch1Amplitude", x = "Ch2Amplitude", color = "call", palette = "aaas",
                  size = 0.3, alpha = 0.5) -> p
facet(p = p, facet.by = "well", nrow = 5, ncol = 2) -> p
ggpubr::ggpar(p, ylim = c(0,15000), xlim = c(0,7500), xlab = "VIC Amplitude (AMELX)", ylab = "FAM Amplitude (AMELY)", legend = "none") 



ggarrange(plotlist = ll, ncol = 3, nrow = 8) -> p

ifelse(test = df$AMELX == 0 & df$AMELY == 1, yes = "XY", no = ifelse(test = df$AMELX == 1 & df$AMELY == 0, yes = "XO", no = "empty")) -> df$call

ggpubr::ggscatter(df, y = "Ch1Amplitude", x = "Ch2Amplitude", color = "call")