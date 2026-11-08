list.files(path = "A:/Dropbox/LOY_paralogs/data/ddPCR/control_test/attempt2/raw_files/", full.names = T) -> files

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
  tmp[,3] -> df$well
  return(df)
}) -> l
data.table::rbindlist(l) -> l
  

ggpubr::ggscatter(l, y = "Ch1Amplitude", x = "Ch2Amplitude", color = "call", palette = "aaas",
                  size = 0.3, alpha = 0.5) -> p
facet(p = p, facet.by = "well", nrow = 5, ncol = 2) -> p
ggpubr::ggpar(p, ylim = c(0,15000), xlim = c(0,7500), xlab = "VIC Amplitude (AMELX)", ylab = "FAM Amplitude (AMELY)", legend = "none") 

ggpubr::ggscatter(l %>% dplyr::filter(well == "F05"), y = "Ch1Amplitude", x = "Ch2Amplitude", color = "call", palette = "aaas",
                  size = 0.3, alpha = 0.5) -> p
facet(p = p, facet.by = "well", nrow = 5, ncol = 2) -> p
ggpubr::ggpar(p, ylim = c(0,15000), xlim = c(0,7500), xlab = "VIC Amplitude (AMELX)", ylab = "FAM Amplitude (AMELY)", legend = "none") 



ggarrange(plotlist = ll, ncol = 2, nrow = 5) -> p

ifelse(test = df$AMELX == 0 & df$AMELY == 1, yes = "XY", no = ifelse(test = df$AMELX == 1 & df$AMELY == 0, yes = "XO", no = "empty")) -> df$call

ggpubr::ggscatter(df, y = "Ch1Amplitude", x = "Ch2Amplitude", color = "call")