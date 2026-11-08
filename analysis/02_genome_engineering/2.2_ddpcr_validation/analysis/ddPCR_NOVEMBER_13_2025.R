library(tidyverse)


## path to assay dir
path = "A:/Dropbox (1)/LOY_paralogs/data/ddPCR/AMELX Craig_NOV_13_2025-11-12-22-50/"

date <- "NOV_13_2025"

## pull from the amplitude export in QX
data.table::fread(paste0(path,"/plate.csv")) -> plate
plate %>% dplyr::distinct(Well, .keep_all = T) -> plate

## pull from the data table tab on QX
data.table::fread(paste0(path,"/data_table.csv")) -> dt
paste0(dt$Well,"_",dt$`Sample description 1`) -> dt$name

##############
list.files(path = path, 
           full.names = T, pattern = "*_Amplitude.csv") -> files

lapply(X = files, FUN = function(y){
  data.table::fread(file = y, skip = 4) %>% dplyr::select(-V5) %>% return()
}) -> l
names(l) <- basename(files)


lapply(X = 1:length(l), FUN = function(x){
  l[[x]] -> df
  names(df)[c(3,4)] <- c("Y","X")
  ifelse(test = df$X== 0 & df$Y== 1, yes = "Y+", 
         no = ifelse(test = df$X == 1 & df$Y == 0, yes = "Y-", 
                     no = ifelse(test = df$X == 1 & df$Y == 1, yes = "X+Y", no = "empty"))) -> df$call
  basename(files)[x] %>% stringr::str_split(pattern = "_", simplify = T) -> tmp
  tmp[,5] -> df$well
  return(df)
}) -> l
data.table::rbindlist(l) -> l
l %>% tidyr::drop_na(call) -> l

l$Ch1Amplitude %>% as.numeric() -> l$Ch1Amplitude
l$Ch2Amplitude %>% as.numeric() -> l$Ch2Amplitude
l$Y %>% as.numeric() -> l$Y
l$X %>% as.numeric() -> l$X

dplyr::left_join(x = l, y = plate %>% dplyr::select(Well, `Sample description 1`,TargetName), by = c("well"="Well")) -> l


l %>% dplyr::filter(!(`Sample description 1` %in%  c("NTC"))) -> l
#l %>% dplyr::filter(!(well %in% c("A01","A07","B08","E09","H07"))) -> l

table(l$well,l$call)

paste0(l$well,"_",l$`Sample description 1`) -> l$name

l %>% dplyr::group_by(name, `Sample description 1`,call, TargetName) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  reshape2::dcast(name + TargetName ~ call, value.var = "n") -> tmp
tmp[is.na(tmp)] <- 0

# Extract second field
stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,2] -> second_field

# Classify: if starts with P1/P2/P3/P4/D0 = HT1080, otherwise use the field itself
tmp$cell_line <- ifelse(grepl("^(P[1-4]|D0)", second_field), 
                        "HT1080", 
                        second_field)
stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,c(1)] -> tmp$well

#dplyr::left_join(x = tmp, y = samples, by = c("cell_line"="ID")) -> tmp

# tmp %>% dplyr::arrange(run) -> tmp
tmp %>% dplyr::mutate(ratio = (`Y+` + `X+Y`) / (`Y-` + `X+Y`)) -> tmp

tmp %>% dplyr::mutate(empty_rate = empty / (`Y+` + `Y-` + `X+Y` + empty)) -> tmp
tmp %>% dplyr::mutate(double_rate = `X+Y` / (`Y+` + `Y-` + `X+Y` + empty)) -> tmp
tmp$date <- date

dt %>% dplyr::select(name, Target, `Conc(copies/µL)`, PoissonConfMax, PoissonConfMin ) -> dt
reshape2::dcast(data = dt, name ~ Target, value.var = "Conc(copies/µL)") -> a
names(a) <- c("name","AMELX_Conc","AMELY_Conc","DDX3Y_Conc", "RNaseP_Conc") 

reshape2::dcast(data = dt, name ~ Target, value.var = "PoissonConfMax") -> b 
names(b) <- c("name","AMELX_PoissonConfMax","AMELY_PoissonConfMax","DDX3Y_PoissonConfMax","RNaseP_PoissonConfMax") 

reshape2::dcast(data = dt, name ~ Target, value.var = "PoissonConfMin") -> c
names(c) <- c("name","AMELX_PoissonConfMin","AMELY_PoissonConfMin","DDX3Y_PoissonConfMin","RNaseP_PoissonConfMin")

dplyr::left_join(x = a, y = b, by = "name") -> a
dplyr::left_join(x = a, y = c, by = "name") -> a

dplyr::left_join(x = tmp, y = a, by = "name") -> tmp

tmp$AMELX_Conc %>% as.numeric() -> tmp$AMELX_Conc
tmp$AMELY_Conc %>% as.numeric() -> tmp$AMELY_Conc

dplyr::mutate(.data = tmp, Conc_ratio = AMELY_Conc / AMELX_Conc) -> tmp

data.table::fwrite(x = tmp, file = paste0(path,"/INFO_TABLE.csv"))

################################
tmp$date <- date

dplyr::left_join(x = tmp %>% dplyr::select(-TargetName), y = plate %>% dplyr::select(Well, TargetName), by = c("well"="Well")) -> l

l[l$TargetName == "AMELY",]$TargetName <- "AMELY/AMELX ratio"
l[l$TargetName == "DDX3Y",]$TargetName <- "DDX3Y/RNaseP ratio"
paste(stringr::str_split(string = l$name, pattern = "_", simplify = T)[,c(2)], 
      stringr::str_split(string = l$name, pattern = "_", simplify = T)[,c(3)]) -> l$name
paste0(l$name,"_",l$Date) -> l$lab
ggpubr::ggbarplot(data = l, x = "lab",facet.by = "TargetName", fill  = "cell_line", add = c("mean_se"),
                  y = "ratio", palette = "aaas", position = position_dodge()) -> p
ggpubr::ggpar(legend.title = "Cell line", p = p, x.text.angle = 45, font.x = c(6), ylim = c(0,1.05)) -> p



ggdotplot(tmp, x = "TargetName", y = "ratio", add = "mean_se", facet.by = "TargetName")

saveRDS(object = l, file = paste0(path,"/DATA.RDS"))