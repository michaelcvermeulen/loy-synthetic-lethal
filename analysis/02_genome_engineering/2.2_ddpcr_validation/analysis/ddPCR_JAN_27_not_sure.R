### 

library(tidyverse)


## call the CH1 AMELX and the CH2 AMELY
## call all NTCs, as NTC

## path to assay dir
path = "A:/Dropbox (1)/LOY_paralogs/data/ddPCR/AMELX Craig_JAN_27_2025_2025-01-27-16-45/"

date <- "JAN_27_2025"

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

l %>% dplyr::group_by(name, `Sample description 1`,call) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  reshape2::dcast(name ~ call, value.var = "n") -> tmp
tmp[is.na(tmp)] <- 0

stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,c(2)] -> tmp$cell_line
stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,c(1)] -> tmp$well
paste0(stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,c(2)],"_",
       stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,c(3)]) -> tmp$run
tmp %>% dplyr::arrange(run) -> tmp
tmp %>% dplyr::mutate(ratio = (`Y+` + `X+Y`) / (`Y-` + `X+Y`)) -> tmp

tmp %>% dplyr::mutate(empty_rate = empty / (`Y+` + `Y-` + `X+Y` + empty)) -> tmp
tmp %>% dplyr::mutate(double_rate = `X+Y` / (`Y+` + `Y-` + `X+Y` + empty)) -> tmp
tmp$date <- date


##########

library(tidyverse)


## call the CH1 AMELX and the CH2 AMELY
## call all NTCs, as NTC

data.table::fread(file = "A:/Dropbox (1)/LOY_paralogs/data/DNA_RNA database - DNA_RNA BOX1.csv") -> samples
dplyr::select(samples, ID, Sample, `Cell line`, Date) -> samples

## path to assay dir
path = "A:/Dropbox (1)/LOY_paralogs/data/ddPCR/AMELX Craig_FEB_24_2025_2025-02-25-10-57/"

date <- "FEB_24_2025"

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

stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,c(2)] -> tmp$cell_line
stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,c(1)] -> tmp$well

dplyr::left_join(x = tmp, y = samples, by = c("cell_line"="ID")) -> tmp

# tmp %>% dplyr::arrange(run) -> tmp
tmp %>% dplyr::mutate(ratio = (`Y+` + `X+Y`) / (`Y-` + `X+Y`)) -> tmp

tmp %>% dplyr::mutate(empty_rate = empty / (`Y+` + `Y-` + `X+Y` + empty)) -> tmp
tmp %>% dplyr::mutate(double_rate = `X+Y` / (`Y+` + `Y-` + `X+Y` + empty)) -> tmp
tmp$date <- date

data.table::fwrite(x = tmp, 
                   file = "A:/Dropbox (1)/LOY_paralogs/data/ddPCR/AMELX Craig_FEB_24_2025_2025-02-25-10-57/INFO_TABLE.csv")

