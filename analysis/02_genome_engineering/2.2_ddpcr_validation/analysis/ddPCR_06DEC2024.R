library(tidyverse)


## call the CH1 AMELX and the CH2 AMELY
## call all NTCs, as NTC

## path to assay dir
path = "A:/Dropbox/LOY_paralogs/data/ddPCR/AMELX Craig_2024-12-06-11-28/"

date <- "06_DEC_2024"

## pull from the plate export in QX
data.table::fread("A:/Dropbox/LOY_paralogs/data/ddPCR/AMELX Craig_2024-12-06-11-28/plate.csv") -> plate
plate %>% dplyr::distinct(Well, .keep_all = T) -> plate

## pull from the data table tab on QX
data.table::fread("a:/Dropbox/LOY_paralogs/data/ddPCR/AMELX Craig_2024-12-06-11-28/data_table.csv") -> dt
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
  ifelse(test = df$AMELX== 0 & df$AMELY == 1, yes = "AMELY", 
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

dplyr::left_join(x = l, y = plate %>% dplyr::select(Well, `Sample description 1`), by = c("well"="Well")) -> l


l %>% dplyr::filter(!(`Sample description 1` %in%  c("NTC","NTC"))) -> l
#l %>% dplyr::filter(!(well %in% c("A01","A07","B08","E09","H07"))) -> l

table(l$well,l$call)

paste0(l$well,"_",l$`Sample description 1`) -> l$name

l %>% dplyr::group_by(name, `Sample description 1`,call) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  reshape2::dcast(name ~ call, value.var = "n") -> tmp
tmp[is.na(tmp)] <- 0

stringr::str_split(string = tmp$name, pattern = "_", simplify = T)[,2] -> tmp$run
tmp %>% dplyr::arrange(run) -> tmp
tmp %>% dplyr::mutate(ratio = (AMELY + `X+Y`) / (AMELX + `X+Y`)) -> tmp

tmp %>% dplyr::mutate(empty_rate = empty / (AMELX + AMELY + `X+Y` + empty)) -> tmp
tmp %>% dplyr::mutate(double_rate = `X+Y` / (AMELX + AMELY + `X+Y` + empty)) -> tmp

tmp$date <- date


dt %>% dplyr::select(name, Target, `Conc(copies/µL)`, PoissonConfMax, PoissonConfMin ) -> dt
reshape2::dcast(data = dt, name ~ Target, value.var = "Conc(copies/µL)") -> a
names(a) <- c("name","AMELX_Conc","AMELY_Conc") 

reshape2::dcast(data = dt, name ~ Target, value.var = "PoissonConfMax") -> b 
names(b) <- c("name","AMELX_PoissonConfMax","AMELY_PoissonConfMax") 

reshape2::dcast(data = dt, name ~ Target, value.var = "PoissonConfMin") -> c
names(c) <- c("name","AMELX_PoissonConfMin","AMELY_PoissonConfMin")

dplyr::left_join(x = a, y = b, by = "name") -> a
dplyr::left_join(x = a, y = c, by = "name") -> a

dplyr::left_join(x = tmp, y = a, by = "name") -> tmp

tmp$AMELX_Conc %>% as.numeric() -> tmp$AMELX_Conc
tmp$AMELY_Conc %>% as.numeric() -> tmp$AMELY_Conc

dplyr::mutate(.data = tmp, Conc_ratio = AMELY_Conc / AMELX_Conc) -> tmp

data.table::fwrite(x = tmp, file = paste0(path,"/INFO_TABLE.csv"))


# Calculate the ratio and propagate confidence intervals
tmp <- tmp %>%
  mutate(
    Ratio_CI_low = (AMELY_PoissonConfMin / AMELX_PoissonConfMax), # Lower bound of the ratio
    Ratio_CI_high = (AMELY_PoissonConfMax / AMELX_PoissonConfMin) # Upper bound of the ratio
  )

# Plot the ratio with error bars
ggplot(tmp, aes(x = name, y = Conc_ratio)) +
  geom_point(size = 3) +  # Points for the ratio
  geom_errorbar(aes(ymin = Ratio_CI_low, ymax = Ratio_CI_high), width = 0.2) + # Error bars
  labs(
    title = "AMELY/AMELX Ratio with Confidence Intervals",
    x = "Sample",
    y = "Ratio (AMELY / AMELX)"
  ) +
  theme_pubr() -> p

ggpar(p, x.text.angle = 45)

