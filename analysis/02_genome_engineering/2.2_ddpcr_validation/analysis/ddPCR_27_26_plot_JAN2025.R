readRDS(file = "A:/Dropbox (1)/LOY_paralogs/data/ddPCR/AMELX Craig_JAN_26_2025_2025-01-26-21-39/DATA.RDS") -> a 
readRDS(file = "A:/Dropbox (1)/LOY_paralogs/data/ddPCR/AMELX Craig_JAN_27_2025_2025-01-27-16-45/DATA.RDS") -> b 


a %>% dplyr::select(name,ratio, cell_line, run,TargetName) -> a
b %>% dplyr::select(name,ratio, cell_line, run, TargetName) -> b

rbind(a,b) -> tmp

tmp[tmp$run == "A549_Q53+CP",]$run <- "A549_CP+Q53" 

ggpubr::ggbarplot(data = tmp %>% dplyr::filter(TargetName == "AMELY/AMELX ratio" & !cell_line == "HT1080" ) , x = "run",facet.by = "TargetName", fill  = "cell_line", add = c("mean_se"),
                  y = "ratio", palette = "aaas", position = position_dodge()) -> p
ggpubr::ggpar(legend.title = "Cell line", p = p, x.text.angle = 45, ylim = c(0,1.05), xlab = "Cell line",
              ylab = "AMELY/AMELX ratio (ddPCR)", legend = "none") -> a

ggpubr::ggbarplot(data = tmp %>% dplyr::filter(TargetName == "DDX3Y/RNaseP ratio" & !cell_line == "HT1080" ) , x = "run",facet.by = "TargetName", fill  = "cell_line", add = c("mean_se"),
                  y = "ratio", palette = "aaas", position = position_dodge()) -> p
ggpubr::ggpar(legend.title = "Cell line", p = p, x.text.angle = 45, ylim = c(0,1.05), xlab = "Cell line",
              ylab = "DDX3Y/RNaseP ratio (ddPCR)", legend = "none") -> b

ggarrange(plotlist = list(a,b))



