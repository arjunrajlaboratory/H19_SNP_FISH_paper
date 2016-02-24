library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(yaml)
library(reshape2)

source('SET_ABSOLUTE_PATH.R')

datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <-  "Data_CxB_MEF_9CG_IGF2"
pixelShift <- F
outSelectdir <- file.path(PROJHOME, 'Plots/Supp Figure 4/')
outdir <- file.path(PROJHOME, 'Plots/Supp Figure 4/')


dataSelectdir <- paste(datadir, "ExtractedCSVs/", dataSelect, sep ="")

yamldir <- paste(datadir, dataSelect, "/readMe.yaml", sep = "")
expData <- yaml.load_file(yamldir)

if (!pixelShift) {
  csv <- '_MasterGuideSpotsTable.csv'
  graphName <- dataSelect
}else {
  csv <- '_MasterGuideSpotsTable_Shifted.csv'
  graphName <- paste(dataSelect, "_Shifted", sep = "")
}

tabletoRead <- paste(datadir, "ExtractedCSVs/", expData$name, csv, sep = "")

dataTable <- tbl_df(read.csv(tabletoRead , stringsAsFactors = T))
dataTable$labels <- revalue(dataTable$labels, c("B6 H19 SNPs"="B6", "C7 H19 SNPs"="C7",
                                                "H19 B6" = "B6", "H19 CAST" = "C7", "CAST" = "C7"))
dataTable$labels <- factor(dataTable$labels,  c("B6", "C7", "undetec", "3-color"))


byCellCounts <- dataTable %>% group_by(cellID, labels) %>% summarize(rnaCount = n(), numGFP = unique(numGFP)) %>%
  mutate(IGF2 = ifelse(numGFP > 100, T, F))


totals <- byCellCounts %>% group_by(labels) %>% summarize(totalCount = sum(rnaCount))
totals$name <- as.factor(dataSelect)
totals <- totals %>% mutate(totalFrac = totalCount/sum(totalCount))

byCellCounts_wide <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
toJoin <- byCellCounts %>% group_by(cellID) %>% summarise(IGF2 = unique(IGF2))
byCellCounts_wide$cellID <- as.numeric(byCellCounts_wide$cellID)
byCellCounts_wide <- left_join(byCellCounts_wide, toJoin)


##Plotting.....

minDetec <- 15
monoFrac <- 0.8
cellCounts_forPlot_Wide <- byCellCounts_wide %>% filter(B6 + C7 > minDetec) %>% mutate(plot_idx = seq(1,n())) %>% 
  mutate(mono = (C7/(C7+B6))) %>% mutate(isMono =  mono >= monoFrac) %>% arrange(isMono, (B6+C7)) %>% mutate(plot_idx = seq(1,n())) 

IGF2Pos <- cellCounts_forPlot_Wide$plot_idx[cellCounts_forPlot_Wide$IGF2]
C7IGF2Pos <- cellCounts_forPlot_Wide$C7[cellCounts_forPlot_Wide$IGF2]
                  
upDownBar2 <- ggplot() +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= -C7, fill = "C7"), stat = "identity",
           position = "identity") +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= B6, fill = "B6"), stat = "identity",
           position = "identity") +
  geom_point(aes(x = IGF2Pos, y=  -C7IGF2Pos  -25), colour = 'green') +
  coord_flip() +
  xlab("Individual Cells") +
  ylab("RNA Count") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("B6" ="deepskyblue1", "C7" = "orange")) +
  ggtitle(expData$name)
plot(upDownBar2)

currdir <- getwd()
setwd(outdir)
ggsave(plot = upDownBar2,  width = 9.89, height = 9,   filename =  paste("9CG_IGF2_BAR.pdf", sep = "_"))
setwd(currdir)

numGFP <- dataTable %>% group_by(cellID) %>% summarise(numGFP = unique(numGFP))
cellCounts_forPlot_Wide<- left_join(cellCounts_forPlot_Wide, numGFP)

ScatterRatio <- ggplot() +
  geom_point(data = cellCounts_forPlot_Wide, aes(x = mono, y = numGFP, size = 3)) +
  geom_line(aes(x = 0.8, y =seq(0,3500,100), colour = 'red'))+
  xlab("C7 Ratio") +
  ylab("Detected IGF2 RNA Count") +
  xlim(0,1) +
  ylim(0,700) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  ggtitle(expData$name)
plot(ScatterRatio)

currdir <- getwd()
setwd(outdir)
ggsave(plot = ScatterRatio,  width = 9.89, height = 9,   filename =  paste("9CG_IGF2_Scatter.pdf", sep = "_"))
setwd(currdir)



datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <-  "Data_CxB_MEF_WT_IGF2"
pixelShift <- F
outSelectdir <- file.path(PROJHOME, 'Plots/Supp Figure 4/')



dataSelectdir <- paste(datadir, "ExtractedCSVs/", dataSelect, sep ="")

yamldir <- paste(datadir, dataSelect, "/readMe.yaml", sep = "")
expData <- yaml.load_file(yamldir)

if (!pixelShift) {
  csv <- '_MasterGuideSpotsTable.csv'
  graphName <- dataSelect
}else {
  csv <- '_MasterGuideSpotsTable_Shifted.csv'
  graphName <- paste(dataSelect, "_Shifted", sep = "")
}

tabletoRead <- paste(datadir, "ExtractedCSVs/", expData$name, csv, sep = "")

dataTable <- tbl_df(read.csv(tabletoRead , stringsAsFactors = T))
dataTable$labels <- revalue(dataTable$labels, c("B6 H19 SNPs"="B6", "C7 H19 SNPs"="C7",
                                                "H19 B6" = "B6", "H19 CAST" = "C7", "CAST" = "C7"))
dataTable$labels <- factor(dataTable$labels,  c("B6", "C7", "undetec", "3-color"))


byCellCounts <- dataTable %>% group_by(cellID, labels) %>% summarize(rnaCount = n(), numGFP = unique(numGFP)) %>%
  mutate(IGF2 = ifelse(numGFP > 100, T, F))


totals <- byCellCounts %>% group_by(labels) %>% summarize(totalCount = sum(rnaCount))
totals$name <- as.factor(dataSelect)
totals <- totals %>% mutate(totalFrac = totalCount/sum(totalCount))

byCellCounts_wide <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
toJoin <- byCellCounts %>% group_by(cellID) %>% summarise(IGF2 = unique(IGF2))
byCellCounts_wide$cellID <- as.numeric(byCellCounts_wide$cellID)
byCellCounts_wide <- left_join(byCellCounts_wide, toJoin)


##Plotting.....

minDetec <- 15
monoFrac <- 0.8
cellCounts_forPlot_Wide <- byCellCounts_wide %>% filter(B6 + C7 > minDetec) %>% mutate(plot_idx = seq(1,n())) %>% 
  mutate(mono = (C7/(C7+B6))) %>% mutate(isMono =  mono >= monoFrac) %>% arrange(isMono, (B6+C7)) %>% mutate(plot_idx = seq(1,n())) 

IGF2Pos <- cellCounts_forPlot_Wide$plot_idx[cellCounts_forPlot_Wide$IGF2]
C7IGF2Pos <- cellCounts_forPlot_Wide$C7[cellCounts_forPlot_Wide$IGF2]

upDownBar2 <- ggplot() +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= -C7, fill = "C7"), stat = "identity",
           position = "identity") +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= B6, fill = "B6"), stat = "identity",
           position = "identity") +
  geom_point(aes(x = IGF2Pos, y=  -C7IGF2Pos  -25), colour = 'green') +
  coord_flip() +
  xlab("Individual Cells") +
  ylab("RNA Count") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("B6" ="deepskyblue1", "C7" = "orange")) +
  ggtitle(expData$name)
plot(upDownBar2)

currdir <- getwd()
setwd(outdir)
ggsave(plot = upDownBar2,  width = 9.89, height = 9,   filename =  paste("WT_IGF2_Bar.pdf", sep = "_"))
setwd(currdir)

numGFP <- dataTable %>% group_by(cellID) %>% summarise(numGFP = unique(numGFP))
cellCounts_forPlot_Wide<- left_join(cellCounts_forPlot_Wide, numGFP)

ScatterRatio <- ggplot() +
  geom_point(data = cellCounts_forPlot_Wide, aes(x = mono, y = numGFP, size = 3)) +
  geom_line(aes(x = 0.8, y =seq(0,3500,100), colour = 'red') )+
  xlab("C7 Ratio") +
  ylab("Detected IGF2 RNA Count") +
  xlim(0,1) +
  ylim(0,3000) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  ggtitle(expData$name)
plot(ScatterRatio)

currdir <- getwd()
setwd(outdir)
ggsave(plot = ScatterRatio,  width = 9.89, height = 9,   filename =  paste("WT_IGF2_Scatter.pdf", sep = "_"))
setwd(currdir)

