library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(yaml)
library(reshape2)

source('SET_ABSOLUTE_PATH.R')


datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <- "Data_CxB_MEF_9CG_R"
pixelShift <- F
outSelectdir <- file.path(PROJHOME, 'Plots/Supp Figure 3/')


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

byCellCounts <- dataTable %>% group_by(cellID, labels) %>% summarize(rnaCount = n())

totals <- byCellCounts %>% group_by(labels) %>% summarize(totalCount = sum(rnaCount))
totals$name <- as.factor(dataSelect)
totals <- totals %>% mutate(totalFrac = totalCount/sum(totalCount))

byCellCounts_wide <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
byCellCounts_wide$cellID <- as.numeric(byCellCounts_wide$cellID)



##Plotting.....

minDetec <- 15
monoRatio <- 0.8
cellCounts_forPlot_Wide <- byCellCounts_wide %>% filter(B6 + C7 > minDetec) %>% 
  mutate(mono = (C7/(C7+B6))) %>% mutate(isMono = mono >= 0.8) %>% arrange(isMono, (B6+C7)) %>% mutate(plot_idx = seq(1,n())) 



totalH19 <- dataTable %>% group_by(cellID) %>% summarise(totalH19 = n())
cellCounts_forPlot_Wide<- left_join(cellCounts_forPlot_Wide, totalH19)

ScatterRatio <- ggplot() +
  geom_point(data = cellCounts_forPlot_Wide, aes(x = mono, y = totalH19, size = 3, colour = "blue")) +
  geom_line(aes(x = 0.8, y =seq(0,4000,100), colour = 'red') )+
  xlab("C7 Ratio") +
  ylab("Total H19 RNA Count") +
  xlim(0,1) +
  ylim(0,4000) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  ggtitle(expData$name)
plot(ScatterRatio)

currdir <- getwd()
setwd(outSelectdir)
ggsave(plot = ScatterRatio,  width = 9.89, height = 9,   filename =  paste("9CG_H19_Threshold_Scatter.pdf", sep = "_"))
setwd(currdir)


meanExpressionPlot <- ggplot() +
  geom_boxplot(data = cellCounts_forPlot_Wide, aes(x = isMono, y = totalH19)) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  ylab("Total H19 RNA Count") 
plot(meanExpressionPlot) 

currdir <- getwd()
setwd(outSelectdir)
ggsave(plot = meanExpressionPlot,  width = 5.89, height = 9,   filename =  paste("Mono_vs_Bi_CellCounts.pdf", sep = "_"))
setwd(currdir)


##Calculate P value

result <- t.test(cellCounts_forPlot_Wide$totalH19 ~ cellCounts_forPlot_Wide$isMono)
result

diff <- 0;
for(runID in seq(1,10000)){
  randPerm <- sample(cellCounts_forPlot_Wide$isMono)
  meanMono <- mean(cellCounts_forPlot_Wide$totalH19[randPerm]);
  meanBi <- mean(cellCounts_forPlot_Wide$totalH19[!randPerm]);
  diff[runID] <- meanBi - meanMono
}

cumDist <- ecdf(diff)
pvalue <- 1 -cumDist(277.4)
summary(cellCounts_forPlot_Wide$isMono)

## WT SCATTER ########
datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <- "Data_CxB_MEF_WT_IGF2_2"
pixelShift <- F
outSelectdir <- file.path(PROJHOME, 'Plots/Supp Figure 3/')


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

byCellCounts <- dataTable %>% group_by(cellID, labels) %>% summarize(rnaCount = n())

totals <- byCellCounts %>% group_by(labels) %>% summarize(totalCount = sum(rnaCount))
totals$name <- as.factor(dataSelect)
totals <- totals %>% mutate(totalFrac = totalCount/sum(totalCount))

byCellCounts_wide <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
byCellCounts_wide$cellID <- as.numeric(byCellCounts_wide$cellID)
byCellCounts_wide <- filter(byCellCounts_wide, cellID < 113) #Photobleaching in final positions- throw out cells


##Plotting.....

minDetec <- 15
monoRatio <- 0.8
cellCounts_forPlot_Wide <- byCellCounts_wide %>% filter(B6 + C7 > minDetec) %>% 
  mutate(mono = (C7/(C7+B6))) %>% mutate(isMono = mono >= 0.8) %>% arrange(isMono, (B6+C7)) %>% mutate(plot_idx = seq(1,n())) 


totalH19 <- dataTable %>% group_by(cellID) %>% summarise(totalH19 = n())
cellCounts_forPlot_Wide<- left_join(cellCounts_forPlot_Wide, totalH19)

ScatterRatio <- ggplot() +
  geom_point(data = cellCounts_forPlot_Wide, aes(x = mono, y = totalH19, size = 3)) +
  geom_line(aes(x = 0.8, y =seq(0,4000,100), colour = 'red'))+
  xlab("C7 Ratio") +
  ylab("Total H19 RNA Count") +
  xlim(0,1) +
  ylim(0,3000) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  ggtitle(expData$name)
plot(ScatterRatio)

currdir <- getwd()
setwd(outSelectdir)
ggsave(plot = ScatterRatio,  width = 9.89, height = 9,   filename =  paste("WT_H19_Threshold_Scatter.pdf", sep = "_"))
setwd(currdir)

plot(ScatterRatio)


