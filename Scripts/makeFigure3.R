#Figure 3


library(ggplot2)
library(plyr)
library(dplyr)


source('SET_ABSOLUTE_PATH.R')

datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <- "Data_CxB_9CG_MEFinCRL_IGF2"
pixelShift <- F
outdir <- file.path(PROJHOME, 'Plots/Figure 3')
outSelectdir <- file.path(PROJHOME, 'Plots/Figure 3')



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
cellArrayObj <- tbl_df(read.csv(paste(datadir, "ExtractedCSVs/", expData$name, '_CellArrayObj.csv', sep = ""), stringsAsFactors = T))

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


grouped_Positions <- c(4,6,8,9,11,13,14,17,19,26)

clusterID <- 1;
cellArrayObj$clusterID = 1;
for (i in seq(2,nrow(cellArrayObj))){
  if ((cellArrayObj$ArrayNum[i-1] == cellArrayObj$ArrayNum[i]) || (cellArrayObj$ArrayNum[i] %in% grouped_Positions)) {
    cellArrayObj$clusterID[i] <- clusterID    
  }
  else{
  clusterID <- clusterID + 1
  cellArrayObj$clusterID[i] <- clusterID 
  }
}


cellCounts_forPlot_Wide <- byCellCounts_wide %>% left_join(cellArrayObj)

#FILTER COLONIES that are high expressing with 3-4 cells
filterID <- c(5:7, 10:13, 24:28, 29:31, 35:39, 46:50, 57:60)
cellCounts_forPlot_Wide <- cellCounts_forPlot_Wide %>% filter(seq(1,nrow(cellCounts_forPlot_Wide)) %in% filterID)


nullRow <- as.data.frame(t(rep(0,9)))
nullRow <- rbind(nullRow, nullRow, nullRow, nullRow)
colnames(nullRow) <- colnames(cellCounts_forPlot_Wide)


newCellCounts_forPlot_Wide <- cellCounts_forPlot_Wide[1,]

for (i in seq(2, nrow(cellCounts_forPlot_Wide))){
  if (cellCounts_forPlot_Wide$clusterID[i] != cellCounts_forPlot_Wide$clusterID[i-1]){
    newCellCounts_forPlot_Wide <- rbind(newCellCounts_forPlot_Wide, nullRow)
  }
  newCellCounts_forPlot_Wide <- rbind(newCellCounts_forPlot_Wide, cellCounts_forPlot_Wide[i,])
}



#Color Values
colorVec <- c("B6" ="deepskyblue1", "C7" = "orange", "undetec" = "grey", "3-color" = "brown")


newCellCounts_forPlot_Wide$numRow = 1:nrow(newCellCounts_forPlot_Wide)

IGF2Plot <- newCellCounts_forPlot_Wide %>% filter(IGF2 == 1) 

upDownBar2 <- ggplot() +
  geom_bar(data = newCellCounts_forPlot_Wide, aes(x = 1:nrow(newCellCounts_forPlot_Wide), y = -C7,
                                                  fill = "C7"), stat = "identity", position = "dodge") +
  geom_bar(data = newCellCounts_forPlot_Wide, aes(x = 1:nrow(newCellCounts_forPlot_Wide), y = B6, fill = "B6"), stat = "identity") +
  geom_point(data = IGF2Plot, aes(x = numRow, y = (-C7 - 25)), colour = 'green' ) +
  coord_flip() +
  xlab("Individual Cells") +
  ylab("RNA Count") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("B6" ="deepskyblue1", "C7" = "orange"), name = 'Allele') +
  ggtitle(expData$name)
plot(upDownBar2)
currdir <- getwd()
setwd(outdir)
ggsave(plot = upDownBar2,  width = 9.89, height = 9,   filename =  paste("Figure3_IGF2_COLONY_PLOT.pdf", sep = "_"))
setwd(currdir)


### COMPARISON SCRIPTS

datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <- "IGF2+ MEF Colony CS"
pixelShift <- F
outSelectdir <- file.path(PROJHOME, 'Plots/Figure 3')


dataSelectdir <- paste(datadir, "ExtractedCSVs/", dataSelect, sep ="")

yamldir <- paste(datadir, "Clone Bisulfite/07-27-2015_IGF2+_1-3_4x/readMe.yaml", sep = "")
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
dataTable$labels <- revalue(dataTable$labels, c("H19 B6"="B6", "H19 CAST"="C7"))
dataTable$labels <- factor(dataTable$labels,  c("B6", "C7", "undetec", "3-color"))


byCellCounts <- dataTable %>% group_by(cellID, labels) %>% summarize(rnaCount = n())

totals <- byCellCounts %>% group_by(labels) %>% summarize(totalCount = sum(rnaCount))
totals$name <- as.factor(dataSelect)
totals <- totals %>% mutate(totalFrac = totalCount/sum(totalCount))

byCellCounts_wide <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
byCellCounts_wide$cellID <- as.numeric(byCellCounts_wide$cellID)

byCellCounts_wide <- byCellCounts_wide[1:30,]



##Plotting.....

minDetec <- 15
monoRatio <- 0.8
cellCounts_forPlot_Wide <- byCellCounts_wide %>% filter(B6 + C7 > minDetec) %>% 
  mutate(mono = (C7/(C7+B6))) %>% mutate(isMono = mono >= 0.8) %>% arrange((B6+C7)) %>% mutate(plot_idx = seq(1,n())) 

fracMono = sum(cellCounts_forPlot_Wide$isMono)/nrow(cellCounts_forPlot_Wide);
fracBi = 1 - fracMono;
#Color Values
colorVec <- c("B6" ="deepskyblue1", "C7" = "orange", "undetec" = "grey", "3-color" = "brown")

upDownBar2 <- ggplot() +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= -C7, fill = "C7"), stat = "identity",
           position = "identity") +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= B6, fill = "B6"), stat = "identity",
           position = "identity") +
  coord_flip() +
  xlab("Individual Cells") +
  ylab("RNA Count") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("B6" ="deepskyblue1", "C7" = "orange")) +
  ggtitle(graphName)
plot(upDownBar2)
currdir <- getwd()
setwd(outSelectdir)
ggsave(plot = upDownBar2, filename =   paste("MONO_UPDOWNBAR.pdf", sep = "_"))
setwd(currdir)




## Biallelic Colony


datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <- "Bi MEF Colony R"
pixelShift <- F
outSelectdir <- file.path(PROJHOME, 'Plots/Figure 3')


dataSelectdir <- paste(datadir, "ExtractedCSVs/", dataSelect, sep ="")

yamldir <- paste(datadir, "Clone Bisulfite/06-04-2015_Bi_5-3_1x_Bi#1_Fig4/readMe.yaml", sep = "")
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
dataTable$labels <- revalue(dataTable$labels, c("H19 B6"="B6", "H19 CAST"="C7"))
dataTable$labels <- factor(dataTable$labels,  c("B6", "C7", "undetec", "3-color"))

byCellCounts <- dataTable %>% group_by(cellID, labels) %>% summarize(rnaCount = n())

totals <- byCellCounts %>% group_by(labels) %>% summarize(totalCount = sum(rnaCount))
totals$name <- as.factor(dataSelect)
totals <- totals %>% mutate(totalFrac = totalCount/sum(totalCount))

byCellCounts_wide <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
byCellCounts_wide$cellID <- as.numeric(byCellCounts_wide$cellID)

byCellCounts_wide <- byCellCounts_wide[1:30,]



##Plotting.....

minDetec <- 15
monoRatio <- 0.8
cellCounts_forPlot_Wide <- byCellCounts_wide %>% filter(B6 + C7 > minDetec) %>% 
  mutate(mono = (C7/(C7+B6))) %>% mutate(isMono = mono >= 0.8) %>% arrange((B6+C7)) %>% mutate(plot_idx = seq(1,n())) 

fracMono = sum(cellCounts_forPlot_Wide$isMono)/nrow(cellCounts_forPlot_Wide);
fracBi = 1 - fracMono;
#Color Values
colorVec <- c("B6" ="deepskyblue1", "C7" = "orange", "undetec" = "grey", "3-color" = "brown")

upDownBar2 <- ggplot() +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= -C7, fill = "C7"), stat = "identity",
           position = "identity") +
  geom_bar(data = cellCounts_forPlot_Wide, aes(x = plot_idx, y= B6, fill = "B6"), stat = "identity",
           position = "identity") +
  coord_flip() +
  xlab("Individual Cells") +
  ylab("RNA Count") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("B6" ="deepskyblue1", "C7" = "orange")) +
  ggtitle(graphName)
plot(upDownBar2)
currdir <- getwd()
setwd(outSelectdir)
ggsave(plot = upDownBar2, filename =   paste("BI_Colony_UPDOWNBAR.pdf", sep = "_"))
setwd(currdir)






