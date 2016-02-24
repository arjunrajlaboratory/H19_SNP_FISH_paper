
library(ggplot2)
library(plyr)
library(dplyr)



datadir <- "/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Raw Data/"
dataSelect <- "Data_CxB_9CG_MEFinCRL"
pixelShift <- F
outSelectdir <- "/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Plots/Figure 3/"


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
cellArrayObj <- tbl_df(read.csv(paste(datadir, "ExtractedCSVs/", expData$name, '_CellArrayObj.csv', sep = ""), stringsAsFactors = T))

dataTable$labels <- revalue(dataTable$labels, c("B6 H19 SNP"="B6", "C7 H19 SNP"="C7",
                                                "H19 B6 SNP" = "B6", "H19 CAST SNP" = "C7", "CAST" = "C7"))
dataTable$labels <- factor(dataTable$labels,  c("B6", "C7", "undetec", "3-color"))

byCellCounts <- dataTable %>% group_by(cellID, labels) %>% summarize(rnaCount = n())


totals <- byCellCounts %>% group_by(labels) %>% summarize(totalCount = sum(rnaCount))
totals$name <- as.factor(dataSelect)
totals <- totals %>% mutate(totalFrac = totalCount/sum(totalCount))

byCellCounts_wide <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
byCellCounts_wide$cellID <- as.numeric(byCellCounts_wide$cellID)


grouped_Positions <- c(8,10,13,16)

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
