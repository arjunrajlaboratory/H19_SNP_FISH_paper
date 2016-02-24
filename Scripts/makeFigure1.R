#Produces 4 BXB, CXB, CXC, and BXB WT Control figures in Figure1. 

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(yaml)

source('SET_ABSOLUTE_PATH.R')


dataSelections <- c("Data_BxB_MEF_WT_Control",
                    "Data_BxC_MEF_WT_Control",
                    "Data_CxB_MEF_WT_Control",
                    "Data_CxC_MEF_WT_Control")


for (data in dataSelections){
  
  datadir <- file.path(PROJHOME, 'Raw Data/')
  dataSelect <- data
  pixelShift <- F
  outdir <- file.path(PROJHOME, 'Plots/Figure 1')
  
  
  dataSelectdir <- paste(datadir, "ExtractedCSVs/", dataSelect, sep ="")
  outSelectdir <- outdir
  yamldir <- paste(datadir, dataSelect, "/readMe.yaml", sep = "")
  expData <- yaml.load_file(yamldir)
  
  if (!pixelShift) {
    csv <- '_MasterGuideSpotsTable.csv'
    graphName23 <- dataSelect
  }else {
    csv <- '_MasterGuideSpotsTable_Shifted.csv'
    graphName23 <- paste(dataSelect, "_Shifted", sep = "")
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
  cellCounts_forPlot_Wide <- byCellCounts_wide %>% filter(B6 + C7 > minDetec)
  cellCounts_forPlot_Wide <- cellCounts_forPlot_Wide[seq(1,50),] %>%
    mutate(mono = (C7/(C7+B6)))  %>% arrange((B6+C7)) %>% mutate(plot_idx = seq(1,n())) %>% filter(plot_idx <= 50)
  
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
    theme_bw(base_size = 20, base_family = "Helvetica") +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_manual(values = c("B6" ="deepskyblue1", "C7" = "orange")) +
    ggtitle(expData$name)
  plot(upDownBar2)
  currdir <- getwd()
  setwd(outSelectdir)
  ggsave(plot = upDownBar2,  width = 9.89, height = 9,   filename =  paste(graphName23, "Figure1_Plot.pdf", sep = "_"))
  setwd(currdir)
  
}

