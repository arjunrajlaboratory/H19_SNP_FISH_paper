
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(yaml)
library(reshape2)

source('SET_ABSOLUTE_PATH.R')


dataSelections <- c("Data_BxB_MEF_WT_Control",
                    "Data_BxC_MEF_WT_Control",
                    "Data_CxB_MEF_WT_Control",
                    "Data_CxC_MEF_WT_Control",
                    "Data_CXB_MEF_9CG_R",
                    "Data_CXB_MEF_9CG_CS")
for (i in c(F)){
  for (data in dataSelections){
    
    datadir <- file.path(PROJHOME, 'Raw Data/')
    dataSelect <- data
    pixelShift <- i
    outSelectdir <- file.path(PROJHOME, 'Plots/Supp Figure 2/')
    
    
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
    
    byCellCounts <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
    
    
    ##Plotting.....
    
    
    #Plot first 50 cells that have at least 50 detected RNAs
    cellCounts_forPlot_Wide <- byCellCounts %>% filter(B6 + C7 >= 50) %>% mutate(plot_idx = seq(1,n()))
    cellCounts_forPlot_Long <- cellCounts_forPlot_Wide %>% gather(labels, rnaCount, -cellID, -plot_idx)
    
    
    
    
    #Color Values
    colorVec <- c("B6" ="deepskyblue1", "C7" = "orange", "undetec" = "grey", "3-color" = "brown")
    
    
    C7Counts <- cellCounts_forPlot_Long %>% filter(labels == "C7")
    B6Counts <- cellCounts_forPlot_Long %>% filter(labels == "B6")

    #All CELLS
    AllcellCounts_forPlot <- byCellCounts %>% mutate(plot_idx = seq(1,n())) %>% gather(labels, rnaCount, -cellID, -plot_idx)
    AllC7Counts <- AllcellCounts_forPlot %>% filter(labels == "C7")
    AllB6Counts <- AllcellCounts_forPlot %>% filter(labels == "B6")
    
    
    ## Plot Pie Graph
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
      )
    
    overallPieChart <- ggplot(data = totals) +
      geom_bar(aes(x="", y = totalFrac, fill = labels), stat = "identity", width = 1) +
      coord_polar(theta="y") + blank_theme + 
      theme(axis.text.x=element_blank()) +
      geom_text(aes(x = 1.8,  y = totalFrac/3 + c(0,cumsum(totalFrac[-length(totalFrac)])),
                    label = percent(totalFrac)), size=5) +
      scale_fill_manual(values = colorVec) +
      ggtitle(expData$name)
    plot(overallPieChart)
    
    currdir <- getwd()
    setwd(outSelectdir)
    ggsave(plot = overallPieChart, filename =   paste(graphName, "OverallPie.pdf", sep = "_"))
    setwd(currdir)
    
    
    
  }
}


## PIXEL SHIFT PIE CHART
    
    datadir <- file.path(PROJHOME, 'Raw Data/')
    dataSelect <- "Data_CxB_MEF_WT_Control"
    pixelShift <- F
    outSelectdir <- file.path(PROJHOME, 'Plots/Supp Figure 2/')
    
    
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
    
    byCellCounts <- dataTable %>% xtabs( ~  cellID + labels, data = .) %>% as.data.frame() %>% spread(labels, Freq)
    
    
    ##Plotting.....
    
    
    #Plot first 50 cells that have at least 50 detected RNAs
    cellCounts_forPlot_Wide <- byCellCounts %>% filter(B6 + C7 >= 50) %>% mutate(plot_idx = seq(1,n()))
    cellCounts_forPlot_Long <- cellCounts_forPlot_Wide %>% gather(labels, rnaCount, -cellID, -plot_idx)
    
    
    
    
    #Color Values
    colorVec <- c("B6" ="deepskyblue1", "C7" = "orange", "undetec" = "grey", "3-color" = "brown")
    
    
    C7Counts <- cellCounts_forPlot_Long %>% filter(labels == "C7")
    B6Counts <- cellCounts_forPlot_Long %>% filter(labels == "B6")
    
    #All CELLS
    AllcellCounts_forPlot <- byCellCounts %>% mutate(plot_idx = seq(1,n())) %>% gather(labels, rnaCount, -cellID, -plot_idx)
    AllC7Counts <- AllcellCounts_forPlot %>% filter(labels == "C7")
    AllB6Counts <- AllcellCounts_forPlot %>% filter(labels == "B6")
    
    
    ## Plot Pie Graph
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold")
      )
    
    overallPieChart <- ggplot(data = totals) +
      geom_bar(aes(x="", y = totalFrac, fill = labels), stat = "identity", width = 1) +
      coord_polar(theta="y") + blank_theme + 
      theme(axis.text.x=element_blank()) +
      geom_text(aes(x = 1.8,  y = totalFrac/3 + c(0,cumsum(totalFrac[-length(totalFrac)])),
                    label = percent(totalFrac)), size=5) +
      scale_fill_manual(values = colorVec) +
      ggtitle(expData$name)
    plot(overallPieChart)
    
    currdir <- getwd()
    setwd(outSelectdir)
    ggsave(plot = overallPieChart, filename =   paste(graphName, "_PIXELSHIFT_OverallPie.pdf", sep = "_"))
    setwd(currdir)
    
    
    
