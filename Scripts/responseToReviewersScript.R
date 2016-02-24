## Make Response to Reviewer 1 Figufres

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(yaml)
library(reshape2)

source('SET_ABSOLUTE_PATH.R')


### PIXELSHIFT Response



datadir <- file.path(PROJHOME, 'Raw Data/')
dataSelect <- "Data_CxB_MEF_9CG_R"
pixelShift <- T
outSelectdir <- file.path(PROJHOME, 'Plots/Other')


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


##Calculation
detecFracTable <- dataTable %>% group_by(cellID) %>% summarise(DetecFraction = (sum(labels == "B6") + sum(labels == "C7"))/n(), 
                                                               TotalH19 = n(), TotalB6 = sum(labels == "B6"), TotalC7 = sum(labels == "C7")) %>%
  filter(TotalB6 + TotalC7 >= 15)


detecFracTable <- detecFracTable %>% mutate(TotalUndetec = TotalH19 - TotalB6 - TotalC7)
detecFracTable <- detecFracTable %>% arrange(desc(TotalH19))


DetecTotalLine <- ggplot(data = detecFracTable) +
  geom_point(aes(x = TotalB6 + TotalC7, y = TotalH19, size = 3)) +
  geom_smooth(aes(x = TotalB6 + TotalC7, y = TotalH19), method=lm, se=F, formula =y~x, color="red", linetype="dashed") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Total H19 Transcript Counts") +
  xlab("Detected H19 Transcript Counts")
plot(DetecTotalLine)

setwd(outSelectdir)
ggsave(plot = DetecTotalLine,  width = 9.89, height = 9,   filename =  paste("ForReviewerDetecShift.pdf", sep = "_"))



  
  
