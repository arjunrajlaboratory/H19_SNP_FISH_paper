
library(ggplot2)
library(plyr)
library(dplyr)

source('SET_ABSOLUTE_PATH.R')

datadir <- file.path(PROJHOME, 'Raw Data/5-AZA/extractedCSV/')
outSelectdir <- file.path(PROJHOME, 'Plots/Figure 4/')

master <- data.frame()
#colnames(master) <- c("objArrayNum","objNum","isGood", "AllelicExpression", "area", "alexa.RNACounts", "isTreated", "Genotype", "RepID")

filenames <- list.files(path = datadir)

for (expName in filenames){
  split <- strsplit(expName, "_")
  toJoin <- read.csv(paste(datadir, expName, sep = ""))
  toJoin$isTreated <- split[[1]][1] == "Drug"
  toJoin$Genotype <- split[[1]][2]
  toJoin$RepID <- split[[1]][3]
  master <- bind_rows(master, toJoin)
}

master$Genotype <- as.factor(master$Genotype)
master$RepID <- as.factor(master$RepID)

calculateFracMono <- function(labels){
  mono <- sum(labels == "monoallelic") + sum(labels == "IGF2+")
  total <- mono + sum(labels == "biallelic")
  return(mono/total)
}

master <- master %>% group_by(isTreated, Genotype, RepID) %>% mutate(fracMono = calculateFracMono(AllelicExpression))

master_9CG <- filter(master, Genotype == "9CG") %>% mutate(Condition = ifelse(isTreated, "Drug", "Control"))
master_9CG_forPlot <- master_9CG %>% group_by(RepID, Condition, Genotype) %>% summarise(fracMono = unique(fracMono))
 

master_WT <- filter(master, Genotype == "WT") %>% mutate(Condition = ifelse(isTreated, "Drug", "Control"))
master_WT_forPlot <- master_WT %>% group_by(RepID, Condition, Genotype) %>% summarise(fracMono = unique(fracMono))

p5 <- ggplot() +
  geom_point(data= master_WT_forPlot, aes(x= Condition, y = fracMono, size = 3)) +
  geom_line(data= master_WT_forPlot, aes(x= Condition, y = fracMono, group = RepID)) +
  ylim(c(0,1)) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("WT MEF Monoallelc fraction in Response to 5-AZA treatment")
plot(p5)
setwd(outSelectdir)
ggsave(plot = p5,  width = 8, height = 9,   filename =  paste("WT_MEF_5_AZA.pdf", sep = "_"))
setwd(currdir)

p5 <- ggplot() +
  geom_point(data= master_9CG_forPlot, aes(x= Condition, y = fracMono, size = 3)) +
  geom_line(data= master_9CG_forPlot, aes(x= Condition, y = fracMono, group = RepID)) +
  ylim(c(0,1)) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("9CG MEF Monoallelc fraction in Response to 5-AZA treatment")
plot(p5)
setwd(outSelectdir)
ggsave(plot = p5,  width = 8, height = 9,   filename =  paste("MUT_MEF_5_AZA.pdf", sep = "_"))
setwd(currdir)


master_forPlot <- master %>% mutate(Condition = ifelse(isTreated, "Drug", "Control")) %>% group_by(RepID, Condition, Genotype) %>%
  summarise(fracMono = unique(fracMono))

p6 <- ggplot() +
  geom_point(data= master_forPlot, aes(x= Condition, y = fracMono, colour = Genotype, shape = RepID, size = 3)) +
  geom_line(data= master_WT_forPlot, aes(x= Condition, y = fracMono, group = RepID, color = Genotype)) +
  geom_line(data= master_9CG_forPlot, aes(x= Condition, y = fracMono, group = RepID, color = Genotype)) +
  ylim(c(0,1)) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("MEF Monoallelc fraction in Response to 5-AZA treatment")
plot(p6)

setwd(outSelectdir)
ggsave(plot = p6,  width = 8, height = 9,   filename =  paste("MEF_5_AZA.pdf", sep = "_"))
setwd(currdir)


p7 <- ggplot() +
  geom_point(data= master_forPlot, aes(x= Condition, y = fracMono, colour = Genotype, shape = RepID, size = 3)) +
  scale_shape_discrete(solid=F) +
  ylim(c(0,1)) +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  ggtitle("MEF Monoallelc fraction in Response to 5-AZA treatment")
plot(p7)

setwd(outSelectdir)
ggsave(plot = p7,  width = 8, height = 9,   filename =  paste("MEF_5_AZA_noLine.pdf", sep = "_"))
write.csv(master_forPlot, "5-azaD_Drug Experiment.csv")
setwd(currdir)

##STATISTICS

t.test(master_9CG_forPlot$fracMono~master_9CG_forPlot$Condition == 'Control')

master_9CG_forPlot$fracMono
master_9CG_forPlot$Condition

