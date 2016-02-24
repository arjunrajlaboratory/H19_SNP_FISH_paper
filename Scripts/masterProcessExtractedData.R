#Runs processing of Extracted Data and fits 3colorSpots

setwd("/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Scripts/")
source("/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Scripts/Functions/processExtractedData_andFit3colorSpots.R")

datadir <- "/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Raw Data/"
dataSelect <- "Data_CxB_MEF_WT_Control"
pixelShift <- FALSE
outdir <- "/Users/paulginart/Dropbox/H19 SNP-FISH Paper/Processed Data/"
minDetecH19 <- 50
ratioforConversion <- 0.3

byCellCounts <- processExtractedData_andFit3colorSpots(datadir, dataSelect, pixelShift, outdir, minDetecH19)

