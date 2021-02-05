## setting the work directory 
## this will need to be changed to where the files are stored on your own computer. 
## "PATH" is currently used as place holder. 

setwd("PATH") 

## make sure that you set the working directory properly 
## the files 'go.obo', 'GO_MWU.R', 'gomwu.functions.R', 'gomwu_a.pl', 'gomwu_b.pl' and 
## 'GO_deltaRanks_correlation.R' need to be in the directory 

list.files()

## load relavent packages 

library(plyr)
library(dplyr)
library(ape)
library(tidyr)
library(tidyverse)
library(readr)
source("gomwu.functions.R")

## imported in as a list of dataframes 

temp = list.files(path = "./FOLDER-NAME/", pattern = "*.csv", full.names = F) 

## creates a list where all the dataframe are stored but only contain the necessary columns 
## that are stored in 'my col'

mycol <- c("X", "log2FoldChange") 
dflist <- lapply(temp, function(f){ 
  read.csv(paste0("./FOLDER-NAME/", f),row.names = NULL, colClasses = c("character",rep("numeric",6)))[mycol]
}) 

### remove NA's from the log2FoldChange column

for (i in 1:length(dflist)) {
  dflist[[i]] <- dflist[[i]] %>% drop_na(log2FoldChange) 
}; head(dflist[[1]])                                     

#puts all of the dataframes in the dflist in the Global Environment 

for (i in 1:length(dflist)) {
  assign(paste0(temp[i]), as.data.frame(dflist[[i]]))} 

#lists the name of all the dataframes that were in the list 

dflist2 <- Filter(function(x) is.data.frame(get(x)), ls()) 

## writes .csv file for each dataframe 
## files will be saved into the set directory unless otherwise specified

for (h in 1:length(dflist)) {
  write.csv(dflist[h], file = paste0(dflist2[h]), row.names=F, quote=F)
} 

## run GO-enrichment analysis in loop 
## remember that the files 'go.obo', 'GO_MWU.R', 'gomwu.functions.R', 'gomwu_a.pl', 'gomwu_b.pl' and 
## 'GO_deltaRanks_correlation.R' need to be in the directory  

my_comparisons <- c("CC","BP","MF") 
auto_results <- for (i in 1:length(dflist2)) {
  for(g in my_comparisons){ 
    input=paste0(dflist2[[i]])
    goAnnotations="apoly_goannot.txt"
    goDatabase="go.obo"
    goDivision=paste0(g)
    # download .R file that contains a number of functions
    #source("gomwu.functions.R")
    #
    # Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
               largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
               smallest=5,   # a GO category should contain at least this many genes to be considered
               clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
               #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
               #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
               #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
    ) 
    # do not continue if the printout shows that no GO terms pass 10% FDR.
  }
}
