# GO_MWU automation

---
title: "Looping GO-MWU Code"
author: "Elliott Schmidt"
date: "05/02/2021"
output: html_document
---
## Packages needed
```{r setup, include=T, warning=F, message=F}
library(plyr)
library(dplyr)
library(ape)
library(tidyr)
library(tidyverse)
library(readr) 
```

### Quick note before starting

Make sure that all of the GO files needed for the analysis are in the same folder and that there are no other documents in the folder that are the same file type (in this example '.csv'). Also, make sure that nothing is currently in your **global environment**. 

### Lets get started. 

```{r setwd, include=F}
setwd("PATH") 
```

First you should set your working directory, the working directory can be set by using the **setwd()** function in R and placing your PATH in quotations within the parenthesis. For example setwd("C:/Users/YourName/Folder/"). 

Keep in mind that it will be easier for you if your GO files are located along the PATH (ex. PATH = "C:/Users/YourName/Folder/SubFolder/") with the go files located in the SubFolder

#### Lets import our data!

```{r importing data, include=T, warning=F, message=F}
temp = list.files(path = "./FOLDER-NAME/", pattern = "*.csv", full.names = F)
```

The line of code above will make a list containing the names of the files in the folder that your GO files are stored. 

*   In the **list.files(path="./go_files/")** section the "./" means follow the current path, which was set earlier when the directory was set, and the go into the **go_files/** folder, although the folder can be called anything you want. 
*   The **pattern = ".csv"** means get all files that end in .csv  
*   The **full.names** is set to **FALSE** because I just want the files names and not the the files names with the entire PATH attached to the file name.  

```{r importing data pt2, include=T, warning=F, message=F}
mycol <- c("X", "log2FoldChange")    
dflist <- lapply(temp, function(f){ 
  read.csv(paste0("./FOLDER-NAME/", f),row.names = NULL, colClasses = c("character",rep("numeric",6)))[mycol]
}) 
```

The list of characters stored in **mycol** refers to the columns that I want the imported files to contain. I only care about these two columns and not the others. If you care about all the columns then feel free to ignore the command - and makes sure to remove it from the end of the lapply command on the line below! 

The **lapply** function imports all the files in the go_files folder into a list called **dflist**. I don't want R to add a column containing **row.names** so I set it to equal **F**. **colClasses** tells R that the first column should be treated as a **character** string, while the next 6 columns should be treated as **numeric**. 

Next up is a short **for** loop that will go over all the files in the **dflist** and remove rows where there are **NA's** in the **log2FoldChange** column. After the loop the first couple of rows are called to make sure that the loop worked, although if there were initially no **NA's** in the first couple rows you may need to fully call or view a dataframe.

```{r importing data pt3, include=T, warning=F, message=F}
for (i in 1:length(dflist)) {
  dflist[[i]] <- dflist[[i]] %>% drop_na(log2FoldChange)
}; head(dflist[[1]]) 
``` 

>Depending on your preference you may choose to use lapply or for loops to carry out the same task. There is no reasoning as to why I chose one over the other. Which ever one is shown is the one I was able to get working first. Some would consider this abandoment of rationale insane...c'est la vie.  

Next we are going to take all of our dataframes out of a list and place them in our **Global environment** in R. 

```{r importing data pt4, include=T, warning=F, message=F}
for (i in 1:length(dflist)) {
  assign(paste0(temp[i]), as.data.frame(dflist[[i]]))}
``` 

Next we are going to make a list that contains all the **names of the dataframes** in the **global environment**. 

```{r importing data pt5, include=T, warning=F, message=F}
dflist2 <- Filter(function(x) is.data.frame(get(x)), ls()); dflist2
```

**dflist** should still be in the global environment but because it is **list** and not a **dataframe** it will be ignored. You can call your list after the loop has run to see if it worked. Notice how it is just the **names** of your dataframes in the list and not the actually dataframes themselves - although they should still be stored in dflist. 

Next we are going to export our dataframes as **.csv** files. 

>IMPORTANT NOTE: Something that once caused me many hours of pondering and ....frustration, is that for some reason the GO_MWU loop calls on files from the working directory and not from R. I am not sure why this is the case. But just an important FYI. This is why we are exporting our dataframes now that they have been cleaned up a bit. 

```{r exporting data, include=T, warning=F, message=F}
for (h in 1:length(dflist)) {
  write.csv(dflist[h], file = paste0(dflist2[h]), row.names=F, quote=F)
} 
```

Setting **row.names** and **quote** to **F** is important for making sure that your **.csv** come out properly formatted. 

Below is an example of how the GO analysis would be run via automation. However, the original code below for the GO enrichment analysis which is everything after the **for** loop statements **does not belong to me**. **All the credit for the GO enrichment code belongs to the authors on the manuscript described below and can be found by clicking** [here].

[here]: http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1540-2
 

>Wright, R. M., Aglyamova, G. V., Meyer, E. and Matz, M. V. Gene expression associated with white syndromes in a reef-building coral, Acropora hyacinthus. BMC Genomics 2015, 16: 371. 

Also, before running the GO enrichment analysis you will need to have the necessary files: 

*   go.obo file 
*   go annotation file 
*   gomwu.functions.R 
*   gomwu_a.pl 
*   gomwu_b.pl

The [gomwu.functions.R], [go.obo], [gomwu_a.pl], and [gomwu_b.pl] files can be found by clicking the links 

[gomwu.functions.R]:    https://github.com/z0on/GO_MWU 
[go.obo]:   http://geneontology.org/docs/download-ontology/
[gomwu_a.pl]:   https://github.com/z0on/GO_MWU 
[gomwu_b.pl]:   https://github.com/z0on/GO_MWU 



```{r GO enrichment analysis, include=T, warning=F, message=F, eval=F}
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
```

**Congratulations!!** Your Finished and hopefully results are starting to show up in your working directory. Time to take a break and relax as your computer does your work for you!! 
