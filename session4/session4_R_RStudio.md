# Session 4: R and RStudio

# Expected learning outcome

To understand the basics of R and RStudio, how to use R in creating plots, doing table manipulation, and some basic computation. This document aims a beginner audience, rather than experts.

# Overview

- [Introduction](#introduction)
- [Getting Started](#getting-started)
- [R Basics](#r-basics)
- [R Tutorials](#r-tutorials)
- [RNA-Seq data analysis with R](#rna-seq_data_analysis_with_r)

## Introduction

Today we will do some basic exercises to start learning R. R is a programming language and software environment for statistical analysis, graphics representation and reporting. You can either download and install R into your computer or use R in the cluster. You can use R with a user interface called RStudio or you might prefer using R from the command line. There are a couple of ways to use R. In this tutorial we will use the UMass Cluster. **However, we also want you to install R and R-Studio into your local computer for the next session.**

## Getting Started

This week you will use R using R-Studio in the cluster. You can access the web interface using the link below from your local campus network or VPN using your cluster username and password.

For UMMS VPN instructions
<https://www.umassmed.edu/globalassets/it/documents/get-connected/working-remotely/umms-vpn-windows.pdf>

For Web interface, please enter cluster username and password:

<https://ood.umassmed.edu>

- From Interactive Apps => Choose RStudio Server and set the settings:

  - Singularity image file: /pi/alper.kucukural-umw/umw_biocore/dnext_data/singularity/rstudio_bootcamp.sif
  - Queue: Interactive	
  - Number of hours: 4
  - Number of processors: 1
  - MB of memory per core: 64000

## R Basics

### Packages

R consists of lots of functions that can be installed in your environment. If a package is not installed in your environment. You can install it using **install.packages("package_name")** function. If it is a Bioconductor package, you can install it using **BiocManager::install("package_name")**. Bioconductor provides tools for the analysis and comprehension of high-throughput genomic data.

### Setting your working directory

getwd command will let you know your existing working directory

    getwd()

**setwd("directory")** command will set the working directory. For now, we will keep your working directory. However, if you need to read/write files from/to a specific directory, this function might help you.

## R Tutorials

After you connect RStudio. Let's use swirl to start learning R basics. First please set the data directory and load the library. 
    
    library("swirl")
    swirl_options(swirl_data_dir = file.path("~/swirl"))


To start the tutorial just run the command below;

    swirl()

Choose a name for yourself and follow the instructions to start learning. Swirl package will print out all the instruction you to follow.

In this section please finish the tutorial below. We have recorded these tutorials, in case you stuck at any point.

#### R Programming:

1: The basics of programming in R

<https://youtu.be/haHoWQhJEDY></br>

2: Workspace and Files

<https://youtu.be/DoxUXdAkRkg></br>

3: Sequences of Numbers

<https://youtu.be/VtF8oK0iFc8></br>

4: Vectors

<https://youtu.be/lTD-DvSYWc8></br>

5: Missing Values  
<https://youtu.be/05HFXKg7iSY></br>

6: Subsetting Vectors

<https://youtu.be/Os3w3N0Y2MA></br>

7: Matrices and Data Frames

<https://youtu.be/gn6WuV1Tyeo></br>

12: Looking at Data

<https://youtu.be/ppDBvq5VY_8></br>

## RNA-Seq data analysis with R

### Before we start

Let's use editor in R-Studio and keep what will do in this session a Rscript file. Go to "File" -> "New File" -> "R Script"

<img src="images/editor.png">

- **Tip: To run any line in this editor in RStudio, please use Ctrl+Enter (for Windows) or Command+Enter (for Mac) to execute the line.**

### Source necessary functions

    source("/pi/alper.kucukural-umw/umw_biocore/class/funcs.R")

This file is also uploaded to the location below.

### 1. Read the file

We already made gene quantifications and merged expected counts into a single table using ViaFoundry with whole reads (not only using reduced reads) to give you an idea about the complete picture of DESeq analysis.

    file <- "/pi/alper.kucukural-umw/umw_biocore/class/data.tsv"

- If you are using your own computer, **please change the directory above** to wherever your file is located. Use the full path.

Read the data with row.names. We have the gene names in the first column in this table. When we set row.names=1 It will remove the first column from the data and put them to the row name section.

    rsem <- read.table(file,sep="\t", header=TRUE, row.names=1)

    head(rsem)

If file is loaded successfully you will see these lines in your "global environment".

<img src="images/rsem.png" width="400">

### 2. Creating the data structure for DESeq Analysis

We need to select the columns that we are going to use in the analysis. Here we are
going to use the experiment and control data that will be 6 columns. Here in this step make sure to include the columns you want to use in your analysis.

    columns <- c("exper_rep1","exper_rep2","exper_rep3",
             "control_rep1","control_rep2","control_rep3")

    data <- data.frame(rsem[, columns])

### 3. Scatter plot with condition averages

Scatter plots are a key approach to assess differences between conditions. We will plot the average expression value (counts or TPMs for example) of two conditions. First calculate the average reads per gene. To find the average we sum experiments per gene and divide it to the number of replicas by using `rowMeans` function. We will calculate the average for control libraries too. After that, we are ready to merge the average values from experiment and control to create a two column data structure using the `cbind` function that merges tables.

    exp <- c("exper_rep1","exper_rep2","exper_rep3")
    cont <- c("control_rep1","control_rep2","control_rep3")
    avgall <- cbind( rowMeans(data[exp]), rowMeans(data[cont]) )

    head(avgall)

We can also change the column names using `colnames` function.

    colnames(avgall) <- c("Treat", "Control")

    head(avgall)

There are several packages in R that can be used to generate scatter plots. The same plot can be generated using the `ggplot` package.

    gdat<-data.frame(avgall)
    ggplot() +
      geom_point(data=gdat, aes_string(x="Treat", y="Control"),
                 colour="black", alpha=6/10, size=3) +
      scale_x_log10() +scale_y_log10()

<img src="images/scatter.png">

To learn more about ggplots you can check the link below.

<http://r-statistics.co/Complete-Ggplot2-Tutorial-Part1-With-R-Code.html>

### 4. Filtering non-expressed genes

In Eukaryotes only a subset of all genes are expressed in a given cell. The expression is, therefore, a bimodal distribution, with non-expressed genes having counts that result from experimental and biological noise. It is important to filter out the genes that are not expressed before doing differential gene expression.

You can decide which cutoff separates expressed vs non-expressed genes by looking at your histogram we created. In our case, a total sum of 10 counts separates well expressed from non-expressed genes.

    sumd <- rowSums(data)
    hist(log10(sumd), breaks=100)
    abline(v=1)

<img src="images/hist.png">

### 5. All2all scatter plots

To make all to all scatter plot, use the function below.

    all2all(data)

<img src="images/all2all.png">

### 6. Differential Expression Analysis

The goal of Differential gene expression analysis is to find genes or transcripts whose difference in expression, when accounting for the variance within the condition, is higher than expected by chance. The first step is to indicate the condition that each column (experiment)
in the table represent. Here we define the correspondence between columns and conditions.
Make sure the order of the columns matches to your table.

    conds <- factor( c( "Treat", "Treat","Treat",
    "Control","Control", "Control") )

    columns
    #[1] "exper_rep1"   "exper_rep2"   "exper_rep3"  "control_rep1" "control_rep2" "control_rep3"

DESeq will compute the probability that a gene is differentially expressed (DE) for ALL genes in the table. It outputs both a nominal and multiple hypotheses corrected p-value (padj).
Because we are testing DE for over 20,000 genes, we do need to correct for multiple hypothesis testing. To find genes that are significantly DE, we select the ones has lower padj values and higher fold changes and visualize them on our scatter plot with a different color. padj values are corrected p-values which are multiplied by the number of comparisons. Here we are going to use 0.01 for padj value and > 1 log2foldchange. (1/2 < foldChange < 2)

    de_res <- runDESeq(data, columns, conds,  padj=0.01, log2FoldChange=1, non_expressed_cutoff=10)

    ## estimating size factors
    ## estimating dispersions
    ## gene-wise dispersion estimates
    ## mean-dispersion relationship
    ## final dispersion estimates
    ## fitting model and testing

### 7. Make ggplot and overlay the significantly DE genes on the scatter plot

    overlaid_data <- overlaySig(gdat, de_res$res_selected)

    ggplot() +
      geom_point(data=overlaid_data, aes_string(x="Treat", y="Control",
           colour="Legend"), alpha=6/10, size=3) +
      scale_colour_manual(values=c("All"="darkgrey","Significant"="red"))+
      scale_x_log10() +scale_y_log10()

<img src="images/overlaid.png">

### 8. MA Plot

The Second way to visualize it, we use MA plots.
For MA Plot there is another builtin function that you can use.

    plotMA(de_res$res_detected,ylim=c(-2,2),main="DESeq2");

<img src="images/MA.png">

### 9. Volcano Plot

The third way of visualizing the data is by making a Volcano Plot. Here on the x-axis you have log2foldChange values and y-axis you have your -log10 padj values. To see how significant genes are distributed. Highlight genes that have an absolute fold change > 2 and a padj < 0.01

    volcanoPlot(de_res,  padj=0.01, log2FoldChange=1)

<img src="images/volcano.png">

### 9. Clustering

The fourth way of visualizing the data that is widely used in this type of analysis is clustering and Heatmaps.

    sel_data<-data[rownames(de_res$res_selected),]

#### a. Normalizing the selected data using TMM

    norm_data<-getNormalizedMatrix(sel_data, method="TMM")

Here we usually use normalized data and add a pseudocount value 0.1.

    ld <- log2(norm_data+0.1)

- Alternatively, you can use the selected data without normalization `ld <- log2(sel_data+0.1)` and see the difference in the heatmap.

#### b. Scaling the value using their mean centers can be good if the data is uniformly distributed in all the samples.

    cldt <- scale(t(ld), center=TRUE, scale=TRUE);

We can define different distance methods to calculate the distance between samples. Here we focus on euclidean distance and correlation

#### c. Euclidean distance

    distance<-dist(cldt, method = "euclidean")

To plot only the cluster you can use the command below

    plot(hclust(distance, method = "complete"),
    main="Euclidean", xlab="")

<img src="images/cluster_eucledian.png">

##### The heatmap

    cld <- t(cldt)
    heatmap.2(cld, Rowv=TRUE,dendrogram="column",
        Colv=TRUE, col=redblue(256),labRow=NA,
        density.info="none",trace="none", cexCol=0.8);

<img src="images/heatmap_eucledian.png">

#### d. Use Correlation between libraries

Here we calculate the correlation between samples.

    dissimilarity <- 1 - cor(cld)

We define it as distance. This will create a square matrix that will include the dissimilarities between samples.

    distance <- as.dist(dissimilarity)

To plot only the cluster you can use the command below

    plot(hclust(distance, method = "complete"),
     main="1-cor", xlab="")

 <img src="images/cluster_cor.png">

##### The heatmap

    heatmap.2(cld, Rowv=TRUE,dendrogram="column",
        Colv=TRUE, col=redblue(256),labRow=NA,
        density.info="none",trace="none", cexCol=0.8,
        hclust=function(x) hclust(x,method="complete"),
        distfun=function(x) as.dist((1-cor(t(x)))/2))

 <img src="images/heatmap_cor.png">
