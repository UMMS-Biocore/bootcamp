# Session 8: Single Cell RNA Analysis

# Expected learning outcome

To understand the basics of scRNA data analysis with R, how to use scRNA packages to creating UMI distribution plots, filtering, normalization, clustering and annotation. 
This is a beginner level lecture in scRNA data analysis.

# Overview

- [Introduction to scRNA-Seq](#introduction-to-scrna-seq)
- [Getting Started](#getting-started)
- [Introducton to scRNA-Seq Analysis](#introducton-to-scrna-seq-analysis)
- [scRNA-Seq Processing](#scrna-seq-processing)
- [Data Structures](#data-structures)
- [Filtering](#filtering)
- [Feature Selection and Dimensionality Reduction](#feature-selection-and-dimensionality-reduction)
- [Clustering](#clustering)
- [Normalization](#normalization)
- [Differential expression](#differential-expression)

## Introduction to scRNA-Seq

In previous sessions, we have covered basics of R programming, drawing graphics and figures, and introductory RNA-Seq data analysis.
Today, we will briefly cover single cell RNA sequencing (scRNA-Seq), processing of scRNA-Seq reads and then we will extend on what he have learned 
on R and start programming for the essential scRNA-Seq data analysis.  

Why Single Cell ?

<img src="images/scRNASummary.png" width="600">

## Getting Started

Please install R 4.0.5 into your local computer. 

1. First please install R to your computer

   <https://cran.r-project.org/>

   For windows;
   <https://cran.r-project.org/bin/windows/base/R-4.0.5-win.exe>

   For mac;
   <https://cran.r-project.org/bin/macosx/R-4.0.5.pkg></br>

2. Second, please install R-Studio;

   <https://rstudio.com/products/rstudio/download/>
   
3. We will now install required R packages for conducting today's scRNA analysis practices before covering the basics of scRNA-Seq. Run the commands below. It should install the 'devtools' package, necessary for installing scRNA analysis package of 'SignallingSingleCell' from GitHub. **Installing devtools may take long, please try to install before the session**.

```
    install.packages(“devtools”)
    library(devtools)
```

If Rstudio asks "Do you want to install from sources the package which needs compilation? (Yes/no/cancel)", just type "n" for no, and hit Enter to continue.

    Do you want to install from sources the package which needs compilation? (Yes/no/cancel) n

If R asks to update old packages; please dont write anything and hit Enter to continue:

```
   These packages have more recent versions available.
   It is recommended to update all of them.
   Which would you like to update?
   
   1: All                            
   2: CRAN packages only             
   3: None                           
   4: tibble  (3.1.1 -> 3.1.2) [CRAN]
   5: stringi (1.5.3 -> 1.6.2) [CRAN]
   
   Enter one or more numbers, or an empty line to skip updates:
```

4. Now that we can use 'devtools' package and install GitHub packages, lets install 'SignallingSingleCell'.

Similar to Step 3, answer "n" or "no" if Rstudio asks you to install from source packages, and just enter if R asks to update packages

```
    devtools::install_github(“kgellatl/SignallingSingleCell”)
    library(SignallingSingleCell)
```
    
## Introducton to scRNA-Seq Analysis

<img src="images/scRNAWorkflow.png" width="600">

## Data Structures

R packages are used to conduct data analysis with built-in (already existing) functions. The idea is to immediately use an algorithm or method for any particular type
of analysis, include scRNA Data Analysis.

We have installed and loaded a library called SignallingSingleCell that include a considerable number of functions and algorithms to conducting scRNA data analysis, and
we will use these functions to cover entire cycle of analyzing single cell data. 

Now, lets download our data, load it into our R environment and start investigating! 

```
   load("~/Documents/UMASS/Garber/Data/inDrop/mDC_UMIClean/1-RawFiles/mDC_0hr_1hr_4hr_CLEAN.Rdata")
   load("mDC_0hr_1hr_4hr_CLEAN.Rdata")
```

Here, "Rdata" is a file format designed for R, and its primary use is to store R objects. We only have a single R object within this Rdata file, and it is called "mDC_0hr_1hr_4hr_CLEAN". We can investiate this R object further and understand its structure. 

```
   class(mDC_0hr_1hr_4hr_CLEAN)
   
   ## [1] "matrix" "array" 
  
   dim(mDC_0hr_1hr_4hr_CLEAN)
   
   ## [1] 14222 24169
   
   head(rownames(mDC_0hr_1hr_4hr_CLEAN))
   
   ## [1] "0610007P14Rik" "0610009B22Rik" "0610009O20Rik" "0610010B08Rik" "0610010F05Rik"
   ## [6] "0610010K14Rik"
   
   head(colnames(mDC_0hr_1hr_4hr_CLEAN))
   
   ## [1] "0hrA_CATTTGTTCTAGACCC" "0hrA_CCTACTAGATCTTTGT" "0hrA_CAACAAATATATAGGA"
   ## [4] "0hrA_AAATCAGACACAACAG" "0hrA_GCGTTGCTTCTGTGGT" "0hrA_TGACGGACAAGTAATC"
   
   mDC_0hr_1hr_4hr_CLEAN[1:5,1:5] 
   
   ##               0hrA_CATTTGTTCTAGACCC 0hrA_CCTACTAGATCTTTGT 0hrA_CAACAAATATATAGGA
   ## 0610007P14Rik                     0                     0                     0
   ## 0610009B22Rik                     0                     0                     0
   ## 0610009O20Rik                     0                     0                     0
   ## 0610010B08Rik                     0                     0                     0
   ## 0610010F05Rik                     0                     0                     0
   ##               0hrA_AAATCAGACACAACAG 0hrA_GCGTTGCTTCTGTGGT
   ## 0610007P14Rik                     0                     0
   ## 0610009B22Rik                     0                     0
   ## 0610009O20Rik                     0                     0
   ## 0610010B08Rik                     0                     0
   ## 0610010F05Rik                     0                     0
```

The class function indicates the type of the R object which, in this case, a "matrix" that stores the UMI counts of each barcode in the scRNA experiment associated to each gene.
Columns are the barcodes, and the rows are the genes. The expression matrix contains 14222 genes and 24169 barcodes. 

For analyzing this expression matrix, we have to build a new R object with its own class (like "matrix"), and these classes are often designed within R packages. We 
will use the "ExpressionSet" object which is required for analyzing single cell expression matrices by the SignallingSingleCell package. Lets start by calling the 
"construct_ex_sc" function for that purpose!

```
   ex_sc <- construct_ex_sc(mDC_0hr_1hr_4hr_CLEAN) # sc_dat == Input expression matrix
   rm(mDC_0hr_1hr_4hr_CLEAN)
   ex_sc
   
   ## ExpressionSet (storageMode: lockedEnvironment)
   ## assayData: 14222 features, 24169 samples 
   ##   element names: exprs 
   ## protocolData: none
   ## phenoData: none
   ## featureData: none
   ## experimentData: use 'experimentData(object)'
   ## Annotation:  
```

Now the expression matrix is stored within an R object with ExpressionSet class, and we wont need the original expression matrix, we may delete it with "rm" function. 
The ExpressionSet class (ex_sc) is an extremely convienient data structure that contains 3 dataframes. These dataframes contain expression data, cell information, and gene information respectivelty. 

exprs(ex_sc) is the expression data, where rows are genes and columns are cells  
pData(ex_sc) is cell information, where rows are cells and columns are metadata  
fData(ex_sc) is gene information, where rows are genes and columns are metadata

Often we have metadata on the experiment that can be valuable in the analysis! Writing that information now may be appropriate. Our experiment consists of a time 
course with LPS stimulation. Now we can begin to take advantage of our faceting! We first create an empty column, then insert associated metadata. 

```
   pData(ex_sc)$Timepoint <- NA 
   pData(ex_sc)[grep("0hr", rownames(pData(ex_sc))),"Timepoint"] <- "0hr"
   pData(ex_sc)[grep("1hr", rownames(pData(ex_sc))),"Timepoint"] <- "1hr"
   pData(ex_sc)[grep("4hr", rownames(pData(ex_sc))),"Timepoint"] <- "4hr"
   head(pData(ex_sc))
   
   ##                          Timepoint
   ## 0hrA_CATTTGTTCTAGACCC          0hr
   ## 0hrA_CCTACTAGATCTTTGT          0hr
   ## 0hrA_CAACAAATATATAGGA          0hr
   ## 0hrA_AAATCAGACACAACAG          0hr
   ## 0hrA_GCGTTGCTTCTGTGGT          0hr
   ## 0hrA_TGACGGACAAGTAATC          0hr
   
   View(pData(ex_sc))
```

## Filtering

The first step is to filter your data to remove low quality cells. Often creating a histogram of the values and assigning cutoffs is simple and effective. Typically we remove all cells lower than 500-1000 UMIs / cell, and we also remove cells with more than 10000 cells. The objective is to remove fragments from barcodes whose droplets are either empty or includes more than one cell. Low or considerably high UMI counts are indications of such technical errors that occur during scRNA sequencing. 

Lets count total UMI counts of all barcodes, and visualize UMI density plots! We will also store these information on metadata as we calculate. 

```
   ex_sc <- calc_libsize(ex_sc, suffix = "raw") # sums counts for each cell
   plot_density(ex_sc, title = "UMI Density", val = "UMI_sum_raw", statistic = "mean")    
```

<img src="images/umi_density_prefilter.png" width="600">

```
   ex_sc <- pre_filter(ex_sc, minUMI = 1000, maxUMI = 10000, threshold = 1, minCells = 10,  print_progress = TRUE) 
   
   ## [1] "Filtering Genes"
   ## [1] "Filtering Low Cells"

   ex_sc <- calc_libsize(ex_sc, suffix = "raw_filtered")
   plot_density(ex_sc, title = "UMI Density",  val = "UMI_sum_raw_filtered", statistic = "mean") 
   
   head(pData(ex_sc))
   
   ##                          Timepoint UMI_sum_raw UMI_sum_raw_filtered
   ## 0hrA_TGACGGACAAGTAATC          0hr        4579                 4572
   ## 0hrA_ATGGGCACACCTTGCC          0hr        2897                 2896
   ## 0hrA_TCGAAGCTGTTGCACG          0hr        2266                 2265
   ## 0hrA_TGTTTGAGTCGGTTCG          0hr        3045                 3042
   ## 0hrA_TAAATAGGCACAAGGC          0hr        2238                 2237
   ## 0hrA_GATTAGACGGGAACCT          0hr        1217                 1217
```

<img src="images/umi_density_postfilter.png" width="600">

## Feature Selection and Dimensionality Reduction

Before normalization dimensionality reduction is necessary to form preliminary clusters. These clusters are used to normalize internal to a cluster before normalizing across clusters. First we can subset the genes, and then use these feature selected genes for dimension reduction.

```
   gene_subset <- subset_genes(ex_sc, method = "PCA", threshold = 1, minCells = 30, nComp = 10, cutoff = 0.85) 
   ex_sc <- dim_reduce(ex_sc, genelist = gene_subset, pre_reduce = "iPCA", nComp = 10, tSNE_perp = 30, iterations = 500, print_progress=TRUE)  
   colnames(pData(ex_sc))

   ##  [1] "Timepoint"            "UMI_sum_raw"          "UMI_sum_raw_filtered"
   ##  [4] "x"                    "y"                    "iPC_Comp1"           
   ##  [7] "iPC_Comp2"            "iPC_Comp3"            "iPC_Comp4"           
   ## [10] "iPC_Comp5"            "iPC_Comp6"            "iPC_Comp7"           
   ## [13] "iPC_Comp8"            "iPC_Comp9"            "iPC_Comp10" 

   plot_tsne_metadata(ex_sc, color_by = "UMI_sum_raw", title = "Total UMIs per cell") 
```

<img src="images/umi_tsne.png" width="600">

## Clustering

Now that we have dimension reduced data we can try clustering it! Clustering algorithms often require you to specify a parameter. This is either the number of clusters, or
parameter that represents some "resolution" which might represents a threshold for the similarity between expression profiles of barcodes. 

```
   ex_sc <- cluster_sc(ex_sc, dimension = "Comp", method = "spectral", num_clust = 6) 
   plot_tsne_metadata(ex_sc, color_by = "Cluster", title = "Spectral Cluster on iPCA components") 
```

<img src="images/clusters_tsne.png" width="600">

```
   plot_density(ex_sc, title = "UMIs per cluster", val = "UMI_sum_raw", color_by = "Cluster", statistic = "mean")
```

<img src="images/umi_density_clusters.png" width="600">

We can also calculate a set of markers for these clusters and store the scores of all genes to fData. This is a quick method to find good markers genes for cell identification. These gene scores get written to fData()

```
   ex_sc <- id_markers(ex_sc, print_progress = TRUE) 
   markers <- return_markers(ex_sc, num_markers = 10) 
```
