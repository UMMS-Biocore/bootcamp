# Session 8: Single Cell RNA Analysis

# Expected learning outcome

To understand the basics of scRNA data analysis with R, how to use scRNA packages to creating UMI distribution plots, filtering, normalization, clustering and annotation. 
This is a beginner level lecture in scRNA data analysis.

# Overview

- [Getting Started](#getting-started)
- [Introducton to scRNA-Seq Analysis](#introducton-to-scrna-seq-analysis)
- [scRNA-Seq Processing with DolphinNext](#scrna-seq-processing-with-dolphinnext)
- [Data Structures](#data-structures)
- [Filtering](#filtering)
- [Feature Selection and Dimensionality Reduction](#feature-selection-and-dimensionality-reduction)
- [Exercise 1: Dimension reduction on Skin Data](#exercise-1-dimension-reduction-on-skin-data)
- [Clustering](#clustering)
- [Exercise 2: Clustering on the Skin Data](#exercise-2-clustering-on-the-skin-data)
- [Cell Type Identification](#cell-type-identification)
- [Exercise 3: Cell Type identification Skin Data](#exercise-3-cell-type-identification-skin-data)
- [Normalization](#normalization)
- [Exercise 4: Process the normalized mouse data](#exercise-4-process-the-normalized-mouse-data)
- [Supervised Analysis](#supervised-analysis)
- [DE Analysis](#de-analysis)
- [Homework](#homework)

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
devtools::install_github("kgellatl/SignallingSingleCell")
library(SignallingSingleCell)
```
    
## Introducton to scRNA-Seq Analysis

In previous sessions, we have covered basics of R programming, drawing graphics and figures, and introductory RNA-Seq data analysis.
Today, we will briefly cover single cell RNA sequencing (scRNA-Seq), processing of scRNA-Seq reads and then we will extend on what he have learned 
on R and start programming for the essential scRNA-Seq data analysis.  

Why Single Cell ?

<img src="images/scRNASummary.png" width="600">

We should know that scRNA Analysis may not be necessary for all problems. It is well suited when populations are heterogeneous, and it is a powerful tool for studying intra and inter cell type variations in gene expression

<img src="images/scRNAWorkflow.png" width="600">

<img src="images/barcodes.png" width="600">

<img src="images/barcodes2.png" width="600">

<img src="images/transcript1.png" width="600">

<img src="images/transcript2.png" width="600">

<img src="images/transcript3.png" width="600">

## scRNA-Seq Processing with DolphinNext

To process the fastq files that the instrument generates into a gene expression matrix involves many steps, which can be run as a continuous pipeline with DolphinNext

<img src="images/dnext_singlecell1.png" width="600">

First, DolphinNext pipeline align these reads to reference genome using STAR, HISAT2 or Tophat2 aligners.
Also, mapped reads are merged by samtools, these bam files are also sorted and indexed. 

<img src="images/dnext_singlecell2.png" width="600">

Then, reads are filtered, some cells with low number of reads are removed, and then ESAT (http://garberlab.umassmed.edu/software/esat/) is used to 
create UMI distribution tables.

<img src="images/dnext_singlecell3.png" width="600">

<img src="images/dnext_singlecell4.png" width="600">

## Data Structures

R packages are used to conduct data analysis with built-in (already existing) functions. The idea is to immediately use an algorithm or method for any particular type
of analysis, include scRNA Data Analysis.

We have installed and loaded a library called SignallingSingleCell that include a considerable number of functions and algorithms to conducting scRNA data analysis, and
we will use these functions to cover entire cycle of analyzing single cell data. 

Now, lets download our datasets, load it into our R environment and start investigating! 

```
load(url("https://galaxyweb.umassmed.edu/pub/class/mDC_UMI_Table.Rdata"))
load(url("https://galaxyweb.umassmed.edu/pub/class/ex_sc_skin.Rdata"))
```

Here, "Rdata" is a file format designed for R, and its primary use is to store R objects. We only have a single R object within this Rdata file, and it is called "mDC_UMI_Table". We can investigate this R object further and understand its structure. 

We will use the "ex_sc_skin" dataset during exercises. It is already defined as an ExpressionSet Object. 

```
class(mDC_UMI_Table)

## [1] "matrix" "array" 

dim(mDC_UMI_Table)
 
## [1] 11584  3814
  
head(rownames(mDC_UMI_Table))
   
## [1] "0610007P14Rik" "0610009B22Rik" "0610009O20Rik" "0610010B08Rik" "0610010F05Rik"
## [6] "0610010K14Rik"
   
head(colnames(mDC_UMI_Table))
  
## [1] "0hrA_TGACGGACAAGTAATC" "0hrA_CACAACAGTAGCCTCG" "0hrA_GTTTGTTTGCACCTCT"
## [4] "0hrA_GCTTACCTTGACCCTC" "0hrA_GGAGAAGCGCTTTGGC" "0hrA_AAATCAGAGATCTCGG"

mDC_0hr_1hr_4hr_CLEAN[1:5,1:5] 
   
##               0hrA_TGACGGACAAGTAATC 0hrA_CACAACAGTAGCCTCG 0hrA_GTTTGTTTGCACCTCT
## 0610007P14Rik                     0                     0                     0
## 0610009B22Rik                     0                     0                     0
## 0610009O20Rik                     0                     0                     0
## 0610010B08Rik                     0                     0                     0
## 0610010F05Rik                     0                     0                     0
##               0hrA_GCTTACCTTGACCCTC 0hrA_GGAGAAGCGCTTTGGC
## 0610007P14Rik                     0                     0
## 0610009B22Rik                     0                     0
## 0610009O20Rik                     0                     0
## 0610010B08Rik                     0                     0
## 0610010F05Rik                     0                     0
```

The class function indicates the type of the R object which, in this case, a "matrix" that stores the UMI counts of each barcode in the scRNA experiment associated to each gene.
Columns are the barcodes, and the rows are the genes. The expression matrix contains 11584 genes and 3814 barcodes. 

For analyzing this expression matrix, we have to build a new R object with its own class (like "matrix"), and these classes are often designed within R packages. We 
will use the "ExpressionSet" object which is required for analyzing single cell expression matrices by the SignallingSingleCell package. Lets start by calling the 
"construct_ex_sc" function for that purpose!

```
ex_sc <- construct_ex_sc(mDC_UMI_Table) 
rm(mDC_UMI_Table)
ex_sc
   
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 11584 features, 3814 samples 
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
   
##                       Timepoint
## 0hrA_TGACGGACAAGTAATC       0hr
## 0hrA_CACAACAGTAGCCTCG       0hr
## 0hrA_GTTTGTTTGCACCTCT       0hr
## 0hrA_GCTTACCTTGACCCTC       0hr
## 0hrA_GGAGAAGCGCTTTGGC       0hr
## 0hrA_AAATCAGAGATCTCGG       0hr
 
View(pData(ex_sc))
```

We can also take a subset of the ExpressionSet object

```
mDC_0hr <- subset_ex_sc(ex_sc, variable = "Timepoint", select = c("0hr"))
```

## Filtering

The first step is to filter your data to remove low quality cells. Often creating a histogram of the values and assigning cutoffs is simple and effective. Typically we remove all cells lower than 500-1000 UMIs / cell, and we also remove cells with more than 10000 cells. 

Lets count total UMI counts of all barcodes, and visualize UMI density plots! We will also store these information on metadata as we calculate. 

```
ex_sc <- calc_libsize(ex_sc, suffix = "raw") # sums counts for each cell
ex_sc <- pre_filter(ex_sc, minUMI = 1000, maxUMI = 10000, threshold = 1, minCells = 10,  print_progress = TRUE) 
   
## [1] "Filtering Genes"
## [1] "Filtering Low Cells"

ex_sc <- calc_libsize(ex_sc, suffix = "raw_filtered")
plot_density(ex_sc, title = "UMI Density",  val = "UMI_sum_raw_filtered", statistic = "mean") 
```

<img src="images/umi_density_postfilter.png" width="600">

## Feature Selection and Dimensionality Reduction

Before normalization dimensionality reduction is necessary to form preliminary clusters. These clusters are used to normalize internal to a cluster before normalizing across clusters. First we can subset the genes, and then use these feature selected genes for dimension reduction.

We use the method of Principal Component Analysis (PCA) reduce the dimensionality of the dataset.

```
gene_subset <- subset_genes(ex_sc, method = "PCA", threshold = 1, minCells = 20, nComp = 10, cutoff = 0.85) 
ex_sc <- dim_reduce(ex_sc, genelist = gene_subset, pre_reduce = "iPCA", nComp = 10, tSNE_perp = 30, iterations = 500, print_progress=TRUE) 

colnames(pData(ex_sc))

plot_tsne_metadata(ex_sc, color_by = "Timepoint") 
plot_tsne_metadata(ex_sc, color_by = "iPC_Comp1", title = "PC1 cell loadings") 
plot_tsne_metadata(ex_sc, color_by = "iPC_Comp2", title = "PC2 cell loadings") 
plot_tsne_metadata(ex_sc, color_by = "iPC_Comp3", title = "PC3 cell loadings") 
```

## Exercise 1: Dimension reduction on Skin Data

Now let us try dimension for the skin data!! First we can inspect the skin data to get a sense of what we are working with. Once we get a sense for the data we can then use the same functions as above to perform gene selection and dimension reduction.

```
dim(ex_sc_skin)
colnames(pData(ex_sc_skin))
colnames(fData(ex_sc_skin))
plyr::count(pData(ex_sc_skin))
# Use the subset_genes function to find variable genes in the skin data. Be sure to provide the input argument and method argument. Use ?subset_genes() to view help pages for functions. Try different methods and compare the number of genes you get out for each method.
gene_subset <- subset_genes(ex_sc_skin, method = "PCA", threshold = 1, minCells = 30)

# Use the dim_reduce function to create a 2D representation of the skin data. Be sure to provide the input argument and a gene list.
ex_sc_skin <- dim_reduce(ex_sc_skin, genelist = gene_subset, pre_reduce = "iPCA", nComp = 10, iterations = 500)

# Now you can plot metadata from pData() or genes of interest, onto the tSNE mapping.
plot_tsne_metadata(ex_sc_skin, color_by = "Patient", facet_by = "Skin")
```


## Clustering

Now that we have dimension reduced data we can try clustering it! For dimensions, both Comp and 2d are supported. There will determine if the clustering is done on principal components, or on the 2D representation. There are also 2 clustering algorithms available, density and spectral. Typically we recommend spectral clustering on PCA components, or density clustering on the 2d representation. Try both!

```
ex_sc <- cluster_sc(ex_sc, dimension = "Comp", method = "spectral", num_clust = 4) 
ex_sc$cluster_spectral <- ex_sc$Cluster

ex_sc <- cluster_sc(ex_sc, dimension = "2d", method = "density", num_clust = 4) 
ex_sc$cluster_density <- ex_sc$Cluster

plot_tsne_metadata(ex_sc, color_by = "cluster_spectral")
plot_tsne_metadata(ex_sc, color_by = "cluster_density")
``` 

## Exercise 2: Clustering on the Skin Data

Now let us try clustering for the skin data!! Try both density based and spectral clustering!

```
ex_sc_skin <- cluster_sc(ex_sc_skin, dimension = "Comp", method = "spectral", num_clust = 4) 
plot_tsne_metadata(ex_sc_skin, color_by = "Cluster")
```

## Cell Type Identification

There are many possible ways to identify cell types based on their gene expression. The id_markers function will identify genes that are highly expressed in a high proportion of a given cluster, relative to the other clusters.

```
ex_sc <- id_markers(ex_sc, print_progress = TRUE) 

ex_sc <- calc_agg_bulk(ex_sc, aggregate_by = "Cluster")

markers <- return_markers(ex_sc, num_markers = 15) 

plot_heatmap(ex_sc, genes = unique(unlist(markers)), type = "bulk")
```

## Exercise 3: Cell Type identification Skin Data

Now try to identify the cell types in the skin data!

```
# ex_sc_skin <- id_markers()

# markers <- return_markers() 

# ex_sc_skin <- calc_agg_bulk()

# plot_heatmap()
```

## Normalization

Now that the data has preliminary clusters, we can normalize. SCRAN normalization will first normalize internally in clusters, before normalizing across clusters. Once the data is normalized we can run the same steps as above before visualization. The first step is to select the genes to be used for normalization. One method would be to first only use genes expressed in more than n cells, and then remove the most variable genes. This method can be computationally expensive, and is currently commented out. A simpler approach, counts per million, is also provided below.

```
table(pData(ex_sc)$Cluster)
# ex_sc_norm <- norm_sc(ex_sc, pool_sizes = c(20,25,30,35,40))
x <- exprs(ex_sc)
cSum <- apply(x,2,sum)    # recompute for remaining cells
x <- as.matrix(sweep(x,2,cSum,FUN='/'))*1e6    # normalize to UMIs per million
ex_sc_norm <- construct_ex_sc(x)
pData(ex_sc_norm) <- pData(ex_sc)
```

## Exercise 4: Process the normalized mouse data

Now that we have normalized, it is time to reprocess the data as before, this time on the normalized counts! Use the ex_sc_norm normalized counts from above to run through the same processing steps as before.

```
# gene_subset <- subset_genes()

# ex_sc_norm <- dim_reduce()

# ex_sc_norm <- cluster_sc()

# plot_tsne_metadata()

# plot_tsne_gene()
```

## Supervised Analysis

From the above analysis, it is clear that some clusters are formed based on their cell type, while others are based on their experimental condition. In these cases it can be useful to incorporate prior information in order to obtain clusters and annotations that are grounded in biological significance. Below, we can assign "panels" similar to flow cytometry, that will enable cell groupings based on the expression of genes that you believe to signify biologically relevenant cell types.

```
panel1 <- c("S100a9", "Lcn2") # Neutrophil Markers
panel2 <- c("Ccr7", "Fscn1") # DC
panel3 <- c("Csf1r", "Mertk") # Mac
panels <- list(panel1, panel2, panel3)
plot_tsne_gene(ex_sc_norm, gene = unlist(panels), title = "", log_scale = T)
names(panels) <- c("Neutrophil", "Dendritic", "Macrophage")
ex_sc_norm <- flow_filter(ex_sc_norm, panels = panels, title = "Flow Pass Cells")
ex_sc_norm <- flow_svm(ex_sc_norm, pcnames = "Comp")
plot_tsne_metadata(ex_sc_norm, color_by = "cluster_spectral", title = "Spectral Cluster on PCA components")
plot_tsne_metadata(ex_sc_norm, color_by = "SVM_Classify", title = "Spectral Cluster on PCA components")
```

## DE analysis

Now that cells are grouped by their cell type, we can run DE in order to determine which genes are change in association with our experimental conditions. 

For simplicity we can subset to 0hr and 4hr, to  find the genes that change between these conditions.

It should be noted that DE should always be run on raw counts, not on the normalized counts!

```{r, error=FALSE, warning=FALSE, cache=FALSE, include=TRUE}
ex_sc_norm_0_4 <- subset_ex_sc(ex_sc_norm, variable = "Timepoint", select = c("0hr", "4hr"))
findDEgenes(input = ex_sc,
            pd = pData(ex_sc_norm_0_4),
            DEgroup = "Timepoint",
            contrastID = "0hr",
            facet_by = "SVM_Classify",
            outdir = "~/Downloads/")
plot_volcano(de_path = "~/Downloads/", de_file = "Macrophage_0hr_DEresults.tsv", fdr_cut = 0.000001, logfc_cut = 2)
plot_violin(ex_sc_norm_0_4, color_by = "Timepoint", facet_by = "SVM_Classify", gene = "Cxcl9")
plot_violin(ex_sc_norm_0_4, color_by = "Timepoint", facet_by = "SVM_Classify", gene = "Rsad2")
```

## Homework

Now try to run DE between the 0hr and 1hr timepoints on the mouse data. Then make a volcano plot of the genes that are significantly changed (FDR < 0.001, logfc_cut >= 2) within Dendritic cells.

```{r, error=FALSE, warning=FALSE, cache=FALSE, include=TRUE}
# ex_sc_norm_0_1 <- subset_ex_sc()
# findDEgenes() 
# plot_volcano()
```


## Some notable remarks 

There are numerous R packages and software tools for conducting scRNA data analysis. 
Perhaps, one of the popular of these tools is "Seurat" using R for scRNA analysis, marker analysis, and some other advanced methods. 
Here is a useful link for accessing tutorials of Seurat and how to use it!

https://satijalab.org/seurat/articles/get_started.html
