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
- [Normalization](#normalization)
- [Feature Selection](#feature-selection)
- [Dimensionality Reduction](#dimensionality-reduction)
- [Clustering](#clustering)
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

TBA

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
   > class(mDC_0hr_1hr_4hr_CLEAN)
   [1] "matrix" "array" 
   > dim(mDC_0hr_1hr_4hr_CLEAN)
   [1] 14222 24169
   > head(rownames(mDC_0hr_1hr_4hr_CLEAN))
   [1] "0610007P14Rik" "0610009B22Rik" "0610009O20Rik" "0610010B08Rik" "0610010F05Rik"
   [6] "0610010K14Rik"
   > head(colnames(mDC_0hr_1hr_4hr_CLEAN))
   [1] "0hrA_CATTTGTTCTAGACCC" "0hrA_CCTACTAGATCTTTGT" "0hrA_CAACAAATATATAGGA"
   [4] "0hrA_AAATCAGACACAACAG" "0hrA_GCGTTGCTTCTGTGGT" "0hrA_TGACGGACAAGTAATC"
   > mDC_0hr_1hr_4hr_CLEAN[1:5,1:5] 
                 0hrA_CATTTGTTCTAGACCC 0hrA_CCTACTAGATCTTTGT 0hrA_CAACAAATATATAGGA
   0610007P14Rik                     0                     0                     0
   0610009B22Rik                     0                     0                     0
   0610009O20Rik                     0                     0                     0
   0610010B08Rik                     0                     0                     0
   0610010F05Rik                     0                     0                     0
                 0hrA_AAATCAGACACAACAG 0hrA_GCGTTGCTTCTGTGGT
   0610007P14Rik                     0                     0
   0610009B22Rik                     0                     0
   0610009O20Rik                     0                     0
   0610010B08Rik                     0                     0
   0610010F05Rik                     0                     0
```

The class function indicates the type of the R object which, in this case, a matrix that stores the UMI counts of each barcode in the scRNA experiment associated to each gene.
Columns are the barcodes, and the rows are the genes. The expression matrix contains 14222 genes and 24169 barcodes. 






