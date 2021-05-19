# Session 8: Single Cell RNA Analysis

# Expected learning outcome

To understand the basics of scRNA data analysis with R, how to use scRNA packages to creating UMI distribution plots, filtering, normalization, clustering and annotation. 
This is a beginner level lecture in scRNA data analysis.

# Overview

- [Introduction to scRNA-Seq](#introduction-to-scrna-seq)
- [Getting Started](#getting-started)
- [scRNA-Seq Processing](#scrna-seq-processing)
- [Introducton to scRNA-Seq Analysis](#introducton-to-scrna-seq_analysis)
- [Data Structures](#data-structures)
- [Normalization](#normalization)
- [Feature Selection](#feature-selection)
- [Dimensionality Reduction](#dimensionality-reduction)
- [Clustering](#clustering)
- [Differential expression](#differential-expression)
- [Supervised analysis](#supervised-analysis)

## Introduction to scRNA-Seq

In previous sessions, we have covered basics of R programming, drawing graphics and figures, and introductory RNA-Seq data analysis.
Today, we will briefly cover single cell RNA sequencing (scRNA-Seq), processing of scRNA-Seq reads and then we will extend on what he learn 
on R programming language in previous sessions and start programming with R for essential scRNA-Seq data analysis.  

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

When it is asked. "Do you want to install from sources the package which needs compilation? (Yes/no/cancel)". If so, just type "n" for no, and hit Enter to continue.

    Do you want to install from sources the package which needs compilation? (Yes/no/cancel) n

If R asks to update old packages; please dont write anything and just enter:

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

Similar to before, answer "n" to if R asks you to install from source packages, and just enter if R asks to update packages

```
    devtools::install_github(“kgellatl/SignallingSingleCell”)
    library(SignallingSingleCell)
```
    
## Introducton to scRNA-Seq Analysis

TBA
