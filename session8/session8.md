# Session 8: R and RStudio

# Expected learning outcome

To understand the basics of scRNA data analysis with R, how to use scRNA packages to creating UMI distribution plots, filtering, normalization, clustering and annotation. 
This is a beginner level lecture in scRNA data analysis

# Overview

- [Introduction to scRNA-Seq](#introduction-to-scrna-seq)
- [Getting Started](#getting_started)
- [scRNA-Seq Processing](#scrna-seq_processing)
- [Introducton to scRNA-Seq Analysis](#introducton_to_scrna-seq_analysis)
- [R Tutorials](#r-tutorials)
- [RNA-Seq data analysis with R](#rna-seq_data_analysis_with_r)

## Introduction to scRNA-Seq

In previous sessions, we have covered basics of R programming, drawing graphics and figures, and introductory RNA-Seq data analysis.
Today, we will briefly cover single cell RNA sequencing (scRNA-Seq), processing of scRNA-Seq reads and then we will extend on what he learn 
on R programming language in previous sessions and start programming with R for essential scRNA-Seq data analysis.  

## Getting Started

We will now install required R packages for conducting today's scRNA analysis practices. 

Run the commands below. It should install the 'devtools' package, necessary for installing scRNA analysis package of 'SignallingSingleCell' from GitHub

    install.packages(“devtools”)
    library(devtools)

Now that we can use 'devtools' package and install GitHub packages, lets install 'SignallingSingleCell'

    devtools::install_github(“garber-lab/SignallingSingleCell”)
    library(SignallingSingleCell)
    
