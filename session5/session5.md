Session 5: R and RStudio
========

Expected learning outcome
========

To understand the basics of R and RStudio, how to use R in creating plots, doing table manuplilation and some basic computation. This document aims beginner audience, rather than experts.

Overview
========

  * [Introduction](#introduction)
  * [Getting Started](#getting-started)


## Introduction
Today we will do some basic excercise to start learning R. R is a programming language and software environment for statistical analysis, graphics representation and reporting. You can either download and install R into your computer or use R in the cluster. You can use R with a user interface called RStudio or you might prefer using R from command line. There are a couple of ways to use R. In this tutorial we will use UMass Cluster. However, we also want you to install R and R-Studio into your local computer for the next session.

## Getting Started

This week you will use R using R-Studio in the cluster. You can access the web interface using the link below from your local campus network or VPN using your cluster username and password.

For UMMS VPN instructions
<https://www.umassmed.edu/globalassets/it/documents/get-connected/working-remotely/umms-vpn-windows.pdf>

For Web interface

<https://www.umassrc.org:444>

To access R-Studio in the cluster. Please use the short video below.

<https://www.youtube.com/watch?v=S6UO-KEy8ew>

### Installing R and R-Studio to your local computer

Please install R 4.0.0 in to your local computer to use it in the next session. 

1. First please install R to your computer

	<https://cran.r-project.org/>
	
	For windows;
	<https://cran.r-project.org/bin/windows/base/R-4.0.0-win.exe>
	
	For mac;
	<https://cran.r-project.org/bin/macosx/R-4.0.0.pkg>

2. Second, please install R-Studio;

	<https://rstudio.com/products/rstudio/download/>


3. Please install debrowser to your local computer. This can take time to install all the packages.

	Open R or R-Studio and run the commands in R, it will include most of the packages we will use;

```	
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	
BiocManager::install("debrowser") 
```


    