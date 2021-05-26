# Session 8: Exercise Solutions

Installing and loading libraries
```
# install devtools package for installing R packages from GitHub
# install.packages("devtools")
library(devtools)

# install SignallingSingleCell from github
# install_github("kgellatl/SignallingSingleCell")
library(SignallingSingleCell)
```

Loading data 

```
load(url("https://galaxyweb.umassmed.edu/pub/class/ex_sc_skin.Rdata"))
```

## Exercise 1

Filtering 

```
ex_sc_skin <- calc_libsize(ex_sc_skin, suffix = "raw") 
plot_density(ex_sc_skin, title = "UMI Density", val = "UMI_sum_raw", statistic = "mean")
ex_sc_skin <- pre_filter(ex_sc_skin, minUMI = 1000, maxUMI = 10000, threshold = 1, minCells = 10,  print_progress = TRUE) 
ex_sc_skin <- calc_libsize(ex_sc_skin, suffix = "raw_filtered")
plot_density(ex_sc_skin, title = "UMI Density",  val = "UMI_sum_raw_filtered", statistic = "mean") 
```

Subsetting Genes

```
gene_subset <- subset_genes(ex_sc_skin, method = "PCA", threshold = 1, minCells = 30, nComp = 10, cutoff = 0.85) 
```

Dimensionaliy reduction

```
ex_sc_skin <- dim_reduce(ex_sc_skin, genelist = gene_subset, pre_reduce = "iPCA", nComp = 10, tSNE_perp = 30, iterations = 500, print_progress=TRUE)  
colnames(pData(ex_sc_skin))
plot_tsne_metadata(ex_sc_skin, color_by = "UMI_sum_raw", title = "Total UMIs per cell") 
```

## Exercise 2

Clustering
```
ex_sc_skin <- cluster_sc(ex_sc_skin, dimension = "Comp", method = "spectral", num_clust = 8) 
plot_tsne_metadata(ex_sc_skin, color_by = "Cluster", title = "Spectral Cluster on iPCA components") 
```

Cell Type Identification

```
ex_sc_skin <- id_markers(ex_sc_skin, print_progress = TRUE) 
markers <- return_markers(ex_sc_skin, num_markers = 10)
ex_sc_skin <- calc_agg_bulk(ex_sc_skin, aggregate_by = "Cluster")
plot_heatmap(ex_sc_skin, genes = unique(unlist(markers)), type = "bulk")
```
