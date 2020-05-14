# Signle cell RNA clusering analysis.

# Yu jinge 

## Data

### Abstract 
Our data comes from NCBI with accession number GSE111672, and we choose PDAC-B inDrop1 as scRNA data, PDAC-B ST1 as ST data. Both data set comes from single-cell suspension of tissue of pancreatic adenocarcinoma from homo sapiens. The scRNA data consists of 19831 genes and 2000 cells, spot data has 16528 genes and 996 spots. 

### Availability 

The data are publicly available for download via the online data portal at <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3405534> and <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3405531>. No registration is required.

### Description 

Notice that during experinments, a tumor is first divided and a single-cell suspension is generated from one ption and processed for scRNA-seq to identify the cell populations present in the tissue. From the remaining tissue, cryosection are processed transcripts across the tissue. 

Since genes in origianl data sets are large, we only select the most variable 500 genes of both data separately. Besides, I have normalized and continuously processed both data, the processed data are in Data file.

## Code

### Abstract

Most of the data processing and analysis for this report were done in R, and the parameter updating functioins are witten in Rcpp since MCMC method needs a lot of time. The corresponding code is provided to take exploratory data analysis on the raw data; conduct various preprocessing steps; fit a set of Bayesian models to the preprocessed data via Markov chain Monte Carlo (MCMC) methods; and generate descriptive plots used in the
report.

### Description

All of the R scripts used in the report are available in a public repository on GitHub [https://github.com/jingeyu/single-cell-RNA-cluster]. The MIT license applies to all code, and no permissions are required to access the code.

### Optional Information

R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin15.6.0 (64-bit) were used in the project. And The necessary R libraries for the code used for data processing and analysis are:

- dplyr, version 0.8.5 (https://cloud.r-project.org/web/packages/dplyr/index.html)
- plotly, version 4.9.2 (https://cloud.r-project.org/web/packages/plotly/index.html)
- ggplot2, version 3.3.0 (https://cloud.r-project.org/web/packages/ggplot2/index.html)
- factoextra, version 1.0.7 (https://cloud.r-project.org/web/packages/factoextra/index.html)
- cluster, version 2.1.0 (https://cloud.r-project.org/web/packages/cluster/index.html)
- stringr, version 1.4.0 (https://cloud.r-project.org/web/packages/stringr/index.html)
- miscTools, version 0.6-26 (https://cloud.r-project.org/web/packages/miscTools/index.html)
- umap, version 0.2.5.0 (https://cloud.r-project.org/web/packages/umap/index.html)

Computer information:
CPU: 2.4 GHz Intel Core i5, 8 GB 2133 MHz LPDDR3. 


## Instructions for Use

### Reproducibility

All data preparation and analyses are reproduced, as well as all Figures in the
report.

All workflow information is contained in the Reproduce.R script. The general steps
are:

1. Exploratory Data Analysis
2. Data processing 
3. Single cell RNA data cluster
4. Spot data cluster

Then all the pictures in my slides and report are reproduced.
5. Run Reproduce.R. Then all the pictures in my slides and report are reproduced.
