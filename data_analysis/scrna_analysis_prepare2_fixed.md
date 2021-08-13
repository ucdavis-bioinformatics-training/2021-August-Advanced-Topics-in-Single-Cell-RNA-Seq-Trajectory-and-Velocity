### Prepare your environment for Velocity

Create an RStudio project

In the R console run the following commands to install the needed packages to run Velocity
```r
if (!any(rownames(installed.packages()) == "devtools")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("devtools")
}
library(devtools)

## install SeuratWrappers
devtools::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

devtools::install_github("velocyto-team/velocyto.R")
library(velocyto.R)

sessionInfo()
```

### Download the template Markdown workshop document for Velocity analysis and open it.

In the R console run the following command
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/master/data_analysis/Velocyto.Rmd", "Velocyto.Rmd")
```
