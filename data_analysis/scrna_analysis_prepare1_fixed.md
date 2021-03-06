### Prepare your environment for Monocle

Create an RStudio project

In the R console run the following commands to install the needed packages to run Monocle3
```r
if (!any(rownames(installed.packages()) == "remotes")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("remotes")
}
library(remotes)

if (!any(rownames(installed.packages()) == "Seurat")){
  BiocManager::install("Seurat")
}
library(Seurat)

if (!any(rownames(installed.packages()) == "R.utils")){
  BiocManager::install("R.utils")
}
library(R.utils)

remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

if (!any(rownames(installed.packages()) == "patchwork")){
  BiocManager::install("patchwork")
}
library(patchwork)

## Monocle3 dependancies
BiocManager::install(c("Biobase", "SingleCellExperiment", "batchelor", "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SummarizedExperiment", "pcaMethods"))

remotes::install_github('cole-trapnell-lab/leidenbase')
remotes::install_github('cole-trapnell-lab/monocle3')

## Test out the installation
library(monocle3)

## Mac users may also experience installation problems due to Xcode or gfortran.
## Instructions: https://cole-trapnell-lab.github.io/monocle3/docs/installation/

# If you had issues with leidenbase and gfortran, trying downloading and installing a newer gfortran
# https://gcc.gnu.org/wiki/GFortranBinaries

## Possible you may have to install Rtools if the above fails
## http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/#1

sessionInfo()
```

### Download the template RMarkdown workshop document for Monocle and open it.

In the R console run the following commands
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/master/data_analysis/monocle.Rmd", "monocle.Rmd")
```

<!-- Additionally, download the following data files for monocle analysis.

```r
download.file("https://github.com/ucdavis-bioinformatics-training/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/raw/master/datasets/monocle3_expression_matrix.rds", "monocle3_expression_matrix.rds")
download.file("https://github.com/ucdavis-bioinformatics-training/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/raw/master/datasets/monocle3_cell_metadata.rds", "monocle3_cell_metadata.rds")
download.file("https://github.com/ucdavis-bioinformatics-training/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/raw/master/datasets/monocle3_gene_metadata.rds", "monocle3_gene_metadata.rds")
``` -->
