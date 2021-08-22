## Project set-up

library(Seurat)
library(ggplot2)
library(velocyto.R)
library(SeuratWrappers)

data_location <- "data_download/"
samples <- c("sample1", "sample2", "sample3")

## colorblind-friendly palettes
tol_high_contrast_palette <- c("#DDAA33", "#BB5566", "#004488")
tol_vibrant_palette <- c("#0077BB", "#33BBEE", "#009988",
                         "#EE7733", "#CC3311", "#EE3377",
                         "#BBBBBB")
tol_muted_palette <- c("#332288", "#88CCEE", "#44AA99",
                       "#117733", "#999933", "#DDCC77",
                       "#CC6677", "#882255", "#AA4499")

# Seurat

raw10x <- lapply(samples, function(i){
  d10x <- Read10X_h5(file.path(data_location, paste0(i, "_raw_feature_bc_matrix.h5")))
  colnames(d10x$`Gene Expression`) <- paste(sapply(strsplit(colnames(d10x$`Gene Expression`),split="-"),'[[',1L),i,sep="-")
  colnames(d10x$`Antibody Capture`) <- paste(sapply(strsplit(colnames(d10x$`Antibody Capture`),split="-"),'[[',1L),i,sep="-")
  d10x
})
names(raw10x) <- samples
gex <- CreateSeuratObject(do.call("cbind", lapply(raw10x,"[[", "Gene Expression")),
                          project = "cellranger multi",
                          min.cells = 0,
                          min.features = 300,
                          names.field = 2,
                          names.delim = "\\-")
gex$percent_mito <- PercentageFeatureSet(gex, pattern = "^MT-")

FeatureScatter(gex,
               feature1 = "nCount_RNA",
               feature2 = "percent_mito",
               shuffle = TRUE) +
  geom_vline(xintercept = 10000) +
  geom_hline(yintercept = 10) +
  scale_color_manual(values = tol_high_contrast_palette)
ggsave("FeatureScatter_nCount_percentMito.png")

FeatureScatter(gex,
               feature1 = "nFeature_RNA",
               feature2 = "percent_mito",
               shuffle = TRUE) +
  geom_vline(xintercept = 2500) +
  geom_hline(yintercept = 10) +
  scale_color_manual(values = tol_high_contrast_palette)
ggsave("FeatureScatter_nFeature_percentMito.png")

FeatureScatter(gex,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               shuffle = TRUE) +
  geom_vline(xintercept = 10000) +
  geom_hline(yintercept = 2500) +
  scale_color_manual(values = tol_high_contrast_palette)
ggsave("FeatureScatter_nCount_nFeature.png")


gex <- subset(gex, percent_mito <= 10)
gex <- subset(gex, nCount_RNA <= 10000)
gex <- subset(gex, nFeature_RNA >= 500)

gex <- NormalizeData(object = gex,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)

gex <- CellCycleScoring(object = gex,
                        s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes,
                        set.ident = FALSE)

gex <- ScaleData(object = gex,
                 vars.to.regress = c("S.Score",
                                     "G2M.Score",
                                     "percent_mito",
                                     "nFeature_RNA"))

gex <- RunPCA(object = gex,
              npcs = 100,
              features = rownames(gex),
              verbose = FALSE)
DimPlot(object = gex,
        reduction = "pca",
        group.by = "orig.ident",
        shuffle = TRUE) +
  scale_color_manual(values = tol_high_contrast_palette)
ggsave("DimPlot_PCA.png")


gex <- FindNeighbors(object = gex,
                     reduction = "pca",
                     dims = c(1:50),
                     verbose = FALSE)
gex <- FindClusters(object = gex,
                    resolution = seq(0.25, 2, 0.25),
                    verbose = FALSE)
sapply(grep("res", colnames(gex@meta.data), value = TRUE),
       function(x) {length(unique(gex@meta.data[,x]))})

gex <- RunUMAP(object = gex,
               reduction = "pca",
               dims = c(1:50),
               verbose = FALSE)

DimPlot(object = gex,
        reduction = "umap",
        group.by = "RNA_snn_res.1.5",
        shuffle = TRUE) +
  scale_color_manual(values = tol_muted_palette)
ggsave("DimPlot_UMAP.png")

# Velocyto

## Read in loom files

reformatLoomColnames <- function(loomObject, assay, id){
  paste(substr(colnames(loomObject[[assay]]), 25, 40), id, sep = "-")
}

loom_data <- lapply(samples, function(sample_id){
  loom_object <- ReadVelocity(file.path(data_location,
                                        paste0(sample_id, ".loom")))
  colnames(loom_object$spliced) <- reformatLoomColnames(
    loomObject = loom_object,
    assay = "spliced",
    id = sample_id)
  colnames(loom_object$unspliced) <- reformatLoomColnames(
    loomObject = loom_object,
    assay = "unspliced",
    id = sample_id)
  colnames(loom_object$ambiguous) <- reformatLoomColnames(
    loomObject = loom_object,
    assay = "ambiguous",
    id = sample_id)
  loom_object
})
names(loom_data) <- samples

loom_aggregate <- lapply(c("spliced", "unspliced"), function(assay){
  do.call("cbind", lapply(loom_data,"[[", assay))
})
names(loom_aggregate) <- c("spliced", "unspliced")
vel <- CreateSeuratObject(counts = loom_aggregate$spliced,
                          project = "Advanced Topics Workshop",
                          assay = "spliced",
                          min.cells = 0,
                          min.features = 0,
                          names.field = 2,
                          names.delim = "\\-")
vel[["unspliced"]] <- CreateAssayObject(counts = loom_aggregate$unspliced)

## Select cells used in gene expression analysis

vel <- vel[rownames(gex), colnames(gex)]

## Run Seurat

vel <- NormalizeData(vel, verbose = FALSE)
vel <- ScaleData(vel, verbose = FALSE)
vel <- RunPCA(vel, features = rownames(vel), verbose = FALSE)

vel_reduction <- FindNeighbors(vel, dims = 1:50, verbose = FALSE)
vel_reduction <- FindClusters(vel_reduction, verbose = FALSE)
vel_reduction <- RunUMAP(vel_reduction, dims = 1:50, verbose = FALSE)

vel_integrated <- vel
vel_integrated@graphs[["RNA_nn"]] <- gex@graphs[["RNA_nn"]][colnames(vel),
                                                            colnames(vel)]
vel_integrated <- AddMetaData(object = vel_integrated,
                              metadata = FetchData(gex,
                                                   vars = "RNA_snn_res.1.5",
                                                   cells = colnames(vel)),
                              col.name = "gex.clusters")
vel_integrated@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(FetchData(gex, vars = c("UMAP_1", "UMAP_2"), cells = colnames(vel))), assay = "RNA")

## Run velocity analysis

vel_reduction <- RunVelocity(vel_reduction,
                             deltaT = 1,
                             kCells = 25,
                             fit.quantile = 0.02,
                             verbose = FALSE)

ident_colors <- tol_vibrant_palette
names(ident_colors) <- levels(vel_reduction)
cell_colors <- ident_colors[Idents(vel_reduction)]
names(cell_colors) <- colnames(vel_reduction)
show.velocity.on.embedding.cor(emb = Embeddings(vel_reduction,
                                                reduction = "umap"),
                               vel = Tool(vel_reduction,
                                          slot = "RunVelocity"),
                               n = 200,
                               scale = "sqrt",
                               xlab = colnames(Embeddings(vel_reduction,
                                                          reduction = "umap"))[1],
                               ylab = colnames(Embeddings(vel_reduction,
                                                          reduction = "umap"))[2],
                               cell.colors = ac(x = cell_colors, alpha = 0.5),
                               cex = 0.8,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5,
                               grid.n = 40,
                               arrow.lwd = 1,
                               do.par = FALSE,
                               cell.border.alpha = 0.1)

vel_integrated <- RunVelocity(vel_integrated,
                             deltaT = 1,
                             kCells = 25,
                             fit.quantile = 0.02,
                             verbose = FALSE)
ident_colors <- tol_muted_palette
names(ident_colors) <- levels(vel_integrated$gex.clusters)
cell_colors <- ident_colors[vel_integrated$gex.clusters]
names(cell_colors) <- colnames(vel_integrated)
show.velocity.on.embedding.cor(emb = Embeddings(vel_integrated,
                                                reduction = "umap"),
                               vel = Tool(vel_integrated,
                                          slot = "RunVelocity"),
                               n = 200,
                               scale = "sqrt",
                               xlab = colnames(Embeddings(vel_integrated,
                                                          reduction = "umap"))[1],
                               ylab = colnames(Embeddings(vel_integrated,
                                                          reduction = "umap"))[2],
                               cell.colors = ac(x = cell_colors, alpha = 0.5),
                               cex = 0.8,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5,
                               grid.n = 40,
                               arrow.lwd = 1,
                               do.par = FALSE,
                               cell.border.alpha = 0.1)

sessionInfo()
