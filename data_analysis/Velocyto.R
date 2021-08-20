library(Seurat)
library(ggplot2)
library(velocyto.R)
library(SeuratWrappers)

# custom palettes (colorblind-friendly)
tol_high_contrast_palette <- c("#DDAA33", "#BB5566", "#004488")
tol_vibrant_palette <- c("#0077BB", "#33BBEE", "#009988",
                         "#EE7733", "#CC3311", "#EE3377",
                         "#BBBBBB")
tol_muted_palette <- c("#332288", "#88CCEE", "#44AA99",
                       "#117733", "#999933", "#DDCC77",
                       "#CC6677", "#882255", "#AA4499")

# basic Seurat on cellranger multi
data_location <- "data_download/"
samples <- c("sample1", "sample2", "sample3")

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
ggsave("featureScatter_count_mito.png")

FeatureScatter(gex,
               feature1 = "nFeature_RNA",
               feature2 = "percent_mito",
               shuffle = TRUE) +
  geom_vline(xintercept = 2500) +
  geom_hline(yintercept = 10) +
  scale_color_manual(values = tol_high_contrast_palette)
ggsave("featureScatter_feature_mito.png")

FeatureScatter(gex,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA",
               shuffle = TRUE) +
  geom_vline(xintercept = 10000) +
  geom_hline(yintercept = 2500) +
  scale_color_manual(values = tol_high_contrast_palette)
ggsave("featureScatter_count_feature.png")

gex <- subset(gex, percent_mito <= 10)
gex <- subset(gex, nCount_RNA <= 10000)
gex <- subset(gex, nFeature_RNA >= 500)
table(gex$orig.ident)

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
ggsave("DimPlot_pca.png")

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
        ncol = 2,
        shuffle = TRUE) +
  scale_color_manual(values = tol_muted_palette)
ggsave("DimPlot_umap.png")

# read in Velocyto loom files and reformat column names

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
loom_seurat <- CreateSeuratObject(counts = loom_aggregate$spliced,
                          project = "Advanced Topics Workshop",
                          assay = "spliced",
                          min.cells = 0,
                          min.features = 0,
                          names.field = 2,
                          names.delim = "\\-")
loom_seurat[["unspliced"]] <- CreateAssayObject(counts = loom_aggregate$unspliced)


## Select cells used in gene expression analysis

#If you did not run the initial Seurat section, this is where you would change the code to subset your cells in some other way, or simply use all cells.

loom_seurat_filtered <- loom_seurat[rownames(gex), colnames(gex)]

## Run Seurat

loom_seurat_filtered <- NormalizeData(loom_seurat_filtered,
                                      verbose = FALSE)
loom_seurat_filtered <- ScaleData(loom_seurat_filtered,
                                  verbose = FALSE)
loom_seurat_filtered <- RunPCA(loom_seurat_filtered,
                               features = rownames(loom_seurat_filtered),
                               verbose = FALSE)
loom_seurat_filtered <- FindNeighbors(loom_seurat_filtered,
                                      dims = 1:50,
                                      verbose = FALSE)
loom_seurat_filtered <- FindClusters(loom_seurat_filtered,
                                     verbose = FALSE)
loom_seurat_filtered <- RunUMAP(loom_seurat_filtered,
                                dims = 1:50,
                                verbose = FALSE)

## Run velocity analysis

loom_seurat_filtered <- RunVelocity(loom_seurat_filtered,
                           deltaT = 1,
                           kCells = 25,
                           fit.quantile = 0.02,
                           verbose = FALSE)

ident_colors <- tol_vibrant_palette
names(ident_colors) <- levels(loom_seurat_filtered)
cell_colors <- ident_colors[Idents(loom_seurat_filtered)]
names(cell_colors) <- colnames(loom_seurat_filtered)

pdf("show_velocity1.pdf")
show.velocity.on.embedding.cor(emb = Embeddings(loom_seurat_filtered,
                                                reduction = "umap"),
                               vel = Tool(loom_seurat_filtered,
                                          slot = "RunVelocity"),
                               n = 200,
                               scale = "sqrt",
                               cell.colors = ac(x = cell_colors, alpha = 0.5),
                               cex = 0.8,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5,
                               grid.n = 40,
                               arrow.lwd = 1,
                               do.par = FALSE,
                               cell.border.alpha = 0.1)
dev.off()

# visualize clusters from initial Seurat object

loom_seurat_filtered <- AddMetaData(object = loom_seurat_filtered,
                                    metadata = gex$RNA_snn_res.1.5,
                                    col.name = "gex.clusters")
table(loom_seurat_filtered$gex.clusters,
      loom_seurat_filtered$spliced_snn_res.0.8)

ident_colors2 <- tol_muted_palette
names(ident_colors2) <- levels(loom_seurat_filtered$gex.clusters)
cell_colors2 <- ident_colors2[loom_seurat_filtered$gex.clusters]
names(cell_colors2) <- colnames(loom_seurat_filtered)

pdf("show_velocity2.pdf")
show.velocity.on.embedding.cor(emb = Embeddings(loom_seurat_filtered,
                                                reduction = "umap"),
                               vel = Tool(loom_seurat_filtered,
                                          slot = "RunVelocity"),
                               n = 200,
                               scale = "sqrt",
                               cell.colors = ac(x = cell_colors2,
                                                alpha = 0.5),
                               cex = 0.8,
                               arrow.scale = 3,
                               show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5,
                               grid.n = 40,
                               arrow.lwd = 1,
                               do.par = FALSE,
                               cell.border.alpha = 0.1)
dev.off()
