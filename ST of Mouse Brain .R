devtools::install_github('satijalab/seurat-data')
remotes::install_github("satijalab/seurat")

devtools::install_github("sqjin/CellChat")

install.packages("httpgd")
install.packages("hdf5r", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(CellChat)

library(SeuratData)

library(Seurat)

library(SeuratData)

library(ggplot2)

library(patchwork)

library(dplyr)

library(httpgd)
httpgd::hgd()



library(hdf5r) 

set.seed(1)

 
### WEEK 1

# read Data/loading the data 

section1 <- Load10X_Spatial(

  data.dir= '/home/leaj00001/Desktop/project3/scbi_p3/project_3_dataset/Section_1',

  filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5",

  slice = "slice1",

  image = NULL

)



section2 <- Load10X_Spatial(
  data.dir = '/home/leaj00001/Desktop/project3/scbi_p3/project_3_dataset/Section_2',
  filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5",
  slice = "slice2"
)


## Inspecting the Seurat-object

section1@assays
# 32285 features across 3355 samples within 1 assay

section2@assays
# 32285 features across 3289 samples within 1 assay

## sample image
section1@images
section2@images


# Visualization

# Section_1
section1[["percent.mt"]] <- PercentageFeatureSet(section1, pattern = "^mt-")
section2[["percent.mt"]] <- PercentageFeatureSet(section2, pattern = "^mt-")

plot1 <- VlnPlot(section1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

plot1_Spatial <- SpatialFeaturePlot(section1, features = "nCount_Spatial") + theme(legend.position = "right")

plot1 + plot1_Spatial

 

# Section_2

plot2 <- VlnPlot(section2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

plot2_Spatial <- SpatialFeaturePlot(section2, features = "nCount_Spatial") + theme(legend.position = "right")

plot2 + plot2_Spatial

 

### Data preprocessing

S1 <- subset(section1,subset = nCount_Spatial > 200 & nCount_Spatial < 27000)

S1_cut <- VlnPlot(S1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
S1_cut

S1_cut_Spatial <- SpatialFeaturePlot(S1, features = "nCount_Spatial") + theme(legend.position = "right")
S1_cut_Spatial

S2 <- subset(section2,subset = nCount_Spatial > 200 & nCount_Spatial < 32000)

S2_cut <- VlnPlot(S2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
S2_cut

S2_cut_Spatial <- SpatialFeaturePlot(S2, features = "nCount_Spatial") + theme(legend.position = "right")
S2_cut_Spatial

### Apply SCTransform
S1 <- SCTransform(S1, assay = "Spatial", verbose = FALSE)

 
S2 <- SCTransform(S2, assay = "Spatial", verbose = FALSE)

 ### plots to check the work 

plot1X <- VlnPlot(S1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

plot1Y <- VlnPlot(S2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

plot1X

plot1Y

 

# Dimensionality reduction

S1 <- RunPCA(S1, assay = "SCT", verbose = FALSE)

S1 <- FindNeighbors(S1, reduction = "pca", dims = 1:30)


S1 <- RunUMAP(S1, reduction = "pca", dims = 1:30)

p1_s1 <- DimPlot(S1, reduction = "umap", label = TRUE)
p1_s1

S1 <- FindClusters(S1, verbose = FALSE)

p1_s1 <- DimPlot(S1, reduction = "umap", label = TRUE)
p1_s1

p2_s1 <- SpatialDimPlot(S1, label = TRUE, label.size = 3)

p1_s1 + p2_s1

 ####s2 
S2 <- RunPCA(S2, assay = "SCT", verbose = FALSE)

S2 <- FindNeighbors(S2, reduction = "pca", dims = 1:30)


S2 <- RunUMAP(S2, reduction = "pca", dims = 1:30)

p1_s2 <- DimPlot(S2, reduction = "umap", label = TRUE)
p1_s2

S2 <- FindClusters(S2, verbose = FALSE)
p1_s2 <- DimPlot(S2, reduction = "umap", label = TRUE)
p1_s2

p2_s2 <- SpatialDimPlot(S2, label = TRUE, label.size = 3)

p1_s2 + p2_s2

# DEG based on the clusters 

DEG.markers <- FindMarkers(S2, ident.1 = 1, 
              ident.2 = c(0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
              min.pct = 0.25)


DEG_plot <- SpatialFeaturePlot(object = S2, features = rownames(DEG.markers)[1:3], alpha = c(0.1, 1))

DEG_plot

###### DEG analysis based on the spatial patterning

S2_featuress <- FindSpatiallyVariableFeatures(S2, assay = "SCT", features = VariableFeatures(S2)[1:100],
    selection.method = "markvariogram")
### top 3 features 
### https://github.com/satijalab/seurat-object/issues/25 (citation for the functions used)
top.clusters <- rownames(S2_featuress[["SCT"]]@meta.features)[which(S2_featuress[["SCT"]]@meta.features$moransi.spatially.variable.rank < 6)]
top.clusters <- head(VariableFeatures(S2_featuress, method = "markvariogram"), 3)
top.clusters
SpatialFeaturePlot(
  object = S2,
  features = top.clusters,
  alpha = c(0.1, 1)   
)
#### based on P_Value ( we did not know which one is more correct so we did both)

#### FindAllMarkers Function (citation)
markers <- FindAllMarkers(
  S2,
  assay = "SCT",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

### dplyr Functions
top_markers <- markers %>%
  arrange(p_val_adj) %>%
  slice_head(n = 3) %>%
  pull(gene)

##### SpatialFeaturePlot Function
SpatialFeaturePlot(
  object = S2,
  features = top_markers,
  alpha = c(0.1, 1)
)
### merging without batch correction
merged <- merge(S1, S2)
DefaultAssay(merged) <- "Spatial"
merged <- SCTransform(merged, assay = "Spatial", verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, assay = "SCT",  verbose = FALSE)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)
merged<- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)
DimPlot(merged, reduction = 'umap')
SpatialDimPlot(merged, label = TRUE, label.size = 3)


##### integration for batch correction
integration <-  FindIntegrationAnchors(c(S1, S2))
integration <- IntegrateData(anchorset = integration)
DefaultAssay(integration) <- "integrated"
integration <- SCTransform(integration, assay = "Spatial",  verbose = FALSE)
integration <- RunPCA(integration, npcs = 30, assay = "SCT", verbose = FALSE)
integration <- RunUMAP(integration, reduction = "pca", dims = 1:30)
integration <- FindNeighbors(integration, reduction = "pca", dims = 1:30)
integration <- FindClusters(integration, resolution = 0.5)
DimPlot(integration, reduction = 'umap')
SpatialDimPlot(integration, label = TRUE, label.size = 3) 



###### annotation 
reference_read <- readRDS("/home/leaj00001/Desktop/project3/allen_cortex.rds")
###Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : The total size of the 19 globals exported for future expression (‘FUN()’) is 2.86 GiB.. This exceeds the maximum 
#allowed size of 500.00 MiB (option 'future.globals.maxSize'). The three largest globals are ‘FUN’ (2.84 GiB of class ‘function’), ‘umi_bin’ (19.21 MiB of class ‘numeric’) and ‘data_step1’ (3.48 MiB of class ‘list’)
## therefore we will maximize the size :
options(future.globals.maxSize = 5 * 1024^3) 
library(future)
plan("sequential") 
reference <- SCTransform(reference_read, verbose = FALSE)
reference <- RunPCA(reference, verbose = FALSE)
reference <- RunUMAP(reference, dims = 1:30)
reference <- FindNeighbors(reference, dims = 1:30)
reference <- FindClusters(reference, resolution = 0.5)
DimPlot(reference, reduction = "umap", label = TRUE) + ggtitle("Reference Dataset Clustering")


####### Manual Annotation
M_Annotation <- RunPCA(reference, assay = "SCT", verbose = FALSE)

M_Annotation <- FindNeighbors(M_Annotation, reduction = "pca", dims = 1:30)


M_Annotation <- RunUMAP(M_Annotation, reduction = "pca", dims = 1:30)

p1_M_Annotation <- DimPlot(M_Annotation, reduction = "umap", label = TRUE)
p1_M_Annotation

M_Annotation <- FindClusters(M_Annotation, verbose = FALSE)
p2_M_Annotation <- DimPlot(M_Annotation, reduction = "umap", label = TRUE)
p2_M_Annotation


### Plot Marker Gene Expression in UMAP
FeaturePlot(
  object = merged,
  features = c("Abi3bp", "Acp5", "Abcb1a", "Abhd2", "LPPR1", "LRRC17"),   
  reduction = "umap"
)

#### Plot on the tissue slides.
plot <- SpatialFeaturePlot(
  object = merged,
  features = c("Abi3bp", "Acp5", "Abcb1a", "Abhd2", "LPPR1", "LRRC17"),   
  alpha = c(0.1, 1)
)
plot + plot_layout(ncol = 4)
