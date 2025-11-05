library(Seurat)
library(slingshot)
library(tidyverse)

options(future.globals.maxSize = 128000 * 1024^2)

output_dir <- "~/Slingshot_Tracy/results/"

fobj <- readRDS("~/Slingshot_Tracy/Fabio_finalobj.rds")
sub <-  readRDS("~/Slingshot_Tracy/Fabio_sub.rds")

DefaultAssay(fobj) <- "RNA"
Idents(fobj) <- "tracy_clusters"
fobj <- JoinLayers(fobj)


DefaultAssay(fobj) <- "RNA"
Idents(fobj) <- "tracy_clusters"

fobj <- SCTransform(fobj, verbose = FALSE)  
fobj <- RunPCA(fobj, npcs = 50, verbose = FALSE)
fobj <- RunUMAP(fobj, dims = 1:30, verbose = FALSE)   
fobj <- FindNeighbors(fobj, dims = 1:30)
fobj <- FindClusters(fobj, resolution = 0.4)  

DefaultAssay(sub) <- "RNA"
Idents(sub) <- "tracy_clusters"

sub <- SCTransform(sub, verbose = FALSE)  
sub <- RunPCA(sub, npcs = 50, verbose = FALSE)
sub <- RunUMAP(sub, dims = 1:30, verbose = FALSE)   
sub <- FindNeighbors(sub, dims = 1:30)
sub <- FindClusters(sub, resolution = 0.4)  


saveRDS(fobj, file.path("~/Slingshot_Tracy/fobj_sct.rds"))
saveRDS(sub, file.path("~/Slingshot_Tracy/sub_sct.rds"))


fobj_sce <- as.SingleCellExperiment(fobj)
sub_sce <- as.SingleCellExperiment(sub)

pca_fobj <- Embeddings(fobj, "pca")                   
pca_sub  <- Embeddings(sub,  "pca")

reducedDims(fobj_sce)$PCA <- pca_fobj[colnames(fobj_sce), , drop = FALSE]
reducedDims(sub_sce)$PCA  <- pca_sub [colnames(sub_sce),  , drop = FALSE]

Idents(fobj) <- "tracy_clusters"
Idents(sub) <- "tracy_clusters"

colData(fobj_sce)$cluster <- factor(Idents(fobj)[colnames(fobj_sce)])
colData(sub_sce)$cluster  <- factor(Idents(sub) [colnames(sub_sce)])

stopifnot(
  identical(rownames(reducedDims(fobj_sce)$PCA), colnames(fobj_sce)),
  identical(rownames(reducedDims(sub_sce)$PCA),  colnames(sub_sce))
)

fobj_sce <- slingshot(
  fobj_sce,
  clusterLabels = "tracy_clusters",
  reducedDim    = "PCA",
  start.clus    = "Pax6_Prog",
  end.clus      = c("CA1","CA3","DG"),
  stretch       = 1.0
)

sub_sce <- slingshot(
  sub_sce,
  clusterLabels = "tracy_clusters",
  reducedDim    = "PCA",
  start.clus    = "Pax6_Prog",
  end.clus      = c("CA1","CA3","DG"),
  stretch       = 1.0
)



saveRDS(fobj_sce, file.path("~/Slingshot_Tracy/fobj_sce.rds"))
saveRDS(sub_sce, file.path("~/Slingshot_Tracy/sub_sce.rds"))
