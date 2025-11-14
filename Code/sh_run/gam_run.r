library(Seurat)
library(slingshot)
library(tradeSeq)
library(tidyverse)
library(qs)

output_dir <- "~/Slingshot_Tracy/results/"

set.seed(83)

fobj_sce <- readRDS("~/Slingshot_Tracy/fobj_sce.rds")
sub_sce  <- readRDS("~/Slingshot_Tracy/sub_sce.rds")

fobj_sds <- readRDS("~/Slingshot_Tracy/fobj_sds.rds")
fobj_pt  <- readRDS("~/Slingshot_Tracy/fobj_pt.rds")
sub_sds  <- readRDS("~/Slingshot_Tracy/sub_sds.rds")
sub_pt   <- readRDS("~/Slingshot_Tracy/sub_pt.rds")

fobj_counts <- as.matrix(assays(fobj_sce)$counts)
sub_counts  <- as.matrix(assays(sub_sce)$counts)

fobj_gam <- fitGAM(fobj_sds,
  counts    = fobj_counts,
  pseudotime = fobj_pt,
  nknots     = 6,
  verbose    = TRUE
)

sub_gam <- fitGAM(sub_sds,
  counts    = sub_counts,
  pseudotime = sub_pt,
  nknots     = 6,
  verbose    = TRUE
)

saveRDS(fobj_gam, file.path("~/Slingshot_Tracy/fobj_gam.rds"))
saveRDS(sub_gam,  file.path("~/Slingshot_Tracy/sub_gam.rds"))