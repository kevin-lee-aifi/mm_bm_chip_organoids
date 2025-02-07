library(ggplot2, quietly = T)
library(tidyverse, quietly = T)
library(patchwork, quietly = T)
library(dplyr, quietly = T)
library(readxl, quietly = T)
library(Seurat, quietly = T)
library(SeuratObject, quietly = T)
# library(SeuratDisk, quietly = T)
library(harmony, quietly = T)

process_h5 <- function(cr_out) {
    files <- list.files(
        path = file.path(cr_out),
        pattern = "sample_filtered_feature_bc_matrix.h5$",
        recursive = TRUE,
        full.names = TRUE
    )

    so_ls <- lapply(seq_along(files), function(i) {
        name <- unlist(strsplit(files[[i]], '/'))[[9]]
        mtx <- Read10X_h5(files[[i]])
        so <- CreateSeuratObject(counts = mtx, assay = "RNA", min.cells = 1)
        so@meta.data[['replicate']] <- name
        so@meta.data[['orig.ident']] <- unlist(strsplit(name, '_'))[1]
        return(so)
    })

    so <- merge(so_ls[[1]], y = so_ls[2:length(so_ls)])
    so <- JoinLayers(so)

    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    so[["percent.ig"]] <- PercentageFeatureSet(so, pattern = "^IG-")
    
    so <- SCTransform(so, assay = "RNA", verbose = FALSE)
    so <- RunPCA(so, verbose = FALSE)
    so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
    so <- FindClusters(so, resolution = c())
    so <- RunUMAP(so, dims = 1:20)

    return(so)
}