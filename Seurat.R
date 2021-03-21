# Droplet-based retina dataset from Macosko et al. (2015),

# BiocManager::install("scRNAseq")
library(scRNAseq)
library(dplyr)
# 1) source("https://z.umn.edu/archived-seurat")
# 2) install.packages('Seurat')
library(Seurat)
library(patchwork)
library(Matrix)

# log2FC matrix comes from log2_FC_of_genes.R script --> log2FC_MatrixCounts
# After identifying SD<1 from Distribution_analysis_of_genes.R --> SDMoreThanOnelog2FC_MatrixCounts

data <- CreateSeuratObject(counts = log2FC_MatrixCounts)
data <- CreateSeuratObject(counts = SDMoreThanOnelog2FC_MatrixCounts)
# data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
# VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 25)

all.genes <- rownames(data)
data <- ScaleData(data, do.scale =  F, do.center = F, features = all.genes)
data <- FindVariableFeatures(object = data, selection.method = 'mvp') #mvp becaue of error in log
data <- RunPCA(data, npcs = 50, features = VariableFeatures(object = data))
# print(data[["pca"]], dims = 1:5, nfeatures = 5)
# DimPlot(data, reduction = "pca")
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.5,  algorithm= 1) # color in Seurat umap output

data <- RunUMAP(data, dims = 1:50, n.components = 3L)
DimPlot(data, reduction = "umap", split.by = "seurat_clusters")
DimPlot(data, reduction = "umap")

data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
# data_umap_coord <- as.data.frame(data[["pca"]]@cell.embeddings[,3])

umap_dist <- data_umap_coord
# umap_dist <- as.data.frame(data[["umap"]]@cell.embeddings)


