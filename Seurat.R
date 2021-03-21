# Droplet-based retina dataset from Macosko et al. (2015),

# BiocManager::install("scRNAseq")
library(scRNAseq)
library(dplyr)
# 1) source("https://z.umn.edu/archived-seurat")
# 2) install.packages('Seurat')
library(Seurat)
library(patchwork)
library(Matrix)

matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_cell_type_ec.csv"))
rownames(matrix_counts) <- matrix_counts[,1]
matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
matrix_counts <- data.matrix(matrix_counts) # must be a matrix object!

data <- CreateSeuratObject(counts = columns_genes_transp_plus_one_log)
data <- CreateSeuratObject(counts = columns_genes_transp_plus_one_log_remove_std_less_than_1)
# data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
# VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 25)

all.genes <- rownames(data)
data <- ScaleData(data, do.scale =  F, do.center = F, features = all.genes)

data <- FindVariableFeatures(object = data, selection.method = 'mvp') #mvp becaue of error in log


data <- RunPCA(data, npcs = 50, features = VariableFeatures(object = data))
# print(data[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(data, reduction = "pca")

data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.5,  algorithm= 1)

data <- RunUMAP(data, dims = 1:50, n.components = 3L)
DimPlot(data, reduction = "umap", split.by = "seurat_clusters")
DimPlot(data, reduction = "umap")

data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
# data_umap_coord <- as.data.frame(data[["pca"]]@cell.embeddings[,3])

umap_dist <- data_umap_coord
# umap_dist <- as.data.frame(data[["umap"]]@cell.embeddings)


