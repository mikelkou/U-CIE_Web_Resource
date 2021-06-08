# Droplet-based retina dataset from Macosko et al. (2015),

# BiocManager::install("scRNAseq")
library(scRNAseq)
library(Matrix)

matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_time_course_ec.csv")) # 1
matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_cell_type_ec.csv")) # 2
matrix_counts <- data.frame(read.delim("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE52529_fpkm_matrix.txt")) # 3
rownames(matrix_counts) <- matrix_counts[,1]
matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
matrix_counts <- data.matrix(matrix_counts) # must be a matrix object!

#--- Osteoclasts ---#
library(Matrix)
matrix_dir = "../osteoclasts/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

matrix_counts <- data.matrix(mat)
# matrix_counts <- matrix_counts[sample(nrow(matrix_counts), 20000), ]
# matrix_counts <- matrix_counts[, sample(ncol(matrix_counts), 2000) ]
#------------------------------------------------------------------------------#

TransposedMatrixCounts <- t(matrix_counts)

cnt <- 0
FC_MatrixCounts <- c()
for(i in 1:ncol(TransposedMatrixCounts)){
  mean_col <- mean(TransposedMatrixCounts[,i])
  FC <-  (TransposedMatrixCounts[,i]+1)/(mean_col+1)
  FC_MatrixCounts <- cbind(FC_MatrixCounts, FC)
  
  cnt <- cnt + 1
  if(cnt %% 1000 == 0){
    print(cnt)
  }
}
colnames(FC_MatrixCounts) <- rownames(matrix_counts)
# FC_MatrixCounts[is.infinite(FC_MatrixCounts)] <- 0

#log the FC
log2FC_MatrixCounts <- log2(FC_MatrixCounts)

# write.table(FC_MatrixCounts, "Desktop/CBPP_22012021/Lars_Lab/single_cells/CIELAB_umap/FC_MatrixCounts.tsv", quote = F, sep = "\t", row.names = T, col.names = T)

# Use log2FC_MatrixCounts as input in Seurat pipeline

# 1 
log2FC_MatrixCounts_tc <- log2FC_MatrixCounts
# log2FC_MatrixCounts <- log2FC_MatrixCounts_tc
# write.table(log2FC_MatrixCounts_tc, "log2FC_MatrixCounts_tc.tsv", quote = F, row.names = T, col.names = T, sep = "\t")
# 2
log2FC_MatrixCounts_ct <- log2FC_MatrixCounts
# log2FC_MatrixCounts <- log2FC_MatrixCounts_ct
# write.table(log2FC_MatrixCounts_ct, "log2FC_MatrixCounts_ct.tsv", quote = F, row.names = T, col.names = T, sep = "\t")

# 3
log2FC_MatrixCounts_myoblasts <- log2FC_MatrixCounts
# write.table(log2FC_MatrixCounts_myoblasts, "log2FC_MatrixCounts_myoblasts.tsv", quote = F, row.names = T, col.names = T, sep = "\t")

# 4
log2FC_MatrixCounts_osteoclasts <- log2FC_MatrixCounts
# write.table(log2FC_MatrixCounts_osteoclasts, "log2FC_MatrixCounts_osteoclasts.tsv", quote = F, row.names = T, col.names = T, sep = "\t")


