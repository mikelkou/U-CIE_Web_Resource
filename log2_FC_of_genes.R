# Droplet-based retina dataset from Macosko et al. (2015),

# BiocManager::install("scRNAseq")
library(scRNAseq)
library(Matrix)

# matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_time_course_ec.csv"))
matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_cell_type_ec.csv"))
rownames(matrix_counts) <- matrix_counts[,1]
matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
matrix_counts <- data.matrix(matrix_counts) # must be a matrix object!

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