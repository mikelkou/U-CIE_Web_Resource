# Droplet-based retina dataset from Macosko et al. (2015),

# BiocManager::install("scRNAseq")
library(scRNAseq)
library(Matrix)

# matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_time_course_ec.csv"))
matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_cell_type_ec.csv"))
rownames(matrix_counts) <- matrix_counts[,1]
matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
matrix_counts <- data.matrix(matrix_counts) # must be a matrix object!

transp <- t(matrix_counts)

cnt <- 0
columns_genes_transp_plus_one <- c()
for(i in 1:ncol(transp)){
  mean_col <- mean(transp[,i])
  FC <-  (transp[,i]+1)/(mean_col+1)
  columns_genes_transp_plus_one <- cbind(columns_genes_transp_plus_one, FC)
  
  cnt <- cnt + 1
  if(cnt %% 1000 == 0){
    print(cnt)
  }
}
colnames(columns_genes_transp_plus_one) <- colnames(t(matrix_counts))
# columns_genes_transp_plus_one[is.infinite(columns_genes_transp_plus_one)] <- 0

write.table(columns_genes_transp_plus_one, "Desktop/CBPP_22012021/Lars_Lab/single_cells/CIELAB_umap/columns_genes_transp_plus_one.tsv", quote = F, sep = "\t", row.names = T, col.names = T)

