# Analysis of genes (sd, bulb genes, etc)

# Matrix with SD of each gene to exclude those with SD<1
SDMoreThanOnelog2FC_MatrixCounts <- c()
for(i in 1:ncol(log2FC_MatrixCounts)){
  std <- sd(log2FC_MatrixCounts[,i])
  SDMoreThanOnelog2FC_MatrixCounts <- cbind(SDMoreThanOnelog2FC_MatrixCounts, std)
}
colnames(SDMoreThanOnelog2FC_MatrixCounts) <- colnames(log2FC_MatrixCounts)
# length(which(SDMoreThanOnelog2FC_MatrixCounts[1,]<1))

SDMoreThanOnelog2FC_MatrixCounts <- select(as.data.frame(log2FC_MatrixCounts), 
                                           -c(names(which(SDMoreThanOnelog2FC_MatrixCounts[1,]<1))))

head(colnames(SDMoreThanOnelog2FC_MatrixCounts))
any(names(which(SDMoreThanOnelog2FC_MatrixCounts[1,]<1)) %in% colnames(SDMoreThanOnelog2FC_MatrixCounts))

# Use SDMoreThanOnelog2FC_MatrixCounts as input in Seurat pipeline
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# Which part of the plot to KEEP based on the coordinates
MatrixToKeep <- umap_dist[which(umap_dist[,1] > (-5) & umap_dist[,2] > (-3.5)), ]
# MatrixToKeep <- umap_dist[which(MatrixToKeep[,2] > (-5) & MatrixToKeep[,2] < 12), ]

VectorOfGenesToKeep <- rownames(MatrixToKeep)
RegulatedGenes <- log2FC_MatrixCounts[,VectorOfGenesToKeep]
# RegulatedGenes <- log2FC_MatrixCounts[,colnames(log2FC_MatrixCounts)[VectorOfGenesToKeep] ]
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# Which part of the plot to REMOVE based on the coordinates
BulbGenes <- which(umap_dist[,1] < 9.0 & umap_dist[,1] > 7.0) #5047 from df and logcounts
BulbGenes <- umap_dist[BulbGenes , ] 
BulbGenes <- which(BulbGenes[,2] < 1 & BulbGenes[,2] >  (-2)) #5006
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# Remove bulb and normalize with std
NotRegulatedGenes <- df[ , colnames(df)[bulb_genes_50_comp]]
# NotRegulatedGenes <- matrix_counts[rownames(matrix_counts)[bulb_genes_seurat], ]
library(dplyr)
log2FC_MatrixCountsWithoutBulbGenes <- select(as.data.frame(log2FC_MatrixCounts), -c(rownames(NotRegulatedGenes)))
any(colnames(NotRegulatedGenes) %in% colnames(log2FC_MatrixCountsWithoutBulbGenes)) # Check
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# Normalization with SD of the matrix AFTER the REMOVAL of the bulb genes
SDNormalizationlog2FC_MatrixCountsWithoutBulbGenes <- c()
for(i in 1:ncol(log2FC_MatrixCountsWithoutBulbGenes)){
  std_col <- sd(log2FC_MatrixCountsWithoutBulbGenes[,i])
  std_norm <-  log2FC_MatrixCountsWithoutBulbGenes[,i]/std_col
  SDNormalizationlog2FC_MatrixCountsWithoutBulbGenes <- cbind(SDNormalizationlog2FC_MatrixCountsWithoutBulbGenes, std_norm)
}
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# Distributions of matrix without bulb genes and bulb genes
log2FC_MatrixCountsWithoutBulbGenes_sd_per_gene <- c()
for(i in 1:ncol(log2FC_MatrixCountsWithoutBulbGenes)){
  std <- sd(log2FC_MatrixCountsWithoutBulbGenes[,i])
  log2FC_MatrixCountsWithoutBulbGenes_sd_per_gene <- cbind(log2FC_MatrixCountsWithoutBulbGenes_sd_per_gene, std)
}

log2FC_BulbGenes_sd_per_gene <- c()
for(i in 1:ncol(NotRegulatedGenes)){
  std <- sd(NotRegulatedGenes[,i])
  log2FC_BulbGenes_sd_per_gene <- cbind(log2FC_BulbGenes_sd_per_gene, std)
}


log2FC_MatrixCountsWithoutBulbGenes_sd_per_gene <- t(log2FC_MatrixCountsWithoutBulbGenes_sd_per_gene)
log2FC_BulbGenes_sd_per_gene <- t(log2FC_BulbGenes_sd_per_gene)

#1
p1 <- hist(log2FC_MatrixCountsWithoutBulbGenes_sd_per_gene)                     # centered at 4
p2 <- hist(log2FC_BulbGenes_sd_per_gene)                     # centered at 6
plot( p1, col=rgb(0,0,1,1/4))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), add=T)  # second

#2
plot(density(log2FC_MatrixCountsWithoutBulbGenes_sd_per_gene[,1]),col="red",xlim=c(-0.5,6),ylim=c(-0.5,2.5))
lines(density(log2FC_BulbGenes_sd_per_gene[,1]))

