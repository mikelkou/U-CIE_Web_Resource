# Droplet-based retina dataset from Macosko et al. (2015),

# BiocManager::install("scRNAseq")
library(scRNAseq)
library(Matrix)
# sce <- MacoskoRetinaData()
#sce <- sce[1:100,1:150]
#matrix_counts <- as.matrix(counts(sce))

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


columns_genes_transp_plus_one_log <- log2(columns_genes_transp_plus_one)
# columns_genes_transp_plus_one_log[is.infinite(columns_genes_transp_plus_one_log)] <- 0

#--- Wrong group meeting implementation ---#
# df_plus_one <- c()
# for(i in 1:nrow(transp)){
#   mean_row <- mean(transp[i,])
#   log2FC <-  transp[i,]+1/mean_row+1
#   df_plus_one <- rbind(df_plus_one, log2FC)
# }
# df_plus_one_log <- log2(df_plus_one)
#--- --------------------------------- ---#

sce <- SingleCellExperiment(assays = list(counts = columns_genes_transp_plus_one_log))
# sce <- sce[, colSums(counts(sce)) > 0] #if error

# Quality control.
library(scater)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization.
# sce <- scater::logNormCounts(sce, log=FALSE) # if log2 already implemented by me
# sce <- scater::logNormCounts(sce)
# normcounts(sce)[1:5,1:5]
logcounts(sce) <- counts(sce)

# Feature selection.
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# Dimensionality reduction.
set.seed(1234)
sce <- runPCA(sce, ncomponents=50, subset_row=hvg)
plotPCA(sce)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

# Clustering.
g <- buildSNNGraph(sce, use.dimred = 'PCA')
colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)

# Visualization.
plotUMAP(sce, colour_by="label")
# head(reducedDim(sce, "UMAP"))

coords_umap_sce <- reducedDim(sce, "PCA")[,3]
# umap_dist <- coords_umap_sce
# umap_dist <- reducedDim(sce, "UMAP")
umap_dist <- cbind(reducedDim(sce, "UMAP"), coords_umap_sce)


#----- Seurat for 3D -----#
#1 initial dataset
sce_seurat_initial_dataset <- RunUMAP(t(columns_genes_transp_plus_one_log), dims = 1:50, n.components = 3L)
umap_dist <- Embeddings(object = sce_seurat_initial_dataset, reduction = "umap")
rownames(umap_dist) <- colnames(columns_genes_transp_plus_one_log)

#2 removed bulb
sce_seurat_new_dataset <- RunUMAP(t(df_new), dims = 1:50, n.components = 3L)
umap_dist <- Embeddings(object = sce_seurat_new_dataset, reduction = "umap")
rownames(umap_dist) <- colnames(df_new)

# removed SD<1
columns_genes_transp_plus_one_log_remove_std_less_than_1_df <- RunUMAP(t(columns_genes_transp_plus_one_log_remove_std_less_than_1), dims = 1:50, n.components = 3L)
umap_dist <- Embeddings(object = columns_genes_transp_plus_one_log_remove_std_less_than_1_df, reduction = "umap")
rownames(umap_dist) <- colnames(columns_genes_transp_plus_one_log_remove_std_less_than_1)



#-- CIELAB coloring --#
# library(uwot)
# # scale = T -->  Scale each column to zero mean and variance 1.
# set.seed(123)
# transp <- t(matrix_counts)
# umap_dist <- umap(transp, n_neighbors = 15, n_components = 3, metric = "euclidean",
#                   learning_rate = 10, init = "spectral",
#                   spread = 10, min_dist = 0.1, nn_method = "fnn", local_connectivity = 2000
# )


# #-- new code T-SNE --#
# x2m <- function(X) {
#   if (!methods::is(X, "matrix")) {
#     m <- as.matrix(X[, which(vapply(X, is.numeric, logical(1)))])
#   }
#   else {
#     m <- X
#   }
#   m
# }
# 
# s1k_tsne <- Rtsne::Rtsne(x2m(umap_dist), perplexity = 15, initial_dims = 759,
#                          partial_pca = F, exaggeration_factor = 4)
# 
# embed_img <- function(X, Y, k = 15, ...) {
#   args <- list(...)
#   args$coords <- Y
#   args$x <- X
#   
#   do.call(vizier::embed_plot, args)
# }
# 
# embed_img(s1k_tsne$Y, s1k_tsne$Y, pc_axes = T, equal_axes = T, alpha_scale = 0.5, title = "time course UMAP", cex = 1)
# 
#-------------------------------------------------------------------------------#

# original umap
UMAP1Size <- max(umap_dist[,1]) - min(umap_dist[,1]) #a
UMAP2Size <- max(umap_dist[,2]) - min(umap_dist[,2]) #b
UMAP3Size <- max(umap_dist[,3]) - min(umap_dist[,3]) #L

MaxScalingFactor_1 <- 256/UMAP1Size
MaxScalingFactor_2 <- 256/UMAP2Size
MaxScalingFactor_3 <- 100/UMAP3Size

MaxScalingFactor <- ifelse(MaxScalingFactor_1 < MaxScalingFactor_2, MaxScalingFactor_1, MaxScalingFactor_2)
MaxScalingFactor <- ifelse(MaxScalingFactor < MaxScalingFactor_3, MaxScalingFactor, MaxScalingFactor_3)

#offset 
UMAP1offset <- ( max(umap_dist[,1]) + min(umap_dist[,1]) )/2
UMAP2offset <- ( max(umap_dist[,2]) + min(umap_dist[,2]) )/2
UMAP3offset <- ( max(umap_dist[,3]) + min(umap_dist[,3]) )/2

umap_dist_scaled <- matrix(nrow=nrow(umap_dist), ncol=ncol(umap_dist))
umap_dist_scaled[,1] <- (umap_dist[,1] - UMAP1offset)*MaxScalingFactor 
umap_dist_scaled[,2] <- (umap_dist[,2] - UMAP2offset)*MaxScalingFactor 
umap_dist_scaled[,3] <- (umap_dist[,3] - UMAP3offset)*MaxScalingFactor + 50
rownames(umap_dist_scaled) <- rownames(umap_dist)

# Lab values
UMAPSize <- ifelse(UMAP1Size < UMAP2Size, UMAP1Size, UMAP2Size)
if(UMAPSize > UMAP3Size){
  print("312")
  Lab <- umap_dist_scaled[,c(3,1,2)]
} else {
    if(UMAP1Size > UMAP2Size){
      print("132")
      Lab <- umap_dist_scaled[,c(1,3,2)]
    } else{
      print("123")
        Lab <- umap_dist_scaled[,c(1,2,3)]
  }
}
colnames(Lab) <- c("L", "a", "b")

Lab <- round(Lab,2)
rawdata=structure(list(Lstar = c(Lab[,1]), Astar = c(Lab[,2]),
                       Bstar = c(Lab[,3])), .Names = c("Lstar","Astar", "Bstar"), row.names = c(rownames(umap_dist)), class = "data.frame")

library(colorspace)
LABdata <- with(rawdata,LAB(Lstar,Astar,Bstar))

ggplot(rawdata, aes(x=umap_dist[,1], y=umap_dist[,2])) +
  geom_point(size=0.5, aes(colour=hex(LABdata,fix = TRUE))) +
  scale_color_identity()


library(plotly)
a <- hex(LABdata,fix = TRUE)

axx <- list(
  nticks = 4,
  title = "x"
)

axy <- list(
  nticks = 4,
  title = "y"
)

axz <- list(
  nticks = 4,
  title = "z"
)

fig <- plot_ly(data = as.data.frame(umap_dist), 
               x = umap_dist[,1], y = umap_dist[,2], z = umap_dist[,3], 
               type = "scatter3d", 
               mode = "markers", 
               text = c(rownames(umap_dist)),
               # hoverinfo = 'text',
               marker = list(color = hex(LABdata,fix = TRUE), size = 3, width=2))
               # marker = list(size = 1, width=1)) # controls size of points
fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

fig

#------------------------------------------------------------------------------#

to_keep <- umap_dist[which(umap_dist[,2] > (-6) & umap_dist[,1] < 2), ]
# to_keep <- umap_dist[which(to_keep[,2] > (-5) & to_keep[,2] < 12), ]

kkk <- rownames(to_keep)
regulated_genes <- columns_genes_transp_plus_one_log[,kkk]
# regulated_genes <- columns_genes_transp_plus_one_log[,colnames(columns_genes_transp_plus_one_log)[kkk] ]



bulb_genes_50_comp_x <- which(umap_dist[,1] < 9.0 & umap_dist[,1] > 7.0) #5047 from df and logcounts
bulb_genes_50_comp_x <- umap_dist[bulb_genes_50_comp_x , ] 
bulb_genes_50_comp <- which(bulb_genes_50_comp_x[,2] < 1 & bulb_genes_50_comp_x[,2] >  (-2)) #5006

not_regulated_genes <- df[,colnames(df)[bulb_genes_50_comp] ]

# Seurat
bulb_genes_seurat <- which(umap_dist[,1] < (-6)) #5047 from df and logcounts
# bulb_genes_seurat <- umap_dist[bulb_genes_seurat , ] 

not_regulated_genes <- matrix_counts[rownames(matrix_counts)[bulb_genes_seurat], ]
library(dplyr)
df_new <- select(as.data.frame(columns_genes_transp_plus_one_log), -c(rownames(not_regulated_genes)))
any(colnames(not_regulated_genes) %in% colnames(df_new))





# Remove bulb and normalize with std
rm <- colnames(not_regulated_genes)
library(dplyr)
df_new <- select(as.data.frame(df_plus_one_log), -c(colnames(not_regulated_genes)))
any(colnames(not_regulated_genes) %in% colnames(df_new))

columns_genes_transp_plus_one_log_sd <- c()
for(i in 1:ncol(columns_genes_transp_plus_one_log)){
  std <- sd(columns_genes_transp_plus_one_log[,i])
  columns_genes_transp_plus_one_log_sd <- cbind(columns_genes_transp_plus_one_log_sd, std)
}
colnames(columns_genes_transp_plus_one_log_sd) <- colnames(columns_genes_transp_plus_one_log)
length(which(columns_genes_transp_plus_one_log_sd[1,]<1))

columns_genes_transp_plus_one_log_remove_std_less_than_1 <- select(as.data.frame(columns_genes_transp_plus_one_log), -c(names(which(columns_genes_transp_plus_one_log_sd[1,]<1))))

# names(which(columns_genes_transp_plus_one_log_sd[1,]<1))
any(names(which(columns_genes_transp_plus_one_log_sd[1,]<1)) %in% colnames(columns_genes_transp_plus_one_log_remove_std_less_than_1))

#---------

# Normalize new df with the sd and re-run umap
# df_new_norm <- c()
# for(i in 1:nrow(df_new)){
#   std_row <- sd(df_new[i,])
#   std_norm <-  df_new[i,]/std_row
#   df_new_norm <- rbind(df_new_norm, std_norm)
# }

# Find distributions of new df and bulb
df_new_sd_per_gene <- c()
for(i in 1:ncol(df_new)){
  std <- sd(df_new[,i])
  df_new_sd_per_gene <- cbind(df_new_sd_per_gene, std)
}

bulb_sd_per_gene <- c()
for(i in 1:ncol(not_regulated_genes)){
  std <- sd(not_regulated_genes[,i])
  bulb_sd_per_gene <- cbind(bulb_sd_per_gene, std)
}


df_new_sd_per_gene <- t(df_new_sd_per_gene)
# df_new_sd_per_gene_2 <- cbind(rep("df_new", nrow(df_new_sd_per_gene)),df_new_sd_per_gene)
bulb_sd_per_gene <- t(bulb_sd_per_gene)
# bulb_sd_per_gene_2 <- cbind(rep("bulb", nrow(bulb_sd_per_gene)),bulb_sd_per_gene)


p1 <- hist(df_new_sd_per_gene)                     # centered at 4
p2 <- hist(bulb_sd_per_gene)                     # centered at 6
plot( p1, col=rgb(0,0,1,1/4))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), add=T)  # second

plot(density(df_new_sd_per_gene[,1]),col="red",xlim=c(-0.5,6),ylim=c(-0.5,2.5))
lines(density(bulb_sd_per_gene[,1]))


#-------------------------------------------------------------------------------#
#-- t-sne
# fig <- plot_ly(data = as.data.frame(umap_dist), 
#                x = s1k_tsne$Y[,1], y = s1k_tsne$Y[,2], z = s1k_tsne$costs, 
#                type = "scatter3d", 
#                mode = "markers", 
#                marker = list(color = hex(LABdata,fix = TRUE)),
#                marker = list(size = 1, width=1)) # controls size of points
# fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
# 
# fig
#-------------------------------------------------------------------------------#
# 
# library(igraph)
# edges <- colnames(matrix_counts)[1:nrow(matrix_counts)]
# nodes <- matrix_counts[,1]
# df_met <- na.omit(data.frame(V1 = edges, V2 = nodes))
# g <- graph_from_data_frame(df_met[1:10,])
# V(g)$color <- hex(LABdata[1:10],fix = TRUE)
# plot(g)
# 
# 
# data(karate, package="igraphdata")
# G <- upgrade_graph(karate)
# L <- layout.circle(G)
# vs <- V(G)
# es <- as.data.frame(get.edgelist(G))
# 
# Nv <- length(vs)
# Ne <- length(es[1]$V1)
# Xn <- L[,1]
# Yn <- L[,2]
# 
# network <- plot_ly(x = ~Xn, y = ~Yn, mode = "markers", text = vs$label, hoverinfo = "text")
# 
# 
# 





