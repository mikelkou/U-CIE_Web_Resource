matrix_counts <- data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_cell_type_ec.csv"))
rownames(matrix_counts) <- matrix_counts[,1]
matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]

pearson <- cor(t(matrix_counts), method = "pearson")
rownames(pearson) <- rownames(matrix_counts)
colnames(pearson) <- rownames(matrix_counts)
pearson[is.na(pearson)] <- 0

# write.table(pearson, "pearson_matrix_GSE75748_sc_cell_type_ec.tsv", quote = F, 
#             sep = "\t", row.names = T, col.names = T)
# pearson <- read.delim("pearson_matrix_GSE75748_sc_cell_type_ec.tsv")
pearson[is.na(pearson)] <- 0


library(reshape2)
pearson_net <- melt(pearson, c("gene_1", "gene_2"))

pearson_net <- subset(pearson_net, pearson_net$value>=0.7)

library(igraph)
g <- graph_from_data_frame(pearson_net)
frucht_3D <- layout.fruchterman.reingold(g, dim=3)
rownames(frucht_3D) <- rownames(pearson)
umap_dist <- fucht_3D
rownames(umap_dist) <- rownames(frucht_3D)

#------------------------------------------------------------------------------#
# Tree networks - MCL problematic
library("igraph")
edges1 = matrix(c(1,3,2,3,3,4,4,5,4,6),byrow = TRUE, ncol = 2)
g1 = graph_from_edgelist(edges1, directed = FALSE)
vertex_attr(g1, name = "name") = 1:6
plot(g1, vertex.size = 25, edge.width = 5, vertex.color = "coral")

library(treemap)
library(data.tree)
data(GNI2014)
GNI2014$pathString <- paste("world", 
                            GNI2014$continent, 
                            GNI2014$country, 
                            sep = "/")
population <- as.Node(GNI2014)
print(population, "iso3", "population", "GNI", limit = 20)
#---#
data(USArrests)
dd <- dist(scale(USArrests), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc, labels = rownames(USArrests), hang = 0.1, 
     main = "Cluster dendrogram", sub = NULL,
     xlab = NULL, ylab = "Height")

g1 <- graph_from_adjacency_matrix(dd)
V(g1)$name <- rownames(USArrests)
plot(g1)

tree_net <- get.edgelist(g1)
colnames(tree_net) <- c("Source", "Target")

write.table(tree_net, "tree_net.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(HashmapDataFrame, "tree_net_cielab_colors.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

plot(g1, layout.fruchterman.reingold(g1, dim=2))

coords <- layout.fruchterman.reingold(g1, dim=3)
umap_dist <- coords
rownames(umap_dist) <- names(V(g1))


#---------#
usa_pearson <- cor(t(USArrests), method = "pearson")







