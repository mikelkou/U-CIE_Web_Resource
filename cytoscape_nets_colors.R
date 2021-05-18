# HexFromLabData[[ UniqueNodes07 ]]
HashmapDataFrame <- HexFromLabData$data.frame()
write.table(HashmapDataFrame, "NodeColorsSD1ALL.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

#-- Coexpression network pearson (computerome) ---#
coexp_net <- read.delim("/Users/tgn531/Documents/GitHub/CIELAB/coexpr_net_pearson_celltype.tsv", header = T)
UniqueNodes_pearson <- unique(c(unique(coexp_net_sd1[,1]), unique(coexp_net_sd1[,2])))
NodesColors_pearson <- HashmapDataFrame[HashmapDataFrame$Keys %in% UniqueNodes_pearson, ]

coexp_net_sd1 <- coexp_net[coexp_net$gene_1 %in% HashmapDataFrame$Keys, ]
coexp_net_sd1 <- coexp_net[coexp_net$gene_2 %in% HashmapDataFrame$Keys, ]
write.table(coexp_net_sd1, "coexp_net_pearson_sd1.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

#---#
network_thrs07 <- read.delim("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs07_ncol100.tsv", header = F)
UniqueNodes07 <- unique(c(unique(network_thrs07$V1), unique(network_thrs07$V2)))
NodesColors07 <- HashmapDataFrame[HashmapDataFrame$Keys %in% UniqueNodes07, ]
#---#
network_thrs04_sd1 <- read.delim("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_GSE_cell_type_sum_duplicated_ensebl_ids/hPSCs_cell_type_network_thrs04_ncol100_SD1.tsv", 
                                 header = F)
UniqueNodes04_SD1 <- unique(c(unique(network_thrs04_sd1$V1), unique(network_thrs04_sd1$V2)))
NodesColors04SD1 <- HashmapDataFrame[HashmapDataFrame$Keys %in% UniqueNodes04_SD1, ]
#---#
network_thrs04_ALL <- read.delim("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs04_ncol100.tsv", 
                                 header = F)
UniqueNodes04_ALL <- unique(c(unique(network_thrs04_ALL$V1), unique(network_thrs04_ALL$V2)))
NodesColors04ALL <- HashmapDataFrame[HashmapDataFrame$Keys %in% UniqueNodes04_ALL, ]
# write.table(network_thrs04_ALL, "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs04_ncol100.tsv",
#             quote = F, row.names = F, col.names = T, sep = "\t")


network_thrs06_ALL <- subset(network_thrs04_ALL, network_thrs04_ALL$Weight>=0.6)
UniqueNodes06_ALL <- unique(c(unique(network_thrs06_ALL$Source), unique(network_thrs06_ALL$Target)))
NodesColors06ALL <- HashmapDataFrame[HashmapDataFrame$Keys %in% UniqueNodes06_ALL, ]

write.table(network_thrs06_ALL, file = "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs06_ALL.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")


write.table(NodesColors06ALL, file = "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs06_colored_nodes_ALL.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")

write.table(a, "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs07_SIZE.tsv", quote = F, 
            col.names = T, row.names = F, sep="\t")
