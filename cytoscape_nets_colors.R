# HexFromLabData[[ UniqueNodes07 ]]
HashmapDataFrame <- HexFromLabData$data.frame()

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



# network_thrs065_ALL <- subset(network_thrs04_ALL, network_thrs04_ALL$Weight>=0.65)
# write.table(network_thrs065_ALL, file = "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs065_ALL.tsv",
#             quote = F, row.names = F, col.names = T, sep = "\t")

network_thrs06_ALL <- subset(network_thrs04_ALL, network_thrs04_ALL$Weight>=0.6)
UniqueNodes06_ALL <- unique(c(unique(network_thrs06_ALL$Source), unique(network_thrs06_ALL$Target)))
NodesColors06ALL <- HashmapDataFrame[HashmapDataFrame$Keys %in% UniqueNodes06_ALL, ]

write.table(network_thrs06_ALL, file = "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs06_ALL.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")



a <- HexFromLabData[[ UniqueNodes07 ]]
a <- cbind(UniqueNodes07, a )

a[,2][!(is.na(a[,2]))] <- 35
a[is.na(a)] <- 0

write.table(NodesColors06ALL, file = "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs06_colored_nodes_ALL.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")

write.table(a, "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs07_SIZE.tsv", quote = F, 
            col.names = T, row.names = F, sep="\t")
