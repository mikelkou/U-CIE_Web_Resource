cielab <- function(umap_dist) {
  library(ggplot2)
  # Size of dimensions
  UMAP1Size <- max(umap_dist[, 1]) - min(umap_dist[, 1])
  UMAP2Size <- max(umap_dist[, 2]) - min(umap_dist[, 2])
  UMAP3Size <- max(umap_dist[, 3]) - min(umap_dist[, 3])
  
  # Find the smallest and use it as L
  
  umap_dist_scaled <-
    matrix(nrow = nrow(umap_dist), ncol = ncol(umap_dist))
  if (UMAP1Size < UMAP2Size & UMAP1Size < UMAP3Size) {
    MaxScalingFactor_1 <- 100 / UMAP1Size
    MaxScalingFactor_2 <- 256 / UMAP2Size
    MaxScalingFactor_3 <- 256 / UMAP3Size
    
    MaxScalingFactor <-
      ifelse(
        MaxScalingFactor_1 < MaxScalingFactor_2,
        MaxScalingFactor_1,
        MaxScalingFactor_2
      )
    MaxScalingFactor <-
      ifelse(MaxScalingFactor < MaxScalingFactor_3,
             MaxScalingFactor,
             MaxScalingFactor_3)
    
    #offset
    UMAP1offset <- (max(umap_dist[, 1]) + min(umap_dist[, 1])) / 2
    UMAP2offset <- (max(umap_dist[, 2]) + min(umap_dist[, 2])) / 2
    UMAP3offset <- (max(umap_dist[, 3]) + min(umap_dist[, 3])) / 2
    
    umap_dist_scaled[, 1] <-
      (umap_dist[, 1] - UMAP1offset) * MaxScalingFactor + 50
    umap_dist_scaled[, 2] <-
      (umap_dist[, 2] - UMAP2offset) * MaxScalingFactor
    umap_dist_scaled[, 3] <-
      (umap_dist[, 3] - UMAP3offset) * MaxScalingFactor
  }
  if (UMAP2Size < UMAP1Size & UMAP2Size < UMAP3Size) {
    MaxScalingFactor_1 <- 256 / UMAP1Size
    MaxScalingFactor_2 <- 100 / UMAP2Size
    MaxScalingFactor_3 <- 256 / UMAP3Size
    
    MaxScalingFactor <-
      ifelse(
        MaxScalingFactor_1 < MaxScalingFactor_2,
        MaxScalingFactor_1,
        MaxScalingFactor_2
      )
    MaxScalingFactor <-
      ifelse(MaxScalingFactor < MaxScalingFactor_3,
             MaxScalingFactor,
             MaxScalingFactor_3)
    
    #offset
    UMAP1offset <- (max(umap_dist[, 1]) + min(umap_dist[, 1])) / 2
    UMAP2offset <- (max(umap_dist[, 2]) + min(umap_dist[, 2])) / 2
    UMAP3offset <- (max(umap_dist[, 3]) + min(umap_dist[, 3])) / 2
    
    umap_dist_scaled[, 1] <-
      (umap_dist[, 2] - UMAP2offset) * MaxScalingFactor + 50
    umap_dist_scaled[, 2] <-
      (umap_dist[, 1] - UMAP1offset) * MaxScalingFactor
    umap_dist_scaled[, 3] <-
      (umap_dist[, 3] - UMAP3offset) * MaxScalingFactor
  }
  if (UMAP3Size < UMAP1Size & UMAP3Size < UMAP2Size) {
    MaxScalingFactor_1 <- 256 / UMAP1Size
    MaxScalingFactor_2 <- 256 / UMAP2Size
    MaxScalingFactor_3 <- 100 / UMAP3Size
    
    MaxScalingFactor <-
      ifelse(
        MaxScalingFactor_1 < MaxScalingFactor_2,
        MaxScalingFactor_1,
        MaxScalingFactor_2
      )
    MaxScalingFactor <-
      ifelse(MaxScalingFactor < MaxScalingFactor_3,
             MaxScalingFactor,
             MaxScalingFactor_3)
    
    #offset
    UMAP1offset <- (max(umap_dist[, 1]) + min(umap_dist[, 1])) / 2
    UMAP2offset <- (max(umap_dist[, 2]) + min(umap_dist[, 2])) / 2
    UMAP3offset <- (max(umap_dist[, 3]) + min(umap_dist[, 3])) / 2
    
    umap_dist_scaled[, 1] <-
      (umap_dist[, 3] - UMAP3offset) * MaxScalingFactor  + 50
    umap_dist_scaled[, 2] <-
      (umap_dist[, 1] - UMAP1offset) * MaxScalingFactor
    umap_dist_scaled[, 3] <-
      (umap_dist[, 2] - UMAP2offset) * MaxScalingFactor
  }
  rownames(umap_dist_scaled) <- rownames(umap_dist)
  colnames(umap_dist_scaled) <- c("L", "a", "b")
  Lab <- umap_dist_scaled
  
  Lab <- round(Lab, 2)
  rawdata = structure(
    list(
      Lstar = c(Lab[, 1]),
      Astar = c(Lab[, 2]),
      Bstar = c(Lab[, 3])
    ),
    .Names = c("Lstar", "Astar", "Bstar"),
    row.names = c(rownames(umap_dist)),
    class = "data.frame"
  )
  
  library(colorspace)
  LABdata <- with(rawdata, LAB(Lstar, Astar, Bstar))
  
  ggplot(rawdata, aes(x = umap_dist[, 1], y = umap_dist[, 2])) +
    geom_point(size = 0.5, aes(colour = hex(LABdata, fix = TRUE))) +
    scale_color_identity()
  
  
  library(plotly)
  a <- hex(LABdata, fix = TRUE)
  
  axx <- list(nticks = 4,
              title = "x")
  
  axy <- list(nticks = 4,
              title = "y")
  
  axz <- list(nticks = 4,
              title = "z")
  
  fig <- plot_ly(
    data = as.data.frame(umap_dist),
    x = umap_dist[, 1],
    y = umap_dist[, 2],
    z = umap_dist[, 3],
    type = "scatter3d",
    mode = "markers",
    text = c(rownames(umap_dist)),
    # hoverinfo = 'text',
    marker = list(
      color = hex(LABdata, fix = TRUE),
      size = 3,
      width = 2
    )
  )
  # marker = list(size = 1, width=1)) # controls size of points
  fig <- fig %>% layout(scene = list(
    xaxis = axx,
    yaxis = axy,
    zaxis = axz
  ))
  
  
  library(hashmap)
  HexFromLabData <- hashmap(keys = rownames(umap_dist_scaled), values = hex(LABdata, fix = TRUE))
  
  return(fig)
}
# HexFromLabData[[ UniqueNodes07 ]]
HashmapDataFrame <- HexFromLabData$data.frame()

network_thrs07 <- read.delim("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs07_ncol100.tsv", header = F)
UniqueNodes07 <- unique(c(unique(network_thrs07$V1), unique(network_thrs07$V2)))
NodesColors07 <- HashmapDataFrame[HashmapDataFrame$Keys %in% UniqueNodes07, ]


write.table(NodesColors07, file = "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs07_colored_nodes.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")

write.table(a, "/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/scdata_results/hPSCs_cell_type_network_thrs07_ncol100.tsv", quote = F, 
            col.names = T, row.names = F, sep="\t")

