# if (!require("DT")) install.packages('DT', dependencies = T) !Done!
# if (!require("config")) install.packages('config', dependencies = T) !Done!
# if (!require("shiny")) install.packages('shiny', dependencies = T) !Done!
# if (!require("shinyjs")) install.packages('shinyjs', dependencies = T) !Done!
# if (!require("shinydashboard")) install.packages('shinydashboard', dependencies = T) !Done!
# if (!require("scRNAseq")) BiocManager::install("scRNAseq") # NOT installed
# if (!require("Matrix")) install.packages('Matrix', dependencies = T) !Done!
# if (!require("dplyr")) install.packages('dplyr', dependencies = T) !Done!
# if (!require("Seurat")) install.packages('Seurat', dependencies = T) !Done!
# if (!require("patchwork")) install.packages('patchwork', dependencies = T) !Done! 
# if (!require("rgl")) install.packages('rgl', dependencies = T) !Done!
# if (!require("geometry")) install.packages('geometry', dependencies = T) !Done!
# if (!require("ptinpoly")) install.packages('ptinpoly', dependencies = T) !Done!
# if (!require("pracma")) install.packages('pracma', dependencies = T) !Done!
# if (!require("plotly")) install.packages('plotly', dependencies = T) !Done!
# if (!require("umap")) install.packages('umap', dependencies = T) !Done!
# if (!require("readxl")) install.packages('readxl', dependencies = T) !Done!
# if (!require("prodlim")) install.packages('prodlim', dependencies = T) !Done!
# if (!require("shinyWidgets")) install.packages('shinyWidgets', dependencies = T) !Done!
library(shinybusy)
library(shiny)
library(shinyjs)
library(shinydashboard)
# library(scRNAseq)
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(rgl)
library(geometry) # faces
library(ptinpoly) # pip3d
# '1' indicates that the point is contained in the polyhedron.
# '0' indicates that the point lies exactly on the surface of the polyhedron.
# '-1' indicates that the point lies outside the polyhedron.
library(pracma) # distmat
library(plotly)
library(DT)
library(umap) # for distance matrix
library(readxl)
library(prodlim)
library(shinyWidgets)
library(colorspace)

library(timecourse) # time course example for high dimensional https://www.bioconductor.org/packages/devel/bioc/vignettes/timecourse/inst/doc/timecourse.pdf

# library(config)
library(future)
library(promises)
# future::plan(multisession, workers = 8)
future::plan(multisession, workers = 16L)
plan(multicore)
# plan(multisession)
# library(ipc)
# 
# library(foreach)
# library(doParallel)

RGB_space <- data.frame("R"= c(seq(0, 255, by=32),255, # K -> R
                               seq(0, 255, by=32),255, # G -> Y
                               rep(255, 9), # Y -> W
                               rep(255, 9), # R -> Y
                               rep(255, 9), # R -> M
                               seq(0, 255, by=32),255, # B -> M
                               rep(0, 9), # B -> C
                               seq(0, 255, by=32),255, # C -> W
                               rep(0, 9), # G -> C
                               rep(0, 9), # K -> G
                               rep(0, 9), # K -> B
                               rep(255, 9)), # M -> W
                        
                        "G"= c(rep(0, 9), # K -> R
                               rep(255, 9), # G -> Y
                               rep(255, 9), # Y -> W
                               seq(0, 255, by=32),255, # R -> Y
                               rep(0, 9),  # R -> M
                               rep(0, 9),  # B -> M
                               seq(0, 255, by=32),255, # B -> C
                               rep(255, 9), # C -> W
                               rep(255, 9), # G -> C
                               seq(0, 255, by=32),255, # K -> G
                               rep(0, 9), # K -> B
                               seq(0, 255, by=32),255), # M -> W
                        
                        "B"= c(rep(0, 9), # K -> R
                               rep(0, 9), # G -> Y
                               seq(0, 255, by=32),255, # Y -> W
                               rep(0, 9), # R -> Y
                               seq(0, 255, by=32),255, # R -> M
                               rep(255, 9),  # B -> M
                               rep(255, 9), # B -> C
                               rep(255, 9), # C -> W
                               seq(0, 255, by=32),255, # G -> C
                               rep(0, 9), # K -> G
                               seq(0, 255, by=32),255, # K -> B
                               rep(255, 9))# M -> W
                        
)

ColorSpacePolygon <- function(ColorSpace){
  RGB2Lab <- function(mat) {
    # cast mat in case of a single numeric vector
    mat <- matrix(mat, ncol=3)
    # input should be a matrix with R,G and B as columns in [0,1], and columnwise pixels as rows.
    # RGB -> XYZ
    thres1 <- 0.04045
    M <- c(0.412453, 0.357580, 0.180423,
           0.212671, 0.715160, 0.072169,
           0.019334, 0.119193, 0.950227)
    M <- matrix(M, nrow=3, byrow=TRUE)
    matthres <- mat > thres1
    mat <- matthres * ((mat + 0.055) / 1.055) ^ 2.2 + (!matthres) * mat / 12.92
    xyz <- mat %*% t(M)
    
    # XYZ -> Lab
    thres2 <- 0.008856
    xyz <- sweep(xyz, 2, c(0.950456,1,1.088754), "/")
    #yalone <- xyz[,2]
    #y3 <- yalone^(1/3)
    xyzthres <- xyz > thres2
    xyz <- xyzthres * xyz^(1/3) + (!xyzthres) * (7.787*xyz+16/116)
    #L <- xyzthres[,2] * (116*y3-16) + (!xyzthres[,2]) * (903.3*yalone)
    L <- 116 * xyz[,2] - 16
    a <- 500 * (xyz[,1] - xyz[,2])
    b <- 200 * (xyz[,2] - xyz[,3])
    return(cbind(L,a,b))
  }
  RGBtoLabCoords <- as.data.frame(RGB2Lab(as.matrix(ColorSpace)/255)) # RGB space
  ch <- chull(as.matrix(RGBtoLabCoords))
  polygon <- as.matrix(RGBtoLabCoords)[c(ch, ch[1]), ] # Convex that cloud should fit in
  return(polygon)
} # Switch colorspaces -- e.g. Fit sRGB (box) into the CIELab color space

UMAPConvex <- function(Query){
  ch_cloud <- chull(as.matrix(Query))
  ConvexCloud <- as.matrix(Query)[c(ch_cloud, ch_cloud[1]), ] # Convex of cloud
  return(ConvexCloud)
} # Create a Convex Hull from the UMAP

Rotation <- function(ConvexCloud, RotL, Rota, Rotb){
  ConvexCloud <- rotate3d(obj = ConvexCloud, angle = RotL, x = 1, y = 0, z = 0)
  ConvexCloud <- rotate3d(obj = ConvexCloud, angle = Rota, x = 0, y = 1, z = 0)
  ConvexCloud <- rotate3d(obj = ConvexCloud, angle = Rotb, x = 0, y = 0, z = 1)
  return(ConvexCloud)
}

Translation <- function(ConvexCloud, TrL, Tra, Trb){
  ConvexCloud <- translate3d(ConvexCloud, TrL, Tra, Trb)
  return(ConvexCloud)
}

Scaling <- function(ConvexCloud, S){
  ConvexCloud <- scale3d(ConvexCloud, S, S, S)
  return(ConvexCloud)
}

NewConvexCloud <- function(S, RotL, Rota, Rotb,  TrL, Tra, Trb, WL, Wa, Wb, Query){
  ConvexCloud <- Rotation(Query, RotL, Rota, Rotb)
  ConvexCloud <- Scaling(ConvexCloud, S)
  ConvexCloud <- Translation(ConvexCloud, TrL, Tra, Trb)
  
  L <- (max(ConvexCloud[,1]) - min(ConvexCloud[,1]))
  a <- (max(ConvexCloud[,2]) - min(ConvexCloud[,2]))
  b <- (max(ConvexCloud[,3]) - min(ConvexCloud[,3]))
  
  if(b > L & b > a){
    ConvexCloud <- ConvexCloud
  } else {
    # print(S)
    ConvexCloud[,1] <- WL*ConvexCloud[,1]
    ConvexCloud[,2] <- Wa*ConvexCloud[,2]
    ConvexCloud[,3] <- Wb*ConvexCloud[,3]
  }
  return(ConvexCloud)
} # Transform the UMAP Convex to otpimize it

Distance <- function(S, RotL, Rota, Rotb,  TrL, Tra, Trb, WL, Wa, Wb, Query, polygon, faces){
  ConvexCloud <- NewConvexCloud(S, RotL, Rota, Rotb,  TrL, Tra, Trb, WL, Wa, Wb, Query)
  point_in_space <- pip3d(polygon, faces, ConvexCloud)
  outside <- ConvexCloud[which(point_in_space==-1), ]
  # print(outside)
  if(length(outside) == 0){
    dist <- 0
  } else {
    dist_mat <- distmat(polygon, outside)
    dist <- as.data.frame(apply(dist_mat,2,min))
  }
  return(dist)
} # Gives me the distance of the points from the Polygon

MasterFunction <- function(param, WL, Wa, Wb, data, polygon, faces){
  X <- Distance(param[1],param[2],param[3],param[4],param[5],param[6],param[7], WL, Wa, Wb, data, polygon, faces) # S, RotL, Rota, Rotb,  TrL, Tra, Trb
  a <- 1
  f <- (a*param[1]) - sum(X^2)
  return(f)
}

FitColorsFunction <- function(umap_dist, polygon, WL, Wa, Wb){
  dat <- UMAPConvex(umap_dist)
  faces <- convhulln(polygon, return.non.triangulated.facets = T)
  
  #------ Initial Guess ---------------------------------------------------------#
  #--- Translation ---#
  centroidval_color_space <- colMeans(polygon) # centroid of color space
  centroidval_cloud <- colMeans(UMAPConvex(umap_dist)) # centroid of cloud
  
  TrL <- (centroidval_color_space - centroidval_cloud)[1]
  Tra <- (centroidval_color_space - centroidval_cloud)[2]
  Trb <- (centroidval_color_space - centroidval_cloud)[3]
  
  #--- Scaling factor ---#
  ConvexCloud <- (UMAPConvex(umap_dist))
  
  Cloud1Size <- max(ConvexCloud[, 1]) - min(ConvexCloud[, 1])
  Cloud2Size <- max(ConvexCloud[, 2]) - min(ConvexCloud[, 2])
  Cloud3Size <- max(ConvexCloud[, 3]) - min(ConvexCloud[, 3])
  
  # Find the smallest and use it as L
  
  Cloud_scaled <- matrix(nrow = nrow(ConvexCloud), ncol = ncol(ConvexCloud))
  if (Cloud1Size < Cloud2Size & Cloud1Size < Cloud3Size) {
    MaxScalingFactor_1 <- max(polygon[,1]) / Cloud1Size
    MaxScalingFactor_2 <- max(polygon[,2]) / Cloud2Size
    MaxScalingFactor_3 <- max(polygon[,3]) / Cloud3Size
    
    MaxScalingFactor <- ifelse(MaxScalingFactor_1 < MaxScalingFactor_2,MaxScalingFactor_1, MaxScalingFactor_2)
    MaxScalingFactor <- ifelse(MaxScalingFactor < MaxScalingFactor_3, MaxScalingFactor, MaxScalingFactor_3)
    
    #offset
    Cloud1offset <- (max(ConvexCloud[, 1]) + min(ConvexCloud[, 1])) / 2
    Cloud2offset <- (max(ConvexCloud[, 2]) + min(ConvexCloud[, 2])) / 2
    Cloud3offset <- (max(ConvexCloud[, 3]) + min(ConvexCloud[, 3])) / 2
    
    Cloud_scaled[, 1] <- (ConvexCloud[, 1] - Cloud1offset) * MaxScalingFactor + 50
    Cloud_scaled[, 2] <- (ConvexCloud[, 2] - Cloud2offset) * MaxScalingFactor
    Cloud_scaled[, 3] <- (ConvexCloud[, 3] - Cloud3offset) * MaxScalingFactor
  }
  if (Cloud2Size < Cloud1Size & Cloud2Size < Cloud3Size) {
    MaxScalingFactor_1 <- max(polygon[,2]) / Cloud1Size
    MaxScalingFactor_2 <- max(polygon[,1]) / Cloud2Size
    MaxScalingFactor_3 <- max(polygon[,3]) / Cloud3Size
    
    MaxScalingFactor <- ifelse(MaxScalingFactor_1 < MaxScalingFactor_2, MaxScalingFactor_1, MaxScalingFactor_2)
    MaxScalingFactor <- ifelse(MaxScalingFactor < MaxScalingFactor_3, MaxScalingFactor, MaxScalingFactor_3)
    
    #offset
    Cloud1offset <- (max(ConvexCloud[, 1]) + min(ConvexCloud[, 1])) / 2
    Cloud2offset <- (max(ConvexCloud[, 2]) + min(ConvexCloud[, 2])) / 2
    Cloud3offset <- (max(ConvexCloud[, 3]) + min(ConvexCloud[, 3])) / 2
    
    Cloud_scaled[, 1] <- (ConvexCloud[, 2] - Cloud2offset) * MaxScalingFactor + 50
    Cloud_scaled[, 2] <- (ConvexCloud[, 1] - Cloud1offset) * MaxScalingFactor
    Cloud_scaled[, 3] <- (ConvexCloud[, 3] - Cloud3offset) * MaxScalingFactor
  }
  if (Cloud3Size < Cloud1Size & Cloud3Size < Cloud2Size) {
    MaxScalingFactor_1 <- max(polygon[,3]) / Cloud1Size
    MaxScalingFactor_2 <- max(polygon[,2]) / Cloud2Size
    MaxScalingFactor_3 <- max(polygon[,1]) / Cloud3Size
    
    MaxScalingFactor <- ifelse( MaxScalingFactor_1 < MaxScalingFactor_2, MaxScalingFactor_1, MaxScalingFactor_2)
    MaxScalingFactor <- ifelse(MaxScalingFactor < MaxScalingFactor_3, MaxScalingFactor, MaxScalingFactor_3)
    
    #offset
    Cloud1offset <- (max(ConvexCloud[, 1]) + min(ConvexCloud[, 1])) / 2
    Cloud2offset <- (max(ConvexCloud[, 2]) + min(ConvexCloud[, 2])) / 2
    Cloud3offset <- (max(ConvexCloud[, 3]) + min(ConvexCloud[, 3])) / 2
    
    Cloud_scaled[, 1] <- (ConvexCloud[, 3] - Cloud3offset) * MaxScalingFactor  + 50
    Cloud_scaled[, 2] <- (ConvexCloud[, 1] - Cloud1offset) * MaxScalingFactor
    Cloud_scaled[, 3] <- (ConvexCloud[, 2] - Cloud2offset) * MaxScalingFactor
  }
  
  S <- MaxScalingFactor
  
  # Simplex optimizer
  set.seed(123)
  simplex_vectors <- c()
  angle <- 1
  start.values <- c(S , pi/4 ,pi/4, pi/4, TrL, Tra , Trb) # S, RotL, Rota, Rotb,  TrL, Tra, Trb
  k <- c()
  
  for(i in 1:25){
    Simplex_optim <- optim(par = start.values,
                           method = "Nelder-Mead",
                           MasterFunction,
                           WL = WL,
                           Wa = Wa,
                           Wb = Wb,
                           data = dat,
                           polygon = polygon,
                           faces = faces,
                           control=list(fnscale=-1, maxit=1000)) #, trace=T
    k[[i]] <- Simplex_optim
    angle <- angle + 0.5
    start.values <- c(start.values[1], start.values[2]+(pi/angle), start.values[3], start.values[4], start.values[5], start.values[6], start.values[7]) # mirror rotation in x
    simplex_vectors <- rbind(simplex_vectors, k[[i]]$par)
  }
  return(simplex_vectors)
}



convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  if(length(mi$attribs$class)>0 && mi$attribs$class=="treeview"){
    mi$attribs$class=NULL
  }
  mi
}

b64 <- base64enc::dataURI(file = "www/favicon-16x16.png", mime = "image/png")
# b64 <- base64enc::dataURI(file = "~/Desktop/cielab-440.png", mime = "image/png")

