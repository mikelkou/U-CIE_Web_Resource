library(shiny)
library(shinyjs)
library(shinydashboard)
library(scRNAseq)
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

# library(future)
# library(promises)
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

Distance <- function(S, RotL, Rota, Rotb,  TrL, Tra, Trb, Query, polygon, faces){
  ConvexCloud <- NewConvexCloud(S, RotL, Rota, Rotb,  TrL, Tra, Trb, Query)
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

MasterFunction <- function(param, data, polygon, faces){
  X <- Distance(param[1],param[2],param[3],param[4],param[5],param[6],param[7], data, polygon, faces) # S, RotL, Rota, Rotb,  TrL, Tra, Trb
  a <- 1
  f <- (a*param[1]) - sum(X^2)
  return(f)
}


