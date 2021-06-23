library(rgl)
library(geometry) # faces
library(ptinpoly) # pip3d
# '1' indicates that the point is contained in the polyhedron.
# '0' indicates that the point lies exactly on the surface of the polyhedron.
# '-1' indicates that the point lies outside the polyhedron.
library(pracma) # distmat

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

NewConvexCloud <- function(S, RotL, Rota, Rotb,  TrL, Tra, Trb, Query){
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
    ConvexCloud[,1] <- 1*ConvexCloud[,1]
    ConvexCloud[,2] <- 1*ConvexCloud[,2]
    ConvexCloud[,3] <- 1.2*ConvexCloud[,3]
  }
  # print(ConvexCloud)
  return(ConvexCloud)
} # Transform the UMAP Convex to otpimize it

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

# Test the Master function
# MasterFunction(param = c(start.values), dat, polygon, faces)

# data points to be fitted
# umap_dist <- read.delim("~/Desktop/umap_dist.tsv")
write.table(umap_dist, "umap_dist.csv", col.names = T, row.names = T, sep = ",")
a <- read.delim("umap_dist.csv")
# cielab(umap_dist)
dat <- UMAPConvex(umap_dist)
polygon <- ColorSpacePolygon(RGB_space)
faces <- convhulln(polygon, return.non.triangulated.facets = T)

# Simplex optimizer

# According to NEW start values every loop
set.seed(123)
start.values <- c(S , pi/4 ,pi/4, pi/4, TrL, Tra , Trb) # S, RotL, Rota, Rotb,  TrL, Tra, Trb
k <- c()
for(i in 1:25){
  Simplex_optim <- optim(par = start.values,
                      method = "Nelder-Mead",
                      MasterFunction,
                      data = dat,
                      polygon = polygon,
                      faces = faces,
                      control=list(fnscale=-1, maxit=1000)) #, trace=T
  k[[i]] <- Simplex_optim
  if((i-1)!=0){
    if(k[[i]]$value >= k[[i-1]]$value){
      start.values <- k[[i]]$par
    } else {
      start.values <- k[[i-1]]$par
      k[[i]] <- k[[i-1]]
    }
    if(i%%2==0){
      start.values <- c(k[[i]]$par[1], k[[i]]$par[2]+(pi/4), k[[i]]$par[3], k[[i]]$par[4], k[[i]]$par[5], k[[i]]$par[6], k[[i]]$par[7]) # mirror rotation in x
    }
  }
  result <- k[[i]]
} # loop around Simplex algorithm
print(result)



# According to INITIAL start values every loop
set.seed(123)
angle <- 1
simplex_vectors <- c()
start.values <- c(S , pi/4 ,pi/4, pi/4, TrL, Tra , Trb) # S, RotL, Rota, Rotb,  TrL, Tra, Trb
k <- c()

# library(parallel)
# library(MASS)

# library(foreach)
# library(doParallel)

numCores <- detectCores()

registerDoParallel(numCores)
# system.time(

# foreach(i=1:25) %dopar% {
for(i in 1:25){
  Simplex_optim <- optim(par = start.values,
                         method = "Nelder-Mead",
                         MasterFunction,
                         data = dat,
                         polygon = polygon,
                         faces = faces,
                         control=list(fnscale=-1, maxit=1000)) #, trace=T
  k[[i]] <- Simplex_optim
  angle <- angle + 0.5
  start.values <- c(start.values[1], start.values[2]+(pi/angle), start.values[3], start.values[4], start.values[5], start.values[6], start.values[7]) # mirror rotation in x
  simplex_vectors <- rbind(simplex_vectors, k[[i]]$par)
# }
}
# stopImplicitCluster()


# )

# library(matrixStats)
start.values <- simplex_vectors[which(simplex_vectors[,1] == max(simplex_vectors[,1]))[1], ]

optional_values <- c()
for(i in 1:nrow(simplex_vectors)){
  vector_dist <- distmat(start.values[1], simplex_vectors[i,1])
  if(vector_dist >= 3) optional_values <- rbind(optional_values, simplex_vectors[i,])
}



# Plot after transformations
# ConvexCloud <- UMAPConvex(umap_dist_scaled)
ConvexCloud <- UMAPConvex(umap_dist)
ConvexCloud <- Scaling(ConvexCloud, start.values[1]*1.2)
ConvexCloud <- Rotation(ConvexCloud,  start.values[2]  , start.values[3] , start.values[4])
ConvexCloud <- Translation(ConvexCloud, start.values[5] , start.values[6] ,start.values[7])

plot(as.matrix(ConvexCloud[,c(1,2)]), pch=19, xlim = c(-100,100), ylim = c(-100,100))
lines(polygon[, c(1,2)], col="blue") # Cielab

colnames(ConvexCloud) <- colnames(RGBtoLabCoords)
a <- rbind(as.data.frame(RGBtoLabCoords), as.data.frame(ConvexCloud))
a$colors <- as.factor(c(rep(1, nrow(RGBtoLabCoords)), rep(0, nrow(ConvexCloud))))

library(plotly)
fig <- plot_ly(a, x = a[,1], y = a[,2], z = a[,3], color = a[,4], colors = c('#BF382A', '#0C4B8E'), mode='lines+markers',
               line = list(width = 6), marker = list(size = 3.5))
# fig <- fig %>% add_markers() # instead of lines
fig <- fig %>% layout(scene = list(xaxis = list(title = 'L'),
                                   yaxis = list(title = 'a'),
                                   zaxis = list(title = 'b')))

fig

#------------------------------------------------------------------------------#

# Apply transformation in the whole cloud
NewUMAP <- Scaling(umap_dist, start.values[1]*1.2)
NewUMAP <- Rotation(as.matrix(NewUMAP), start.values[2]  , start.values[3] , start.values[4])
NewUMAP <- Translation(NewUMAP, start.values[5] , start.values[6] ,start.values[7])

# Colors to the New UMAP cloud
Lab <- NewUMAP
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
# a <- hex(LABdata, fix = TRUE)

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
fig <- fig %>% layout(scene = list(
  xaxis = axx,
  yaxis = axy,
  zaxis = axz
)) 

# fig <- fig %>% add_trace(y=umap_dist[,1], name='L*')
# fig <- fig %>% add_trace(y=umap_dist[,2], name='a*')
# fig <- fig %>% add_trace(y=umap_dist[,3], name='b*')
fig

cielab(umap_dist)
#



# a <- subset(umap_dist, umap_dist[,2]<5)
# a <- subset(umap_dist, umap_dist[,2]>(-5))

plot(a)

#------ Initial Guess ---------------------------------------------------------#
#--- Translation ---#
# umap_dist <- read.delim("~/Desktop/umap_dist.tsv")
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





#------ Comment out -----------------------------------------------------------#

# # SANN optimizer
# set.seed(123)
# start.values <- c(11.2795 , pi/4 ,pi/4, pi/4, 49.25131 , 4.680764 ,2.274956 ) # S, RotL, Rota, Rotb,  TrL, Tra, Trb
# k <- c()
# for(i in 1:25){
#   SANN_optim <- optim(par = start.values,
#                       method = "SANN", 
#                       MasterFunction, 
#                       data = dat, 
#                       polygon = polygon, 
#                       faces = faces, 
#                       control=list(fnscale=-1, maxit = 10000, temp = 20, tmax = 100)) #, trace=T
#   k[[i]] <- SANN_optim
#   if((i-1)!=0){
#     if(k[[i]]$value >= k[[i-1]]$value){
#     start.values <- k[[i]]$par
#   } else {
#     start.values <- k[[i-1]]$par
#     k[[i]] <- k[[i-1]]
#   }
#     if(i%%2==0){
#       start.values <- c(k[[i]]$par[1], k[[i]]$par[2]+(pi), k[[i]]$par[3], k[[i]]$par[4], k[[i]]$par[5], k[[i]]$par[6], k[[i]]$par[7]) # mirror rotation in x
#     }
#     print(start.values)
#   }
#   result <- k[[i]]
# }


