ConvexHullCIELab <- function(ColorSpace,Query){
  # error == 0 --> point inside or on the shape
  # error == 1 --> point outside with distance less than 20
  # error == 2 --> point outside with distance less than 50 and more than 20
  # error == 3 --> point outside with distance more than 50 
  
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
  RGB_space <- as.matrix(RGB_space)
  RGBtoLabCoords <- as.data.frame(RGB2Lab(RGB_space/255)) # RGB space
  ch <- chull(as.matrix(RGBtoLabCoords))
  coords <- as.matrix(RGBtoLabCoords)[c(ch, ch[1]), ]
  plot(as.matrix(RGBtoLabCoords), pch=19)
  lines(coords, col="red")
  
  library(geometry)
  faces <- convhulln(coords,return.non.triangulated.facets = T)
  
  library(ptinpoly)
  # '1' indicates that the point is contained in the polyhedron.
  # '0' indicates that the point lies exactly on the surface of the polyhedron.
  # '-1' indicates that the point lies outside the polyhedron.
  
  # The values in the Faces matrix must be integers with values running from 1 to N, where N is the number of vertices. 
  # A value of '1' in this matrix, for example, represents the 1st vertex, i.e., the vertex defined by the first row in the matrix Vertices.
  query <- data.frame("x" = Query[,1],"y" = Query[,2],"z" = Query[,3])
  
  df <- c()
  for(i in 1:nrow(query)){
    point_in_space <- pip3d(coords, faces, query[i,])
    if(point_in_space == -1){
      closest_point_of_convex <- cbind(coords, query[i,])
      closest_point_of_convex$distance <- apply(closest_point_of_convex, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
      distance <- min(closest_point_of_convex$distance)
    if(distance <= 20){
      error <- 1
    } 
    if(distance <= 50 && distance > 20){
      error <- 2
    } 
    if(distance > 50){
      error <- 3
    }
  }
  if(point_in_space != -1){
    error <- 0
  }
    
  df <- rbind(df, error)
  rownames(df)[i] <- paste(query[i,1], query[i,2], query[i,3], sep = ",")
  }
  colnames(df) <- "error"
  
  sum_error <- sum(df[,1])
  return(sum_error)
}

rgb2hex <- function(r,g,b) sprintf('#%s',paste(as.hexmode(c(r,g,b)),collapse = ''))
crgb <- t(col2rgb(cc <- colors()))

hex_colors <- c()
for(i in 1:nrow(crgb)){
  hex_one <- rgb2hex(crgb[i,1],crgb[i,2],crgb[i,3])
  hex_colors <- rbind(hex_colors, hex_one)
}

# nrow(hex_colors)
# 
# which(hex_colors=="#333")
# crgb[262,]
# 
# rgb2hex(200,15,15)

a <- ConvexHullCIELab(Lab[,1], Lab[,2], Lab[,3])
k <- as.matrix(LABdata)

RGB2Lab(crgb/100)


plot(as.matrix(RGBtoLabCoords), pch=19)
plot(x = 49.761335 , y = 5.330474, xlim = c(-100,100), ylim = c(-100,100))
lines(convex_coords, col="blue") # Cielab

plot(x = 0.5515779 , y = 0.3082405, xlim = c(-0,1), ylim = c(-0,1))
plot(as.matrix(ConvexCloud), pch=19, xlim = c(-100,100), ylim = c(-100,100))
plot(as.matrix(Cloud_scaled), pch=19, xlim = c(-100,100), ylim = c(-100,100))
lines(Cloud_scaled, col="red")
lines(ConvexCloud, col="red")

ConvexCloud <- rotate3d(obj = ConvexCloud, angle = 1, x = 1, y = 0, z = 0)
ConvexCloud <- rotate3d(obj = ConvexCloud, angle = 10, x = 0, y = 1, z = 0)
ConvexCloud <- rotate3d(obj = ConvexCloud, angle = 10, x = 0, y = 0, z = 1)
ConvexCloud <- Cloud_scaled




