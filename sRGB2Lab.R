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
Lab2RGB <- function(mat) {
  # cast mat in case of a single numeric vector
  mat <- matrix(mat, ncol=3)
  # same spec as RGBtoLab
  # Lab -> XYZ
  thres1 <- 0.008856
  thres2 <- 0.0031308
  
  y <- (mat[,1] + 16) / 116
  x <- mat[,2] / 500 + y
  z <- y - mat[,3] / 200
  xyz <- cbind(x,y,z)
  xyzthres <- xyz > thres1
  xyz <- xyzthres * (xyz^3) + (!xyzthres) * ((xyz - 16 / 116) / 7.787)
  xyz <- sweep(xyz, 2, c(0.950456,1,1.088754), "*")
  
  # XYZ -> RGB
  Minv = c(3.240479,-1.537150,-0.498535,
           -0.969256, 1.875992, 0.041556,
           0.055648,-0.204043, 1.057311)
  Minv <- matrix(Minv, nrow=3, byrow=TRUE)
  RGB <- xyz %*% t(Minv)
  RGBthres <- RGB > thres2
  # manage NaN risk with x^a and a<1
  RGB[RGBthres] <- 1.055 * RGB[RGBthres]^(1/2.4) - 0.055
  RGB[!RGBthres] <- 12.92 * RGB[!RGBthres]
  # bound in case of very small absolute <0 values, and resp. for >1 values.
  RGB[RGB<0] <- 0
  RGB[RGB>1] <- 1
  return(RGB)
}

# RGB space
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

library(colorspace)
# x <- polarLUV(L = 70, C = 50, H = c(0, 120, 240))
# y <- as(x, "sRGB")
# 
# 
x <- LABdata[1:10,]
y <- as(x, "sRGB")
# 
# RGB2Lab(c(0,0,255))
# lrgb <- Lab2RGB(Lab)
# 


library(plotly)

# RGB gamut - cube
fig <- plot_ly(type = 'mesh3d',
               x = c(0,   0,    255,    255,    0,     0,     255, 255),
               y = c(0,   255,  255,    0,      0,     255,   255,   0),
               z = c(0,   0,    0,      0,      255,   255,   255, 255),
               
               # x = c(0, 0, 1, 1, 0, 0, 1, 1),
               # y = c(0, 1, 1, 0, 0, 1, 1, 0),
               # z = c(0, 0, 0, 0, 1, 1, 1, 1),
               
               # Cube
               i = c(7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2),
               j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3),
               k = c(0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6),
               # intensity = seq(0, 1, length = 8),
               # color = seq(0, 1, length = 8)
               )
fig

# a == matrix transformed into CIELAB 

axx <- list(
  nticks = 4,
  title = "a",
  range = c(-128,128)
)

axy <- list(
  nticks = 4,
  title = "b",
  range = c(-128,128)
)

axz <- list(
  nticks = 4,
  title = "L",
  range = c(0,100)
)

fig <- plot_ly(type = 'mesh3d',
               x = RGBtoLabCoords[,3],
               y = RGBtoLabCoords[,2],
               z = RGBtoLabCoords[,1],
              
               # Cube
               # i = c(7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2),
               # j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3),
               # k = c(0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6),
               # intensity = seq(0, 1, length = 8),
               # color = seq(0, 1, length = 8),
               # colors = c(colors())
               # colors = c(hex(LABdata, fix = TRUE)[100:767])
               # colorscale='Rainbow'
)
fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig

# sRGB color colors
fig <- plot_ly(type = 'mesh3d',
               x = crgb[,3],
               y = crgb[,2],
               z = crgb[,1],
               
               # Cube
               i = c(7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2),
               j = c(3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3),
               k = c(0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6),
               # intensity = seq(0, 1, length = 8),
               # color = seq(0, 1, length = 8),
               # colors = c(colors())
               # colors = c(hex(LABdata, fix = TRUE)[100:767])
               # colorscale='Rainbow'
)
# install.packages("scatterplot3d")
library(scatterplot3d)
# a <- data.frame("L" = rep(80, nrow(a)), "a"= a[,1], "b"= a[,2])
scatterplot3d(RGBtoLabCoords, highlight.3d = F, col.axis = "black", angle = 290, box = F,
              pch = 20)


# rgb <- data.frame("R"=r, "G"=g, "B"=b)

# Normalize to 0-1 --> R/255, G/255, B/255
mat <- as.matrix(data.frame(R = c(0, 0, 1, 1, 0, 0, 1, 1),
                            G = c(0, 1, 1, 0, 0, 1, 1, 0),
                            B = c(0, 0, 0, 0, 1, 1, 1, 1)))

RGB_space <- as.matrix(RGB_space)
RGBtoLabCoords <- as.data.frame(RGB2Lab(RGB_space/255)) # RGB space

a <- as.data.frame(RGB2Lab(mat)) # RGB corners
a <- as.data.frame(RGB2Lab(crgb/255)) # All RGB colors in R
a <- as.data.frame(RGB2Lab(Rrb/255)) # Rainbow colors



