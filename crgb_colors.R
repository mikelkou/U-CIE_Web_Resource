cubedraw <- function(res3d, min = 0, max = 255, cex = 2)
{
  cube01 <- rbind(0,c(1,0,0),c(1,1,0),1,c(0,1,1),c(0,0,1),c(1,0,1),
                  c(1,0,0),c(1,0,1),1,c(1,1,0),
                  c(0,1,0),c(0,1,1), c(0,1,0),0)
  cub <- min + (max-min)* cube01
  res3d$points3d(cub[ 1:11,], cex = cex, type = 'b', lty = 1)
  res3d$points3d(cub[11:15,], cex = cex, type = 'b', lty = 3)
}
crgb <- t(col2rgb(cc <- colors()))
rr <- scatterplot3d(crgb, color = cc, box = FALSE, angle = 24)
cubedraw(rr)
Rrb <- t(col2rgb(rbc <- rainbow(201)))
rR <- scatterplot3d(Rrb, color = rbc, box = FALSE, angle = 24)
cubedraw(rR)
rR$points3d(Rrb, col = rbc, pch = 16)
