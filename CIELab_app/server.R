#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
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

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
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
    
    myvals <- reactiveValues(
        umap_dist = NULL,
        start.values = NULL
    )
    
    datasetInput <- reactive({
        switch(input$dataset,
               # "GSE75748_time_course" = data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_time_course_ec.csv")),
               # "GSE75748_cell_type" = data.frame(read.csv("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_cell_type_ec.csv"))
               # "GSE52529" = data.frame(read.delim("/Users/tgn531/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE52529_fpkm_matrix.txt"))
               
               # "GSE75748_time_course" = data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_tc.tsv")),
               # "GSE75748_cell_type" = data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
             
               "GSE75748_time_course" = 0,
               "GSE75748_cell_type" = 1
             

               )
    })
    
    
    observeEvent(input$dataset, {
        if(input$dataset == "GSE75748_time_course"){

            #     rownames(matrix_counts) <- matrix_counts[,1]
            #     matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
            #     matrix_counts <- data.matrix(matrix_counts) # must be a matrix object!
            #
            #     TransposedMatrixCounts <- t(matrix_counts)
            #
            #     cnt <- 0
            #     FC_MatrixCounts <- c()
            #     for(i in 1:ncol(TransposedMatrixCounts)){
            #         mean_col <- mean(TransposedMatrixCounts[,i])
            #         FC <-  (TransposedMatrixCounts[,i]+1)/(mean_col+1)
            #         FC_MatrixCounts <- cbind(FC_MatrixCounts, FC)
            #
            #         cnt <- cnt + 1
            #         if(cnt %% 1000 == 0){
            #             print(cnt)
            #     }
            # }
            # colnames(FC_MatrixCounts) <- rownames(matrix_counts)
            # log2FC_MatrixCounts <- log2(FC_MatrixCounts)

            # 1
            # log2FC_MatrixCounts <- log2FC_MatrixCounts_tc[1:500,1:500]
            log2FC_MatrixCounts <- log2FC_MatrixCounts_tc
            # log2FC_MatrixCounts <- data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_tc.tsv"))
            # print(log2FC_MatrixCounts[1:10,1:10])
            print("1")
            SDMoreThanOnelog2FC_MatrixCounts <- c()
            for(i in 1:ncol(log2FC_MatrixCounts)){
                std <- sd(log2FC_MatrixCounts[,i])
                SDMoreThanOnelog2FC_MatrixCounts <- cbind(SDMoreThanOnelog2FC_MatrixCounts, std)
            }
            colnames(SDMoreThanOnelog2FC_MatrixCounts) <- colnames(log2FC_MatrixCounts)


            SDMoreThanOnelog2FC_MatrixCounts <- select(as.data.frame(log2FC_MatrixCounts),
                                                       -c(names(which(SDMoreThanOnelog2FC_MatrixCounts[1,]<1))))

            data <- CreateSeuratObject(counts = SDMoreThanOnelog2FC_MatrixCounts)

            all.genes <- rownames(data)
            data <- ScaleData(data, do.scale =  F, do.center = F, features = all.genes)
            data <- FindVariableFeatures(object = data, selection.method = 'mvp') #mvp because of error in log

            data <- RunPCA(data, npcs = 50, features = VariableFeatures(object = data))
            data <- FindNeighbors(data, dims = 1:15)
            data <- FindClusters(data, resolution = 0.5, algorithm= 1) # color in Seurat umap output

            data <- RunUMAP(data, dims = 1:50, n.components = 3L)
            # DimPlot(data, reduction = "umap")

            data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
            umap_dist <- data_umap_coord
            myvals$umap_dist <- umap_dist
            # print(myvals$umap_dist)
        }
    })
    observe({
        if(input$dataset == "GSE75748_cell_type"){
            # 2
            # log2FC_MatrixCounts <- log2FC_MatrixCounts_ct
            # log2FC_MatrixCounts_ct <- data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
            log2FC_MatrixCounts <- data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
            print("2")
            SDMoreThanOnelog2FC_MatrixCounts <- c()
            for(i in 1:ncol(log2FC_MatrixCounts)){
                std <- sd(log2FC_MatrixCounts[,i])
                SDMoreThanOnelog2FC_MatrixCounts <- cbind(SDMoreThanOnelog2FC_MatrixCounts, std)
            }
            colnames(SDMoreThanOnelog2FC_MatrixCounts) <- colnames(log2FC_MatrixCounts)


            SDMoreThanOnelog2FC_MatrixCounts <- select(as.data.frame(log2FC_MatrixCounts),
                                                       -c(names(which(SDMoreThanOnelog2FC_MatrixCounts[1,]<1))))

            data <- CreateSeuratObject(counts = SDMoreThanOnelog2FC_MatrixCounts)

            all.genes <- rownames(data)
            data <- ScaleData(data, do.scale =  F, do.center = F, features = all.genes)
            data <- FindVariableFeatures(object = data, selection.method = 'mvp') #mvp because of error in log

            data <- RunPCA(data, npcs = 50, features = VariableFeatures(object = data))
            data <- FindNeighbors(data, dims = 1:15)
            data <- FindClusters(data, resolution = 0.5, algorithm= 1) # color in Seurat umap output

            data <- RunUMAP(data, dims = 1:50, n.components = 3L)
            # DimPlot(data, reduction = "umap")

            data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
            umap_dist <- data_umap_coord
            
            myvals$umap_dist <- umap_dist
        }
    })
    
    
    observe({
        umap_dist <- myvals$umap_dist 
        WL <- input$weightL
        Wa <- input$weightA
        Wb <- input$weightB
        
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
                ConvexCloud[,1] <- WL*ConvexCloud[,1]
                ConvexCloud[,2] <- Wa*ConvexCloud[,2]
                ConvexCloud[,3] <- Wb*ConvexCloud[,3]
            }
            return(ConvexCloud)
        } # Transform the UMAP Convex to otpimize it
        
        dat <- UMAPConvex(umap_dist)
        polygon <- ColorSpacePolygon(RGB_space)
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
        
        myvals$start.values <- result$par
    })
    
    
    
    output$plotly_plot <- renderPlotly({
        start.values <- myvals$start.values
        umap_dist <- myvals$umap_dist 
        
    NewUMAP <- Scaling(umap_dist, start.values[1]*input$scaling)
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
    
    fig
    })
    
    
    
    output$satellite1 <- renderPlot({
        umap_dist <- myvals$umap_dist 
        start.values <- myvals$start.values
        ConvexCloud <- UMAPConvex(umap_dist)
        ConvexCloud <- Scaling(ConvexCloud, start.values[1]*input$scaling)
        ConvexCloud <- Rotation(ConvexCloud,  start.values[2]  , start.values[3] , start.values[4])
        ConvexCloud <- Translation(ConvexCloud, start.values[5] , start.values[6] ,start.values[7])

        plot(as.matrix(ConvexCloud[,c(1,2)]), pch=19, xlim = c(-100,100), ylim = c(-100,100))
        lines(polygon[, c(1,2)], col="blue") # Cielab

    })
    
    output$satellite2 <- renderPlot({
        umap_dist <- myvals$umap_dist 
        start.values <- myvals$start.values
        ConvexCloud <- UMAPConvex(umap_dist)
        ConvexCloud <- Scaling(ConvexCloud, start.values[1]*input$scaling)
        ConvexCloud <- Rotation(ConvexCloud,  start.values[2]  , start.values[3] , start.values[4])
        ConvexCloud <- Translation(ConvexCloud, start.values[5] , start.values[6] ,start.values[7])

        plot(as.matrix(ConvexCloud[,c(1,3)]), pch=19, xlim = c(-100,100), ylim = c(-100,100))
        lines(polygon[, c(1,3)], col="blue") # Cielab

    })


    # ConvexCloud <- UMAPConvex(umap_dist)
    # ConvexCloud <- Scaling(ConvexCloud, start.values[1])
    # ConvexCloud <- Rotation(ConvexCloud,  start.values[2]  , start.values[3] , start.values[4])
    # ConvexCloud <- Translation(ConvexCloud, start.values[5] , start.values[6] ,start.values[7])
    # 
    # plot(as.matrix(ConvexCloud[,c(1,2)]), pch=19, xlim = c(-100,100), ylim = c(-100,100))
    # lines(polygon[, c(1,2)], col="blue") # Cielab
    # 
    # colnames(ConvexCloud) <- colnames(RGBtoLabCoords)
    # a <- rbind(as.data.frame(RGBtoLabCoords), as.data.frame(ConvexCloud))
    # a$colors <- as.factor(c(rep(1, nrow(RGBtoLabCoords)), rep(0, nrow(ConvexCloud))))
    # 
    # library(plotly)
    # fig <- plot_ly(a, x = a[,1], y = a[,2], z = a[,3], color = a[,4], colors = c('#BF382A', '#0C4B8E'), mode='lines+markers',
    #                line = list(width = 6), marker = list(size = 3.5))
    # # fig <- fig %>% add_markers() # instead of lines
    # fig <- fig %>% layout(scene = list(xaxis = list(title = 'L'),
    #                                    yaxis = list(title = 'a'),
    #                                    zaxis = list(title = 'b')))
    # 
    # fig

})
