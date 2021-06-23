#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    source('~/Documents/Documents – SUN1012692/GitHub/CIELAB/global.R', local = TRUE)
    polygon <- ColorSpacePolygon(RGB_space)
    
    myvals <- reactiveValues(
        umap_dist = NULL,
        start.values = NULL,
        optional.start.values = NULL,
        simplex_vectors = NULL
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
            log2FC_MatrixCounts <- log2FC_MatrixCounts_tc[1:500,1:500]
            # log2FC_MatrixCounts <- log2FC_MatrixCounts_tc
            # log2FC_MatrixCounts_tc <- data.frame(read.delim("~/Documents/Documents – SUN1012692/GitHub/CIELAB/log2FC_MatrixCounts_tc.tsv"))
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
            # print(SDMoreThanOnelog2FC_MatrixCounts)
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
        
        if(input$dataset == "GSE75748_cell_type"){
            # 2
            log2FC_MatrixCounts <- log2FC_MatrixCounts_ct
            # log2FC_MatrixCounts <- log2FC_MatrixCounts_ct[1:1000,1:1000]
            # log2FC_MatrixCounts_ct <- data.frame(read.delim("~/Documents/Documents – SUN1012692/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
            # log2FC_MatrixCounts <- data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
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
        input$weightButton
        
        WL <- isolate(input$weightL)
        Wa <- isolate(input$weightA)
        Wb <- isolate(input$weightB)
        
        #--- Functions ---#
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
        #---#
        
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
        angle <- 1
        start.values <- c(S , pi/4 ,pi/4, pi/4, TrL, Tra , Trb) # S, RotL, Rota, Rotb,  TrL, Tra, Trb
        simplex_vectors <- c()
        k <- c()
        
        
        # cores=8
        # cores=detectCores()-1
        # # library(doFuture)
        # registerDoFuture()
        # plan(multiprocess)
        # 
        # # cl <- makeCluster(cores[1]-1) #not to overload your computer
        # # registerDoParallel(cl)
        # stopCluster(cl)
        # # library(future.apply)
        # # plan(multiprocess) ## => parallelize on your local computer
        # # 
        # start_time <- Sys.time()
        # print("start_time")
        # print(start_time)
        # 
        # Simplex_optim <- foreach(x = cores) %dopar% {
        #   for(i in 1:25){
        #     tmp <- optim(par = start.values,
        #                method = "Nelder-Mead",
        #                MasterFunction,
        #                data = dat,
        #                polygon = polygon,
        #                faces = faces,
        #                control=list(fnscale=-1, maxit=1000)) #, trace=T
        # 
        #   k[[i]] <- tmp
        #   angle <- angle + 0.5
        #   start.values <- c(start.values[1], start.values[2]+(pi/angle), start.values[3], start.values[4], start.values[5], start.values[6], start.values[7]) # mirror rotation in x
        # 
        #   simplex_vectors <- rbind(simplex_vectors, k[[i]]$par)
        # 
        #   }
        #   return(simplex_vectors)
        # }
        # 
        # Simplex_optim <- Simplex_optim[[1]]
        # print(Simplex_optim[1])
        # end_time <- Sys.time()
        # print(end_time)
        # myvals$start.values <- Simplex_optim[which(Simplex_optim[,1] == max(Simplex_optim[,1]))[1], ]
        # myvals$simplex_vectors <- Simplex_optim
         
        
        # start_time <- Sys.time()
        # print("start_time")
        # print(start_time)

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
        }
        # print(simplex_vectors)
        myvals$start.values <- simplex_vectors[which(simplex_vectors[,1] == max(simplex_vectors[,1]))[1], ]
        myvals$simplex_vectors <- simplex_vectors

        # end_time <- Sys.time()
        # print(end_time)
      
        print("~End~")  

    })
    
    
    output$table <- renderDataTable({
      req(input$weightButton)
      start.values <- myvals$start.values
      simplex_vectors <- myvals$simplex_vectors

      optional_values <- c()
      for(i in 1:nrow(simplex_vectors)){
        vector_dist <- distmat(start.values[1], simplex_vectors[i,1])
        if(vector_dist >= 5) optional_values <- rbind(optional_values, simplex_vectors[i,])
      }
      
      
      if(length(optional_values)!=0){
        vector_final_ordered <- as.data.frame(optional_values)
        if(nrow(vector_final_ordered)==1){
          vector_final_ordered <- as.data.frame(optional_values)
        } else vector_final_ordered <- as.data.frame(optional_values[order(optional_values[,1], decreasing = T),])
        
        vector_final_ordered <- round(vector_final_ordered,2)
        datatable(vector_final_ordered, selection = c("single"), colnames = c("Scaling factor", "Rotation L*", "Rotation a*", "Rotation b*", "Translation L*", "Translation a*", "Translation b*"))
      }
      
      
        })

    output$list_of_parameters <- renderPrint({
      req(input$table_rows_selected)
      start.values <- myvals$start.values
      simplex_vectors <- myvals$simplex_vectors
      
      s <- input$table_rows_selected
      optional_values <- c()
      for(i in 1:nrow(simplex_vectors)){
        vector_dist <- distmat(start.values[1], simplex_vectors[i,1])
        if(vector_dist >= 5) optional_values <- rbind(optional_values, simplex_vectors[i,])
      }
      
      if(length(optional_values)!=0){
        vector_final_ordered <- as.data.frame(optional_values)
        if(nrow(vector_final_ordered)==1){
          vector_final_ordered <- as.data.frame(optional_values)
        } else vector_final_ordered <- as.data.frame(optional_values[order(optional_values[,1], decreasing = T),])
        vector_final_ordered <- round(vector_final_ordered,2)  
        myvals$optional.start.values <- vector_final_ordered[s,]
      }
      
    })
    

    output$plotly_plot <- renderPlotly({
        input$scalingButton
      
      if(!is.null(input$table_rows_selected)){
        start.values <- as.numeric(myvals$optional.start.values)
      } else {
        start.values <- myvals$start.values
      }
      
        umap_dist <- myvals$umap_dist 
        
    NewUMAP <- Scaling(umap_dist, start.values[1]*isolate(input$scaling))
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
                title = "UMAP1")
    
    axy <- list(nticks = 4,
                title = "UMAP2")
    
    axz <- list(nticks = 4,
                title = "UMAP3")
    
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
        input$scalingButton
        umap_dist <- myvals$umap_dist 
        
        if(!is.null(input$table_rows_selected)){
          start.values <- as.numeric(myvals$optional.start.values)
        } else {
          start.values <- myvals$start.values
        }
        
        ConvexCloud <- UMAPConvex(umap_dist)
        ConvexCloud <- Scaling(ConvexCloud, start.values[1]*isolate(input$scaling))
        ConvexCloud <- Rotation(ConvexCloud,  start.values[2]  , start.values[3] , start.values[4])
        ConvexCloud <- Translation(ConvexCloud, start.values[5] , start.values[6] ,start.values[7])
        plot(as.matrix(ConvexCloud[,c(1,2)]), pch=19, xlim = c(-100,100), ylim = c(-100,100))
        lines(polygon[, c(1,2)], col="blue") # Cielab


    })
    
    output$satellite2 <- renderPlot({
        input$scalingButton
        umap_dist <- myvals$umap_dist 
        
        if(!is.null(input$table_rows_selected)){
          start.values <- as.numeric(myvals$optional.start.values)
        } else {
          start.values <- myvals$start.values
        }
        
        ConvexCloud <- UMAPConvex(umap_dist)
        ConvexCloud <- Scaling(ConvexCloud, start.values[1]*isolate(input$scaling))
        ConvexCloud <- Rotation(ConvexCloud,  start.values[2]  , start.values[3] , start.values[4])
        ConvexCloud <- Translation(ConvexCloud, start.values[5] , start.values[6] ,start.values[7])
        
        
        
        
        plot(as.matrix(ConvexCloud[,c(1,3)]), pch=19, xlim = c(-100,100), ylim = c(-100,100))
        lines(polygon[, c(1,3)], col="blue") # Cielab

    })


})
