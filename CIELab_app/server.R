# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "output")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2)
  # path_conf <- config$global
  source('global.R', local = TRUE)
  
  polygon <- ColorSpacePolygon(RGB_space)
    
    myvals <- reactiveValues(
        umap_dist = NULL,
        start.values = NULL,
        optional.start.values = NULL,
        simplex_vectors = NULL,
        uploaded_df = NULL,
        NewUMAP = NULL,
        RemoveGenesFromConvexCloud = NULL,
        select_data = NULL,
        deselect = NULL,
        NewUMAP_stored = NULL,
        umap_dist_stored = NULL,
        DownloadFile = NULL,
        RemovedGenes = NULL,
        select_data_colors = NULL,
        colors = NULL
    )
    
    read_data <- function(datapath, type = c("txt"), header = T, sep = "\t", quote = "\"") ({
      if(endsWith(input$file1$name, 'xlsx')){
        dataset1 <- as.data.frame(read_excel(input$file1$datapath, 1, col_names = ifelse(input$header==T, T, F)))
      } 
      if(endsWith(input$file1$name, 'csv')){
        dataset1 <- read.table(input$file1$datapath, header = input$header, sep = ",", quote = input$quote)
      } 
      if(endsWith(input$file1$name, 'txt')){
        dataset1 <- read.table(input$file1$datapath, header = input$header, sep = "\t", quote = input$quote)
      }
      if(endsWith(input$file1$name, 'tsv')){
        dataset1 <- read.table(input$file1$datapath, header = input$header, sep = "\t", quote = input$quote)
      }
      
      if(ncol(dataset1)<3){
        showModal(modalDialog(title = "Not even 3D? Please check again the dataset.", easyClose = T, fade = T))
        dataset1 <- NULL
      }
      
      return(dataset1)
    })
    
    loadNetworkFromFile <- function() {
      dataset1 <- NULL
      set.seed(123)
      
      switch(
        input$uiLoadGraphOptionsInput,
        oF = {
          if (!is.null(input$file1)) {
            if(input$file1$type == "application/pdf"){
              # "This is PDF! Please upload tsv/csv files."
              showModal(modalDialog(title = "This is PDF! Please upload tsv/csv files.", easyClose = T, fade = T))
            } else {
            dataset1 <- read_data(input$file1$datapath)
            }
          }
        },
        oR_Example1 = {
          dataset1 <- data.frame(read.delim("log2FC_MatrixCounts_tc_cropped.tsv"))
          # log2FC_MatrixCounts <- log2FC_MatrixCounts_tc[1:500,1:500]
          # log2FC_MatrixCounts <- log2FC_MatrixCounts_tc
          # log2FC_MatrixCounts_tc <- data.frame(read.delim("~/Documents/Documents – SUN1012692/GitHub/CIELAB/log2FC_MatrixCounts_tc.tsv"))
          # log2FC_MatrixCounts <- data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_tc.tsv"))
          # print(log2FC_MatrixCounts[1:10,1:10])
        },
        oR_Example2 = {
          dataset1 <- data.frame(read.delim("log2FC_MatrixCounts_ct.tsv"))
          # dataset1 <- log2FC_MatrixCounts_ct[1:500,1:500]
          # log2FC_MatrixCounts <- log2FC_MatrixCounts_ct
          # log2FC_MatrixCounts <- log2FC_MatrixCounts_ct[1:1000,1:1000]
          # log2FC_MatrixCounts_ct <- data.frame(read.delim("~/Documents/Documents – SUN1012692/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
          # log2FC_MatrixCounts <- data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
        }
      )
      return(dataset1)
    }
    
    # write.table(log2FC_MatrixCounts_tc[1:500,1:500], "log2FC_MatrixCounts_tc_cropped.tsv",
    #             quote = F,sep = "\t", row.names = T, col.names = T
    #             )
    output$uiLoadGraphOptionsOutput <- renderUI({
      if (is.null(input$uiLoadGraphOptionsInput))
        return()
      
      # Depending on input$input_type, we'll generate a different UI
      # component and send it to the client.
      if (input$uiLoadGraphOptionsInput == "oF") {
        wellPanel(fileInput(
          "file1",
          "Choose file to upload",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv",
            ".xlsx"
          ))
        )
      } else {}
      
    })
    
    ###### Upload #######
    doAddNetwork <- observeEvent(input$btnanalysis, {
      myvals$uploaded_df <- loadNetworkFromFile()
      # print(head(myvals$uploaded_df))
      # if (!is.null(dataset)) {
      #     nid <- UUIDgenerate(T)      #time-base UUID is generated
      #     nn <- input$networkName
      #     cnt <- 1                    #count
      #     while (nn %in% reactiveVars$StoredNetworks$name) {
      #       #reactiveVars: represents a single reactive variable.
      #       cnt <- cnt + 1
      #       nn <-
      #         paste(input$networkName, cnt)      #paste: converts its arguments (via as.character) to character strings
      #     }
      #     df <- data.frame(id = nid, name = nn, stringsAsFactors = F)
      # }
        
        
      #   if (nrow(dataset) >= 10000) {
      #     createAlert(
      #       session,
      #       "tabUpload_up_to_10000_rows",
      #       "fileUploadAlert_up_to_10000_rows",
      #       title = "Warning !",
      #       style = "danger",
      #       content = paste0(
      #         "Please make sure that your network has less than 10,000 connections."
      #       ),
      #       append = FALSE
      #     )
      #     dataset <- NULL
      #   }
      #   
      # } else
      #   createAlert(
      #     session,
      #     "tabUploadSideAlert",
      #     "fileUploadAlert",
      #     title = "ERROR !",
      #     style = "danger",
      #     content = paste0(
      #       "An error occurred while trying to read your file. Please make sure that it is formatted according to the requirements."
      #     ),
      #     append = FALSE
      #   )
      
    })
    
    
    
    output$contents <- renderDataTable({
      uploaded_df <- loadNetworkFromFile()
      tryCatch(
        {
          df <- uploaded_df
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      if(input$disp == "head") {
        datatable(df[1:10,], options = list(scrollX = TRUE))
      }
      else {
        datatable(df, options = list(scrollX = TRUE))
      }
      
    })
    
    
    observe({
      req(input$btnanalysis)
        if (!is.null(myvals$uploaded_df) & !is.null(input$matrix)) {
          print("ho")
        if(input$matrix == 'Single-cells'){
          print("sc")
          # matrix_counts <- read.csv("~/Desktop/CBPP_22012021/Lars_Lab/single_cells/GSE75748_sc_cell_type_ec.csv")
          matrix_counts <- myvals$uploaded_df
          
          rownames(matrix_counts) <- matrix_counts[,1]
              matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
              matrix_counts <- data.matrix(matrix_counts) # must be a matrix object!

              TransposedMatrixCounts <- t(matrix_counts)
              
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
          # 
          
          # SDMoreThanOnelog2FC_MatrixCounts <- c()
          # for(i in 1:ncol(log2FC_MatrixCounts)){
          #   std <- sd(log2FC_MatrixCounts[,i])
          #   SDMoreThanOnelog2FC_MatrixCounts <- cbind(SDMoreThanOnelog2FC_MatrixCounts, std)
          # }
          # colnames(SDMoreThanOnelog2FC_MatrixCounts) <- colnames(log2FC_MatrixCounts)
          # 
          # SDMoreThanOnelog2FC_MatrixCounts <- select(as.data.frame(log2FC_MatrixCounts),
          #                                            -c(names(which(SDMoreThanOnelog2FC_MatrixCounts[1,]<1))))
         
         # data <- CreateSeuratObject(counts = SDMoreThanOnelog2FC_MatrixCounts)
              withConsoleRedirect("console", {
              data <- CreateSeuratObject(counts = TransposedMatrixCounts)
         all.genes <- rownames(data)
         data <- ScaleData(data, do.scale =  F, do.center = F, features = all.genes)
         
         data <- NormalizeData(data, normalization.method = "LogNormalize")
         
         data <- FindVariableFeatures(object = data) #, selection.method = 'mvp' because of error in log

         data <- RunPCA(data, npcs = 50, features = VariableFeatures(object = data))
         data <- FindNeighbors(data, dims = 1:length(data@reductions$pca))
         data <- FindClusters(data, resolution = 0.5, algorithm= 1) # color in Seurat umap output
         data <- RunUMAP(data, dims = 1:50, n.components = 3L)
         
         data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
              })
         umap_dist <- data_umap_coord
         myvals$umap_dist <- umap_dist
         # myvals$umap_dist <- myvals$uploaded_df
         
         print("Single cells done!")
        } 
      
      if(input$matrix == 'High Dimensional'){
        # data <- CreateSeuratObject(counts = myvals$uploaded_df)
        # all.genes <- rownames(data)
        # data <- ScaleData(data, do.scale =  F, do.center = F, features = all.genes)
        # data <- FindVariableFeatures(object = data, selection.method = 'mvp') #mvp because of error in log
        # data <- RunPCA(data, npcs = 50, features = VariableFeatures(object = data))
        # data <- FindNeighbors(data, dims = 1:length(data@reductions$pca))
        # data <- FindClusters(data, resolution = 0.5, algorithm= 1) # color in Seurat umap output
        # data <- RunUMAP(data, dims = 1:50, n.components = 3L)
        # data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
        withConsoleRedirect("console", {
        data <- umap(myvals$uploaded_df, ret_nn = TRUE, n_neighbors = 5, n_components = 3) # library(uwot)
        data_umap_coord <- as.data.frame(data$embedding)
        })
        
        umap_dist <- data_umap_coord
        myvals$umap_dist <- umap_dist
      }
      if(input$matrix == 'Distance matrix'){
        withConsoleRedirect("console", {
        data <- as.matrix(dist(myvals$uploaded_df))
        data <- umap(data, input="dist", n_components = 3)
        data_umap_coord <- as.data.frame(data$layout)
        })
        umap_dist <- data_umap_coord
        myvals$umap_dist <- umap_dist
        print("Distance")
      }
      if(input$matrix == "3D data"){
        myvals$umap_dist <- myvals$uploaded_df
        print("3D")
      }
      }
      # }
        if(is.null(input$matrix)){
          print("hi")
      if(input$uiLoadGraphOptionsInput == "oR_Example1"){
        # print(input$btnanalysis)
        input$btnanalysis
        
          withConsoleRedirect("console", {
            
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
            log2FC_MatrixCounts <- loadNetworkFromFile()
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

            data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
          })
            umap_dist <- data_umap_coord
            myvals$umap_dist <- umap_dist
            # print(myvals$umap_dist)
        }
        
        if(input$uiLoadGraphOptionsInput == "oR_Example2"){
            # 2
          input$btnanalysis
          withConsoleRedirect("console", {
          log2FC_MatrixCounts <- loadNetworkFromFile()
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
            
            data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
          })
            umap_dist <- data_umap_coord
            
            myvals$umap_dist <- umap_dist
        }
      # myvals$umap_dist_stored <- myvals$umap_dist
        
      }
    })
    
    
    observe({
      withProgress(min = 0, max = 1, {
        incProgress(message = "Calculation in progress",
                    # detail = "This may take a while...",
                    amount = .1)
        
        umap_dist <- myvals$umap_dist 
        
        if(is.null(umap_dist)){}
        else{
          if(!is.null(myvals$RemoveGenesFromConvexCloud)){
            req(input$remove_genes)
            umap_dist <- myvals$RemoveGenesFromConvexCloud
          }
          
          
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
        }
      })
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
      umap_dist <- myvals$umap_dist
      if(is.null(umap_dist)){}
      else{
        if(!is.null(myvals$RemoveGenesFromConvexCloud)){
          req(input$remove_genes)
          umap_dist <- myvals$RemoveGenesFromConvexCloud
        }
        
        
      input$scalingButton
      
      if(!is.null(input$table_rows_selected)){
        start.values <- as.numeric(myvals$optional.start.values)
      } else {
        start.values <- myvals$start.values
      }
     
      
    NewUMAP <- Scaling(umap_dist, start.values[1]*isolate(input$scaling))
    NewUMAP <- Rotation(as.matrix(NewUMAP), start.values[2]  , start.values[3] , start.values[4])
    NewUMAP <- Translation(NewUMAP, start.values[5] , start.values[6] ,start.values[7])
    # myvals$NewUMAP <- NewUMAP
    # NewUMAP <- myvals$NewUMAP
    
    if(input$remove_genes){
      print("remove genes")
      remove_genes()
    }
    if(input$reset_genes){
      print("reset")
      if(input$remove_genes[1]<=input$reset_genes[1]){
        reset()
        # NewUMAP <- myvals$NewUMAP_stored
      }
    } else {
      myvals$NewUMAP <- NewUMAP
      
    }
    NewUMAP <- myvals$NewUMAP
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
        row.names = c(rownames(NewUMAP)),
        class = "data.frame"
    )
    
    library(colorspace)
    LABdata <- with(rawdata, LAB(Lstar, Astar, Bstar))
    
    #--- Download file ---#
    DownloadFile <- cbind(rownames(umap_dist), hex(LABdata, fix = TRUE))
    if(!is.null(myvals$RemovedGenes)){
      DownloadFile <- rbind(DownloadFile, myvals$RemovedGenes)
    }
    myvals$DownloadFile <- DownloadFile
    
    axx <- list(nticks = 4,
                title = "UMAP1")
    
    axy <- list(nticks = 4,
                title = "UMAP2")
    
    axz <- list(nticks = 4,
                title = "UMAP3")
    
    
    myvals$colors <- as.data.frame(cbind(NewUMAP,hex(LABdata, fix = TRUE)))
    
    fig <- plot_ly(
        data = as.data.frame(NewUMAP),
        x = NewUMAP[, 1],
        y = NewUMAP[, 2],
        z = NewUMAP[, 3],
        source = "A",
        type = "scatter3d",
        mode = "markers",
        text = c(rownames(NewUMAP)),
        # hoverinfo = 'text',
        marker = list(
            color = hex(LABdata, fix = TRUE),
            size = 3,
            width = 2
        )
        ,
        width = 735, height = 720
    )
    fig <- fig %>% layout(scene = list(
        xaxis = axx,
        yaxis = axy,
        zaxis = axz
    )) 
    fig
      } # else
    })
    
    remove_genes <- eventReactive(input$remove_genes, {
      
      if (length(myvals$select_data)==0) {
        return()
      } else {
        ConvexCloud <- myvals$NewUMAP
        ConvexCloud <- as.data.frame(ConvexCloud)
        myvals$NewUMAP_stored <- ConvexCloud
      
        if(myvals$select_data[,1]==2){
          comparison <- row.match(ConvexCloud[,c(1,3)],myvals$select_data[,3:4])
      } else comparison <- row.match(ConvexCloud[,c(1,2)],myvals$select_data[,3:4])
      
        myvals$RemoveGenesFromConvexCloud <- ConvexCloud[-c(which(!is.na(comparison))), ] # Not showing these points
        
        RemovedGenes <- ConvexCloud[c(which(!is.na(comparison))), ] # Color them gray in the output file
        RemovedGenes <- cbind(rownames(RemovedGenes), rep("#808080", nrow(RemovedGenes)))
        myvals$RemovedGenes <- RemovedGenes
        
        myvals$NewUMAP <- myvals$RemoveGenesFromConvexCloud
        myvals$umap_dist <- myvals$RemoveGenesFromConvexCloud
      }
    })
    
    
    reset <- eventReactive(input$reset_genes,{
      # print(input$reset_genes)
      myvals$NewUMAP <- myvals$NewUMAP_stored
    })
    
    
    
    output$satellite <- renderPlotly({
      # myvals$select_data <- event_data("plotly_selected")
      
      myvals$select_data <- event_data("plotly_selected", source = "A")

      if(input$remove_genes){
        print("remove genes")
        remove_genes()
      }
      if(input$reset_genes){
        print("reset")
        if(input$remove_genes[1]<=input$reset_genes[1]){
          reset()
          # NewUMAP <- myvals$NewUMAP_stored
        }
      } else {
        # myvals$NewUMAP <- NewUMAP
        umap_dist <- myvals$NewUMAP
      }
        
      umap_dist <- myvals$NewUMAP
        if(is.null(umap_dist)){}
        else{
          input$scalingButton
          
        if(!is.null(input$table_rows_selected)){
          start.values <- as.numeric(myvals$optional.start.values)
        } else {
          start.values <- myvals$start.values
        }
          
          
          ConvexCloud <- myvals$NewUMAP
          

        colordata = structure(
          list(
            Lstar = c(polygon[, 1]),
            Astar = c(polygon[, 2]),
            Bstar = c(polygon[, 3])
          ),
          .Names = c("Lstar", "Astar", "Bstar"),
          row.names = c(rownames(umap_dist)[1:nrow(polygon)]),
          class = "data.frame"
        )

        LABdata_polygon <- with(colordata, LAB(Lstar, Astar, Bstar))
        polygon_colors <- as.data.frame(cbind(polygon, hex(LABdata_polygon, fix = TRUE)))
        
        # Show colors of selected genes
        ConvexCloud <- myvals$NewUMAP
        ConvexCloud <- as.data.frame(ConvexCloud)
        
        if(length(myvals$select_data)==0) {
          myvals$select_data_colors <- as.data.frame(cbind(myvals$colors, rep("gray50", nrow(myvals$colors))))
          }
        if(!(length(myvals$select_data)==0)){
          if(myvals$select_data[,1]==2){
            comparison <- row.match(myvals$colors[,c(1,3)],myvals$select_data[,3:4])
          } else comparison <- row.match(myvals$colors[,c(1,2)],myvals$select_data[,3:4])
          
          a <- as.data.frame(cbind(myvals$colors, as.data.frame(comparison)))
          # print(a)
          a[,4] <- ifelse(is.na(a[,5]),  "gray50",  a[,4])
          myvals$select_data_colors <- a
          # print(a)
        }
        
        
        
        fig1 <- plot_ly(
          data = as.data.frame(ConvexCloud),
          x = ConvexCloud[, 1],
          y = ConvexCloud[, 2],
          source = "A",
          type = "scatter",
          mode = "markers",
          # symbol="square",
          legendgroup = 'UMAP point cloud', showlegend = T,
          name = 'UMAP point cloud',
          text = c(rownames(ConvexCloud)),
          hoverinfo = 'text',
          # marker = list(
          #   size = 5,
          #   width = 2,
          #   line = list(
          #     # color = 'gray50',
          #     color = myvals$select_data_colors[,4],
          #     width = 0.5)
          # )
          marker = list(color = I(myvals$select_data_colors[,4]), opacity = 1, size = 8)
          ,
          width = 735, height = 720
        )
        
        # fig1 <- fig1 %>% add_polygons(polygon, x = polygon[, 1], y = polygon[, 3], color=I("red"), opacity=0.2)
        fig1 <- fig1 %>% add_polygons(polygon, x = polygon[, 1], y = polygon[, 2], 
                                      fillcolor = 'rgba(206, 211, 214, 0.1)', 
                                      name = 'CIE L* a* b*',
                                      legendgroup = 'CIE L* a* b*', showlegend = T,
                                      marker = list(color = I(polygon_colors[,4]), opacity = 1, size = 10, symbol= 'cross', line = list(
                                        color = 'black', width=0.5)),
                                      text = c(polygon_colors[,4]), hoverinfo = 'text', 
                                      line = list(color = '#a5c0cf'))
        fig1 <- fig1 %>% 
          layout(dragmode = "select")
        
        fig1 <- fig1 %>% layout(
          xaxis = list(title = "L"),
          yaxis = list(title = "a"))
        
        fig2 <- plot_ly(
          data = as.data.frame(ConvexCloud),
          x = ConvexCloud[, 1],
          y = ConvexCloud[, 3],
          source = "A",
          type = "scatter",
          mode = "markers",
          # symbol="square",
          # color = I("black"),
          # color = I(colors),
          legendgroup = 'UMAP point cloud', showlegend = F,
          name = 'UMAP point cloud',
          text = c(rownames(ConvexCloud)),
          hoverinfo = 'text',
          hoveron = 'points+fills',
          # marker = list(
          #   size = 5,
          #   width = 2,
          #   line = list(
          #     # color = 'gray50',
          #     color = I(myvals$select_data_colors[,4]),
          #     width = 0.5)
          # )
          marker = list(color = I(myvals$select_data_colors[,4]), opacity = 1, size = 8)
          ,
          width = 735, height = 720
        )
        
        # fig2 <- fig2 %>% add_polygons(polygon, x = polygon[, 1], y = polygon[, 3], color=I("red"), opacity=0.2)
        fig2 <- fig2 %>% add_polygons(polygon, x = polygon[, 1], y = polygon[, 3],  
                                      fillcolor = 'rgba(206, 211, 214, 0.1)', hoveron = 'points+fills',
                                      legendgroup = 'CIE L* a* b*', showlegend = F ,
                                      name = 'CIE L* a* b*',
                                      marker = list(color = I(polygon_colors[,4]), opacity = 1, size = 10, symbol= 'cross', line = list(
                                        color = 'black', width=0.5)), 
                                      text = c(polygon_colors[,4]), hoverinfo = 'text', 
                                      line = list(color = '#a5c0cf'))
        
        fig2 <- fig2 %>% 
          layout(dragmode = "select")
        
        fig2 <- fig2 %>% layout(
          xaxis = list(title = "L"),
          yaxis = list(title = "b"))
        
        
        
          
        fig <- subplot(fig1, fig2, nrows = 2, margin = 0.07, titleY = TRUE, shareY = F, shareX = TRUE)
        
        fig
        
        
}
    })
    
    
    
    
    
    
    
    #--- Legend ---#
    output$legend <- renderDataTable({
      ConvexCloud <- myvals$NewUMAP
      ConvexCloud <- as.data.frame(ConvexCloud)
      
      Lab <- ConvexCloud
      Lab <- round(Lab, 2)
      rawdata = structure(
        list(
          Lstar = c(Lab[, 1]),
          Astar = c(Lab[, 2]),
          Bstar = c(Lab[, 3])
        ),
        .Names = c("Lstar", "Astar", "Bstar"),
        row.names = c(rownames(ConvexCloud)),
        class = "data.frame"
      )
      
      library(colorspace)
      LABdata <- with(rawdata, LAB(Lstar, Astar, Bstar))
      
      legend_colors <- as.data.frame(cbind(ConvexCloud,hex(LABdata, fix = TRUE)))
      # print(legend_colors)
      
      if (length(myvals$select_data)==0) {
        myvals$legend_genes <- myvals$NewUMAP
      } else {
        # req(input$remove_genes)
        if(myvals$select_data[,1]==2){
          comparison <- row.match(ConvexCloud[,c(1,3)],myvals$select_data[,3:4])
        } else comparison <- row.match(ConvexCloud[,c(1,2)],myvals$select_data[,3:4])
        
        myvals$legend_genes <- ConvexCloud[c(which(!is.na(comparison))), ]
        legend_colors <- as.data.frame(cbind(legend_colors, as.data.frame(comparison)))
        legend_colors <- subset(legend_colors, !is.na(legend_colors[,5]))
      }
      if(input$remove_genes){
        remove_genes()
        myvals$legend_genes <- myvals$RemoveGenesFromConvexCloud
        # legend_colors <- legend_colors[c(row.match(legend_colors[,c(1:3)], myvals$legend_genes)), ]
        
      }
      if(input$reset_genes){
        if(input$remove_genes[1]<=input$reset_genes[1]){
          reset()
          myvals$legend_genes <- myvals$NewUMAP
        }
        
      }
      
      
      convex_colors <- as.data.frame(cbind(rownames(legend_colors), legend_colors[,4]))
      
      # sweep(t(col2rgb(c('#fff000', '#000fff', '#45738a'))), MARGIN=2, c(0.2126, 0.7152, 0.0722), `*`)
      
      convex_colors <- cbind(convex_colors, brightness = rowSums(sweep(t(col2rgb(c(legend_colors[,4]))), MARGIN=2, c(0.2126, 0.7152, 0.0722), `*`)))
      # print(convex_colors)
      
      options(DT.options = list(pageLength = 25))
      df = as.data.frame(convex_colors)
      colnames(df) <- c("Names", "Colors", "brightness")
      
      datatable(df, rownames = FALSE, extensions = 'Responsive', selection = 'none', options = list(columnDefs = list(list(targets = 2, visible = FALSE)))) %>% formatStyle(colnames(df), 'Names', # target = 'row',
        backgroundColor = styleEqual(c(convex_colors[,1]), c(convex_colors[,2])), fontWeight = "bold"
        # , fontSize = '200%'
      ) %>%
        formatStyle(
          'brightness',
          target = 'row',
          color = styleInterval(40, c('gray', 'black'))
        )
      
      
    })
    
    
    #--- Downloads ---#
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('Hex_codes-', input$file1, "-" ,Sys.Date(), '.tsv', sep='')
      },
      content = function(con) {
        write.table(myvals$DownloadFile, con, quote = F, row.names = F, col.names = F, sep = "\t")
      }
    )
    
    output$download_table <- renderDataTable({
      datatable(myvals$DownloadFile, selection = c("none"), colnames = c("Genes", "Colors - Hex codes"))
    })
    
    observeEvent(input$openModal, {
      showModal(
        modalDialog(title = "Some title",
                    p("Some information"))
      )
    })
})
