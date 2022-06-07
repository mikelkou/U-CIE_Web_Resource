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

# future_promise({
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
      if(input$matrix == 'Single-cells'){
        if (is.null(input$file1)) {
          return()
        } else {
          
        if(grepl('xls', input$file1$name, fixed = TRUE)){
          dataset1 <- as.data.frame(read_excel(input$file1$datapath, 1, col_names = ifelse(input$header1==T, T, F)))
        } 
        if(grepl('csv', input$file1$name, fixed = TRUE)){
          dataset1 <- read.table(input$file1$datapath, header = input$header1, sep = ",")
        } 
        if(grepl('txt', input$file1$name, fixed = TRUE)){
          dataset1 <- read.table(input$file1$datapath, header = input$header1, sep = "\t")
        }
        if(grepl('tsv', input$file1$name, fixed = TRUE)){
          dataset1 <- read.table(input$file1$datapath, header = input$header1, sep = "\t")
        }
          #Missing values
          if(any(is.na(dataset1)) || any(dataset1 =="")){
            showModal(modalDialog(title = "There are missing values in your dataset. Please check again!", easyClose = T, fade = T))
            dataset1 <- NULL
          }
        return(dataset1)
        }
        
      }
      
      if(input$matrix == 'High Dimensional'){
        if (is.null(input$file2)) {
          return()
        } else {
        if(grepl('xls', input$file2$name, fixed = TRUE)){
          dataset2 <- as.data.frame(read_excel(input$file2$datapath, 1, col_names = ifelse(input$header2==T, T, F)))
        } 
        if(grepl('csv', input$file2$name, fixed = TRUE)){
          dataset2 <- read.table(input$file2$datapath, header = input$header2, sep = ",")
        } 
        if(grepl('txt', input$file2$name, fixed = TRUE)){
          dataset2 <- read.table(input$file2$datapath, header = input$header2, sep = "\t")
        }
        if(grepl('tsv', input$file2$name, fixed = TRUE)){
          dataset2 <- read.table(input$file2$datapath, header = input$header2, sep = "\t")
        }
          #Missing values
          if(any(is.na(dataset2)) || any(dataset2 =="")){
            showModal(modalDialog(title = "There are missing values in your dataset. Please check again!", easyClose = T, fade = T))
            dataset2 <- NULL
          }
          return(dataset2)
        } 
      }
      
      if(input$matrix == 'Distance matrix'){
        if (is.null(input$file3)) {
          return()
        } else {
        if(grepl('xls', input$file3$name, fixed = TRUE)){
          dataset3 <- as.data.frame(read_excel(input$file3$datapath, 1, col_names = ifelse(input$header3==T, T, F)))
        } 
        if(grepl('csv', input$file3$name, fixed = TRUE)){
          dataset3 <- read.table(input$file3$datapath, header = input$header3, sep = ",")
        } 
        if(grepl('txt', input$file3$name, fixed = TRUE)){
          dataset3 <- read.table(input$file3$datapath, header = input$header3, sep = "\t")
        }
        if(grepl('tsv', input$file3$name, fixed = TRUE)){
          dataset3 <- read.table(input$file3$datapath, header = input$header3, sep = "\t")
        }
          #Missing values
          if(any(is.na(dataset3)) || any(dataset3 =="")){
            showModal(modalDialog(title = "There are missing values in your dataset. Please check again!", easyClose = T, fade = T))
            dataset3 <- NULL
          }
        return(dataset3)
        }
      }
      
      if(input$matrix == '3D data'){
        if (is.null(input$file4)) {
          return()
        } else {
          
        if(grepl('xls', input$file4$name, fixed = TRUE)){
          dataset4 <- as.data.frame(read_excel(input$file4$datapath, 1, col_names = ifelse(input$header4==T, T, F)))
        } 
        if(grepl('csv', input$file4$name, fixed = TRUE)){
          dataset4 <- read.table(input$file4$datapath, header = input$header4, sep = ",")
        } 
        if(grepl('txt', input$file4$name, fixed = TRUE)){
          dataset4 <- read.table(input$file4$datapath, header = input$header4, sep = "\t")
        }
        if(grepl('tsv', input$file4$name, fixed = TRUE)){
          dataset4 <- read.table(input$file4$datapath, header = input$header4, sep = "\t")
        }
        
      if(ncol(dataset4)<3){
        showModal(modalDialog(title = "Not even 3D? Please check again the dataset.", easyClose = T, fade = T))
        dataset4 <- NULL
      }
          #Missing values
          if(any(is.na(dataset4)) || any(dataset4 =="")){
            showModal(modalDialog(title = "There are missing values in your dataset. Please check again!", easyClose = T, fade = T))
            dataset4 <- NULL
          }
        return(dataset4)
        }
      }
      
    })
    
    loadNetworkFromFile <- function() {
      dataset1 <- NULL
      dataset2 <- NULL
      dataset3 <- NULL
      dataset4 <- NULL
      set.seed(123)
      
      switch(
        input$LoadFileSingleCellsInput,
        oF = {
          if (!is.null(input$file1)) {
            if(input$file1$type == "application/pdf"){
              # "This is PDF! Please upload tsv/csv files."
              showModal(modalDialog(title = "This is PDF! Please upload tsv/txt/csv/xls files.", easyClose = T, fade = T))
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
          # dataset1 <- data.frame(log2FC_MatrixCounts_ct)
          # dataset1 <- log2FC_MatrixCounts_ct[1:500,1:500]
          # log2FC_MatrixCounts <- log2FC_MatrixCounts_ct
          # log2FC_MatrixCounts <- log2FC_MatrixCounts_ct[1:1000,1:1000]
          # log2FC_MatrixCounts_ct <- data.frame(read.delim("~/Documents/Documents – SUN1012692/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
          # log2FC_MatrixCounts <- data.frame(read.delim("~/Documents/GitHub/CIELAB/log2FC_MatrixCounts_ct.tsv"))
        }
      )
      
      switch(
        input$LoadFileHighDInput,
        oF = {
          if (!is.null(input$file2)) {
            if(input$file2$type == "application/pdf"){
              # "This is PDF! Please upload tsv/csv files."
              showModal(modalDialog(title = "This is PDF! Please upload tsv/txt/csv/xls files.", easyClose = T, fade = T))
            } else {
              dataset2 <- read_data(input$file2$datapath)
            }
          }
        }
        ,
        oR_Example_highd = {
          # dataset2 <- data.frame(read.delim("GSE75748_cropped_sc_cell_type_ec.tsv"))
          data(fruitfly)
          dataset2 <- as.data.frame(fruitfly)
        }
      )
      
      switch(
        input$LoadFileDistInput,
        oF = {
          if (!is.null(input$file3)) {
            if(input$file3$type == "application/pdf"){
              # "This is PDF! Please upload tsv/csv files."
              showModal(modalDialog(title = "This is PDF! Please upload tsv/txt/csv/xls files.", easyClose = T, fade = T))
            } else {
            dataset3 <- read_data(input$file3$datapath)
            }
          }
        },
        oR_Example_distance = {
          dataset3 <- data.frame(read.delim("TreeOfLife.tsv"))
        }
      )
      
      switch(
        input$LoadFile3DInput,
        oF = {
          if (!is.null(input$file4)) {
            if(input$file4$type == "application/pdf"){
              # "This is PDF! Please upload tsv/csv files."
              showModal(modalDialog(title = "This is PDF! Please upload tsv/txt/csv/xls files.", easyClose = T, fade = T))
            } else {
            dataset4 <- read_data(input$file4$datapath)
            }
          }
        },
        oR_Example3 = {
          dataset4 <- data.frame(read.delim("fruitfly_lle_dim_red.tsv", header = T))
        }
      )
      
      enable("btnanalysis")
      
      if(input$matrix == 'Single-cells'){
        if (input$LoadFileSingleCellsInput == "oF" && is.null(input$file1)) {
          return()
        } else return(dataset1)
      }
      
      if(input$matrix == 'High Dimensional'){
        if (input$LoadFileHighDInput == "oF" && is.null(input$file2)) {
          return()
        } else return(dataset2)
      }
        
      if(input$matrix == 'Distance matrix'){
        if (input$LoadFileDistInput == "oF" && is.null(input$file3)) {
          return()
        } else return(dataset3)
      }
        
      if(input$matrix == '3D data'){
        if (input$LoadFile3DInput == "oF" && is.null(input$file4)) {
          return()
        } else return(dataset4)
      }

      
      # return(dataset1)
    }
    
    
    output$LoadFileSingleCellsOutput <- renderUI({
      if (input$LoadFileSingleCellsInput == "oF" ){
        wellPanel(fileInput(
          "file1",
          "Choose file to upload",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv",
            ".xls"
          ))
        )
      } else{}
    })
    
    output$LoadFileHighDOutput <- renderUI({  
      if (input$LoadFileHighDInput == "oF" ){
        wellPanel(fileInput(
          "file2",
          "Choose file to upload",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv",
            ".xls"
          ))
        )
      } else{}
    })
    
    output$LoadFileDistOutput <- renderUI({
      if (input$LoadFileDistInput == "oF" ){
        wellPanel(fileInput(
          "file3",
          "Choose file to upload",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv",
            ".xls"
          ))
        )
      } else{}
    })
    
    output$LoadFile3DOutput <- renderUI({
      if (input$LoadFile3DInput == "oF" ){
        wellPanel(fileInput(
          "file4",
          "Choose file to upload",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv",
            ".xls"
          ))
        )
      } else{}
    })
    
    ###### Upload #######
    doAddNetwork <- observeEvent(input$btnanalysis, {
      req(input$matrix)
      myvals$uploaded_df <- loadNetworkFromFile()
      
    })
    
    
    output$contents <- renderDataTable({
      req(input$matrix)
      
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
        datatable(df[1:10,], caption = paste("Total number of rows: ", nrow(df), "  |  ", "Total number of columns: ", ncol(df)), 
                  options = list(pageLength = 10, scrollX = TRUE))
      
    })
    
    
    observe({
      disable("btnanalysis")
      req(input$btnanalysis)
      enable("btnanalysis")
      
      isolate({
      if(is.null(myvals$uploaded_df)) {
        return()
      }
      
        # if(input$uiLoadGraphOptionsInput=="oF") {
          req(input$matrix)
          isolate({
            show_modal_spinner(spin = "circle", text = "Please wait..." )
        if(is.null(input$matrix)){
          return()
        }
          
        if(input$matrix == 'Single-cells'){
          
          if(input$LoadFileSingleCellsInput != "oR_Example1"){
          matrix_counts <- loadNetworkFromFile()
          
          if(nrow(matrix_counts) < 50){
            print(nrow(matrix_counts))
            showModal(modalDialog(title = "The limitation of 50 cells is coming from the downstream of the package.", easyClose = T, fade = T))
            return() 
          }
          
          if(!is.character(matrix_counts[,1]) && !is.na(as.integer(rownames(matrix_counts))) ){
            showModal(modalDialog(title = "No cell names (rownames) names present in the input matrix.", easyClose = T, fade = T))
            return()
          }

          if(is.character(matrix_counts[,1])){
            if(length(unique(matrix_counts[,1])) != length(matrix_counts[,1])){
              matrix_counts[,1] = paste(matrix_counts[,1], 1:nrow(matrix_counts), sep = "_")
              rownames(matrix_counts) <- matrix_counts[,1]
              matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
            } else{
            rownames(matrix_counts) <- matrix_counts[,1]
            matrix_counts <- matrix_counts[,2:ncol(matrix_counts)]
            }
          }
          
          
          if(length(matrix_counts[matrix_counts < 0])){
            print(matrix_counts[matrix_counts < 0])
            showModal(modalDialog(title = "The 'Single cells' track expects expression matrix with counts and therefore, 0 or positive values. If your matrix contains negative values (e.g. after log normalization), 
            please consider the 'High-dimensional' track and transpose your matrix based on your interest.
                                  )", easyClose = T, fade = T))
            return() 
          }
          
          withConsoleRedirect("console", {
          matrix_counts <- data.matrix(matrix_counts) # must be a matrix object!
          TransposedMatrixCounts <- t(matrix_counts) # Correct!
              
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
          
          if(input$LoadFileSingleCellsInput == "oR_Example1"){
            input$btnanalysis
            
            isolate({
              show_modal_spinner(spin = "circle", text = "Please wait..." )
              
              withConsoleRedirect("console", {
                log2FC_MatrixCounts <- loadNetworkFromFile()
                
                # print("1")
                SDMoreThanOnelog2FC_MatrixCounts <- c()
                for(i in 1:ncol(log2FC_MatrixCounts)){
                  std <- sd(log2FC_MatrixCounts[,i])
                  SDMoreThanOnelog2FC_MatrixCounts <- cbind(SDMoreThanOnelog2FC_MatrixCounts, std)
                }
                colnames(SDMoreThanOnelog2FC_MatrixCounts) <- colnames(log2FC_MatrixCounts)
                
                SDMoreThanOnelog2FC_MatrixCounts <- dplyr::select(as.data.frame(log2FC_MatrixCounts),
                                                           -c(names(which(SDMoreThanOnelog2FC_MatrixCounts[1,]<1))))
                
                # print(SDMoreThanOnelog2FC_MatrixCounts)
                
                data <- CreateSeuratObject(counts = SDMoreThanOnelog2FC_MatrixCounts)
                all.genes <- rownames(data)
                data <- ScaleData(data, do.scale =  F, do.center = F, features = all.genes)
                data <- FindVariableFeatures(object = data, selection.method = 'mvp') #mvp because of error in log
                data <- RunPCA(data, npcs = 50, features = VariableFeatures(object = data))
                #print(data)
                data <- FindNeighbors(data, dims = 1:15)
                data <- FindClusters(data, resolution = 0.5, algorithm= 1) # color in Seurat umap output
                data <- RunUMAP(data, dims = 1:50, n.components = 3L)
                
                data_umap_coord <- as.data.frame(data[["umap"]]@cell.embeddings)
                # data_umap_coord <- data[["umap"]]@cell.embeddings
              }) # console
              umap_dist <- data_umap_coord
              myvals$umap_dist <- umap_dist
            }) # isolate
          }
        } # if input$matrix == 'Single-cells'
      
      if(input$matrix == 'High Dimensional'){
        if(input$LoadFileHighDInput == "oR_Example_highd"){
          # 2
          input$btnanalysis
          
          isolate({
            show_modal_spinner(spin = "circle", text = "Please wait..." )
            withConsoleRedirect("console", {
              data <- uwot::umap(loadNetworkFromFile(), ret_nn = TRUE, n_neighbors = 5, n_components = 3) # library(uwot)
              data_umap_coord <- as.data.frame(data$embedding)
              rownames(data_umap_coord) <- rownames(loadNetworkFromFile())
            })
            print("oR_Example_highd")
          }) # isolate
        } #if if example
        
        else{
          isolate({
            df = loadNetworkFromFile()
            
            if(is.character(df[,1])){
              if(length(unique(df[,1])) != length(df[,1])){
                df[,1] = paste(df[,1], 1:nrow(df), sep = "_")
                rownames(df) <- df[,1]
                df <- df[,2:ncol(df)]
              } else{
                rownames(df) <- df[,1]
                df <- df[,2:ncol(df)]
              }
            }
            
            # if(is.character(df[,1])){
            #   rownames(df) <- df[,1]
            #   df <- df[,2:ncol(df)]
            # }
            
            show_modal_spinner(spin = "circle", text = "Please wait..." )
        withConsoleRedirect("console", {
        data <- uwot::umap(loadNetworkFromFile(), ret_nn = TRUE, n_neighbors = 5, n_components = 3) # library(uwot)
        data_umap_coord <- as.data.frame(data$embedding)
        rownames(data_umap_coord) <- rownames(df)
        })
          }) # isolate
        } #else
        umap_dist <- data_umap_coord
        myvals$umap_dist <- umap_dist
        # }
      }
            
      if(input$matrix == 'Distance matrix'){
        if(ncol(loadNetworkFromFile())!=nrow(loadNetworkFromFile())){
          withConsoleRedirect("console", {
          print(paste("Number of rows: ", nrow(loadNetworkFromFile())))
          print(paste("Number of columns: ", ncol(loadNetworkFromFile())))
          })
          
          showModal(modalDialog(title = paste("Distance matrix should have the same number of rows and columns (NxN matrix). ", 
                                "Number of rows: ", nrow(loadNetworkFromFile()), 
                                "Number of columns: ", ncol(loadNetworkFromFile()))
                                , 
                                easyClose = T, fade = T))
          
          
          
          return()
         
        } else {
          withConsoleRedirect("console", {
            data <- as.matrix(dist(loadNetworkFromFile()))
            data <- umap(data, input="dist", n_components = 3)
            data_umap_coord <- as.data.frame(data$layout)
            umap_dist <- data_umap_coord
            myvals$umap_dist <- umap_dist
          })
          print("Distance")
        }
        # }
      }
            
            if(input$matrix == "3D data"){
              if(input$LoadFile3DInput != "oR_Example3"){
                
                df <- c()
                for(i in 1:ncol(loadNetworkFromFile())){

                  if(is.character(loadNetworkFromFile()[,i])){
                    df <- loadNetworkFromFile()[,-i]
                    
                    if(length(unique(loadNetworkFromFile()[,i])) != length(loadNetworkFromFile()[,i])){
                      new_names = paste(loadNetworkFromFile()[,i], 1:nrow(loadNetworkFromFile()), sep = "_")
                      rownames(df) <- new_names
                    } else{
                      rownames(df) <- loadNetworkFromFile()[,1]
                    }
                  }
                }
                
                if(is.null(df)){
                  df <- loadNetworkFromFile()
                }
                
                if(ncol(df)!=3){
                  showModal(modalDialog(title = "If more than 3 numeric columns (+1 optional column with names), please check again the dataset.", easyClose = T, fade = T))
                  return()
                }
                
                myvals$umap_dist <- df
                myvals$NewUMAP <- df
                print("3D")
                
              }
              
              if(input$LoadFile3DInput == "oR_Example3"){
                # 2
                input$btnanalysis
                
                isolate({
                  show_modal_spinner(spin = "circle", text = "Please wait..." )
                  withConsoleRedirect("console", {
                    df <- loadNetworkFromFile()
                    # rownames(df) <- loadNetworkFromFile()[,1]
                    # df <- df[,2:ncol(df)]
                    myvals$umap_dist <- df
                    myvals$NewUMAP <- df
                  }) 
                }) # isolate
                
              } # "oR_example3
              
            } # 3D data
          }) # isolate
          # } # of
        
        
      # myvals$umap_dist_stored <- myvals$umap_dist
        enable("btnanalysis")
        remove_modal_spinner()
    
      
    }) # isolate analysis
    })
    
    
    output$css_upload <- renderUI({
      if(input$tabs=='upload'){
      tags$style(".fa-file-upload {color:#fcd049}")
      }
    })
    output$css_umap <- renderUI({
      if(input$tabs=='umap'){
      tags$style(".fa-cube {color:#fcd049}")
      }
    })
    output$css_satellites <- renderUI({
      if(input$tabs=='satellites'){
      tags$style(".fa-bar-chart-o {color:#fcd049}")
      }
    })
    output$css_console <- renderUI({
      if(input$tabs=='console'){
      tags$style(".fa-terminal {color:#fcd049}")
      }
    })
    output$css_download <- renderUI({
      if(input$tabs=='Download'){
      tags$style(".fa-file-export {color:#fcd049}")
      }
    })
    
    
    
    
    # Change tab when UMAP is ready
    switch_tabs <- observeEvent(input$btnanalysis, {
      if(is.null(loadNetworkFromFile())){
        return()
      }
      newtab <- switch(input$tabs,
                       "upload" = "umap",
                       "umap" = "upload"
      )
      updateTabItems(session, "tabs", newtab)
      
    })
    
    
    # Disable weight button
    ReactiveValuesOfAxesWeightL <- reactiveValues(prev_bins = NULL)
    ReactiveValuesOfAxesWeightA <- reactiveValues(prev_bins = NULL)
    ReactiveValuesOfAxesWeightB <- reactiveValues(prev_bins = NULL)
    # Append new value to previous values when input$bins changes 
    observeEvent(input$weightButton, {
      ReactiveValuesOfAxesWeightL$prev_bins <- c(1, ReactiveValuesOfAxesWeightL$prev_bins, input$weightL)
      ReactiveValuesOfAxesWeightA$prev_bins <- c(1, ReactiveValuesOfAxesWeightA$prev_bins, input$weightA)
      ReactiveValuesOfAxesWeightB$prev_bins <- c(1, ReactiveValuesOfAxesWeightB$prev_bins, input$weightB)
    })
    
    observe({
      if(length(ReactiveValuesOfAxesWeightL$prev_bins)==0){
        return()
      }
      if(tail(ReactiveValuesOfAxesWeightL$prev_bins,1) == input$weightL
         && tail(ReactiveValuesOfAxesWeightA$prev_bins,1) == input$weightA
         && tail(ReactiveValuesOfAxesWeightB$prev_bins,1) == input$weightB
      ) {
        disable("weightButton")
      } 
      if(tail(ReactiveValuesOfAxesWeightL$prev_bins,1) != input$weightL
         || tail(ReactiveValuesOfAxesWeightA$prev_bins,1) != input$weightA
         || tail(ReactiveValuesOfAxesWeightB$prev_bins,1) != input$weightB
         ) {
        # print(input$weightL)
        enable("weightButton")
      }
    })
     
    # Disable scale button
    ReactiveValuesOfAxesScale <- reactiveValues(prev_bins = NULL)
    observeEvent(input$scalingButton, {
      ReactiveValuesOfAxesScale$prev_bins <- c(ReactiveValuesOfAxesScale$prev_bins, input$scaling)
    })

    observe({
      if(length(ReactiveValuesOfAxesScale$prev_bins)==0){
        return()
      }
      if(tail(ReactiveValuesOfAxesScale$prev_bins,1) == input$scaling) {
        disable("scalingButton")
      }
      if(tail(ReactiveValuesOfAxesScale$prev_bins,1) != input$scaling) {
        enable("scalingButton")
      }
    })
    
    
    observe({
      umap_dist <- myvals$umap_dist 
      if(is.null(umap_dist)){}
      else{
        show_modal_spinner(spin = "circle", text = "Please wait..." )
        
        if(!is.null(myvals$RemoveGenesFromConvexCloud)){
          req(input$remove_genes)
          umap_dist <- myvals$RemoveGenesFromConvexCloud
        }
        
        input$weightButton
        WL <- isolate(input$weightL)
        Wa <- isolate(input$weightA)
        Wb <- isolate(input$weightB)
        
        # simplex_vectors <- FitColorsFunction(umap_dist, polygon, WL, Wa, Wb)
        # plan(multisession, workers = 16L)
        promise <- future_promise(FitColorsFunction(umap_dist, polygon, WL, Wa, Wb))
        then(promise, function(simplex_vectors) {
          myvals$start.values <- simplex_vectors[which(simplex_vectors[,1] == max(simplex_vectors[,1]))[1], ]
          myvals$simplex_vectors <- simplex_vectors
          remove_modal_spinner()
        }) # then
      } # else
    })
    
    
    
    
    output$table <- renderDataTable({
      # req(input$weightButton)
      start.values <- myvals$start.values
      if(is.null(start.values)){
        return()
      }
      
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
      input$table_rows_selected
      
      if(is.null(input$table_rows_selected)){
        # stop("")
        tryCatch({
          stop("")
        }, error = function(e) {
          shiny:::reactiveStop(conditionMessage(e))
        })
      }
      
      # input$table2_rows_selected
      start.values <- myvals$start.values
      simplex_vectors <- myvals$simplex_vectors
      
      s <- input$table_rows_selected
      # s2 <- input$table2_rows_selected
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
        # Table on satellites as well
        # if(!is.null(s2)){
        #   s <- s2
        # }
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
    DownloadFile <- as.data.frame(cbind(DownloadFile, LABdata@coords))
    rownames(DownloadFile) <- NULL
    colnames(DownloadFile) <- c("Genes", "Colors - Hex codes", "Lstar", "Astar", "Bstar")
    
    if(!is.null(myvals$RemovedGenes)){
      DownloadFile <- rbind(DownloadFile, myvals$RemovedGenes)
    }
    myvals$DownloadFile <- DownloadFile
    
    axx <- list(nticks = 4,
                title = "L*")
    
    axy <- list(nticks = 4,
                title = "a*")
    
    axz <- list(nticks = 4,
                title = "b*")
    
    
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
        width = 1000, height = 720
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
      select_data <- event_data("plotly_selected", source = "A")
      
      myvals$select_data <- select_data
      
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
        
        if(length(myvals$select_data)==0 && length(myvals$select_data2)==0) {
          myvals$select_data_colors <- as.data.frame(cbind(myvals$colors, rep("gray50", nrow(myvals$colors))))
          }
        if(!(length(myvals$select_data)==0) && length(myvals$select_data2)==0){
          # if(myvals$select_data[,1]==2){
            comparison <- row.match(myvals$colors[,c(1,2)],myvals$select_data[,3:4])
          # } else comparison <- row.match(myvals$colors[,c(1,2)],myvals$select_data[,3:4])
          
          a <- as.data.frame(cbind(myvals$colors, as.data.frame(comparison)))
          a[,4] <- ifelse(is.na(a[,5]),  "gray50",  a[,4])
          myvals$select_data_colors <- a
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
          # ,
          # width = 735, height = 720
        )
        
        fig1 <- fig1 %>% add_polygons(polygon, x = polygon[, 1], y = polygon[, 2], 
                                      fillcolor = 'rgba(206, 211, 214, 0.1)', 
                                      name = 'CIE L* a* b*',
                                      legendgroup = 'CIE L* a* b*', showlegend = T,
                                      marker = list(color = I(polygon_colors[,4]), opacity = 1, size = 10, symbol= 'cross', line = list(
                                        color = 'black', width=0.5)),
                                      text = c(polygon_colors[,4]), hoverinfo = 'text', 
                                      line = list(color = '#a5c0cf'))
        
        fig1 <- fig1 %>% add_trace(ConvexCloud, x = ConvexCloud[, 1], y = ConvexCloud[, 2], showlegend = F)
        
        fig1 <- fig1 %>% 
          layout(dragmode = "select")
        
        fig1 <- fig1 %>% layout(
          xaxis = list(title = "L"),
          yaxis = list(title = "a"))
        
        # fig2 <- plot_ly(
        #   data = as.data.frame(ConvexCloud),
        #   x = ConvexCloud[, 1],
        #   y = ConvexCloud[, 3],
        #   source = "A",
        #   type = "scatter",
        #   mode = "markers",
        #   legendgroup = 'UMAP point cloud', showlegend = F,
        #   name = 'UMAP point cloud',
        #   text = c(rownames(ConvexCloud)),
        #   hoverinfo = 'text',
        #   hoveron = 'points+fills',
        #   marker = list(color = I(myvals$select_data_colors[,4]), opacity = 1, size = 8)
        #   # ,
        #   # width = 735, height = 720
        # )
        # 
        # fig2 <- fig2 %>% add_polygons(polygon, x = polygon[, 1], y = polygon[, 3],  
        #                               fillcolor = 'rgba(206, 211, 214, 0.1)', hoveron = 'points+fills',
        #                               legendgroup = 'CIE L* a* b*', showlegend = F ,
        #                               name = 'CIE L* a* b*',
        #                               marker = list(color = I(polygon_colors[,4]), opacity = 1, size = 10, symbol= 'cross', line = list(
        #                                 color = 'black', width=0.5)), 
        #                               text = c(polygon_colors[,4]), hoverinfo = 'text', 
        #                               line = list(color = '#a5c0cf'))
        # 
        # fig2 <- fig2 %>% add_trace(ConvexCloud, x = ConvexCloud[, 1], y = ConvexCloud[, 3], showlegend = F)
        # 
        # fig2 <- fig2 %>%
        #   layout(dragmode = "select")
        # 
        # fig2 <- fig2 %>% layout(
        #   xaxis = list(title = "L"),
        #   yaxis = list(title = "b"))
        #   
        # fig <- subplot(fig1, fig2, nrows = 2, margin = 0.07, titleY = TRUE, shareY = F, shareX = TRUE)
        fig <- fig1
        fig
}
    })
    
    
    output$satellite2 <- renderPlotly({
      select_data <- event_data("plotly_selected", source = "B")
      
      myvals$select_data2 <- select_data

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
        
        # if(length(myvals$select_data2)==0 && length(myvals$select_data)==0) {
        #   myvals$select_data_colors <- as.data.frame(cbind(myvals$colors, rep("gray50", nrow(myvals$colors))))
        # }
        
        if(!(length(myvals$select_data2)==0) && length(myvals$select_data)==0){
          # if(myvals$select_data[,1]==2){
            comparison <- row.match(myvals$colors[,c(1,3)],myvals$select_data2[,3:4])
          # } else comparison <- row.match(myvals$colors[,c(1,2)],myvals$select_data[,3:4])
          
          a <- as.data.frame(cbind(myvals$colors, as.data.frame(comparison)))
          a[,4] <- ifelse(is.na(a[,5]),  "gray50",  a[,4])
          myvals$select_data_colors <- a
        } 
        
        
        fig2 <- plot_ly(
          data = as.data.frame(ConvexCloud),
          x = ConvexCloud[, 1],
          y = ConvexCloud[, 3],
          source = "B",
          type = "scatter",
          mode = "markers",
          legendgroup = 'UMAP point cloud', showlegend = T,
          name = 'UMAP point cloud',
          text = c(rownames(ConvexCloud)),
          hoverinfo = 'text',
          hoveron = 'points+fills',
          marker = list(color = I(myvals$select_data_colors[,4]), opacity = 1, size = 8)
          # ,
          # width = 735, height = 720
        )
        
        fig2 <- fig2 %>% add_polygons(polygon, x = polygon[, 1], y = polygon[, 3],  
                                      fillcolor = 'rgba(206, 211, 214, 0.1)', hoveron = 'points+fills',
                                      legendgroup = 'CIE L* a* b*', showlegend = T ,
                                      name = 'CIE L* a* b*',
                                      marker = list(color = I(polygon_colors[,4]), opacity = 1, size = 10, symbol= 'cross', line = list(
                                        color = 'black', width=0.5)), 
                                      text = c(polygon_colors[,4]), hoverinfo = 'text', 
                                      line = list(color = '#a5c0cf'))
        
        fig2 <- fig2 %>% add_trace(ConvexCloud, x = ConvexCloud[, 1], y = ConvexCloud[, 3], showlegend = F)
        
        fig2 <- fig2 %>%
          layout(dragmode = "select")
        
        fig2 <- fig2 %>% layout(
          xaxis = list(title = "L"),
          yaxis = list(title = "b"))
        
        # fig <- subplot(fig1, fig2, nrows = 2, margin = 0.07, titleY = TRUE, shareY = F, shareX = TRUE)
        fig <- fig2
        fig
      }
    })
    
    
    
    #--- Legend ---#
    output$legend <- renderDataTable({
      if(length(myvals$NewUMAP)==0){
        return()
      }
      else{
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
      
      LABdata <- with(rawdata, LAB(Lstar, Astar, Bstar))
      
      legend_colors <- as.data.frame(cbind(ConvexCloud,hex(LABdata, fix = TRUE)))
      # print(legend_colors)
      
      
      if (length(myvals$select_data)==0 && length(myvals$select_data2)==0) {
        myvals$legend_genes <- myvals$NewUMAP
      } else {
        if(!is.null(myvals$select_data) && is.null(myvals$select_data2)){
          comparison <- row.match(ConvexCloud[,c(1,2)],myvals$select_data[,3:4])
        } 
        if(is.null(myvals$select_data) && !is.null(myvals$select_data2)){
         comparison <- row.match(ConvexCloud[,c(1,3)],myvals$select_data2[,3:4])
        }
        
        if(!is.null(myvals$select_data) && !is.null(myvals$select_data2)){
          comparison <- ConvexCloud
        }
        
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
      convex_colors <- cbind(convex_colors, brightness = rowSums(sweep(t(col2rgb(c(legend_colors[,4]))), MARGIN=2, c(0.2126, 0.7152, 0.0722), `*`))) # font color based on brightness
      
      options(DT.options = list(pageLength = 15))
      df = as.data.frame(convex_colors)
      colnames(df) <- c("Names", "Colors", "brightness")
      df <- na.omit(df)
      
      if(nrow(df)>=2000){
        df <- df[1:2000,]
      }
      datatable(df, rownames = FALSE, extensions = 'Responsive', selection = 'none', 
                options = list(columnDefs = list(list(targets = c(1,2), visible = FALSE)))) %>% formatStyle(colnames(df), 'Names', # target = 'row',
        backgroundColor = styleEqual(c(convex_colors[1:nrow(df),1]), c(convex_colors[1:nrow(df),2])), fontWeight = "bold"
        # , fontSize = '200%'
      ) %>%
        formatStyle(
          'brightness',
          target = 'row',
          color = styleInterval(50, c('gray', 'black'))
        )
      }
    })
    
    
    #--- Downloads ---#
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('Hex_codes-Lab_coords-', input$file1, "-" , Sys.Date(), '.tsv', sep='')
      },
      content = function(con) {
        write.table(myvals$DownloadFile, con, quote = F, row.names = F, col.names = T, sep = "\t")
      }
    )
    
    output$download_table <- renderDataTable({
      datatable(myvals$DownloadFile, selection = c("none"), colnames = c("Names", "Colors - Hex codes", "Lstar", "Astar", "Bstar"))
    })
    
    observeEvent(input$openModal, {
      showModal(
        modalDialog(title = "Contact info:",
                    p("Mikaela Koutrouli: mikaela.koutrouli@cpr.ku.dk"),
                      p("Lars Juhl Jensen: lars.juhl.jensen@cpr.ku.dk"), 
                      p("For better explanation of the idea, here you can watch a 7min talk about it:"), 
                    tags$a(href = 'https://youtu.be/V6KgC5KJ-3g', "U-CIE: Color encoding of high-dimensional data using the CIELAB color space and UMAP")
                    )
      )
    })
}) # THE END

# }) #future promises
