#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Define UI for application that draws a histogram
source('~/Documents/Documents – SUN1012692/GitHub/CIELAB/global.R', local = TRUE)
options(shiny.maxRequestSize = 100*1024^2)
# shinyUI(

sidebar <- dashboardSidebar(
    # sidebarMenu(id="tabs", sidebarMenuOutput("menu"))
    sidebarMenu(menuItem("Upload Files", tabName = "upload", icon = icon("dashboard")),
                menuItem("UMAP", tabName = "umap", icon = icon("dashboard")),
                menuItem("Satellites", icon = icon("th"), tabName = "satellites",
                         badgeLabel = "new", badgeColor = "green"),
                menuItem("Download", icon = icon("th"), tabName = "Download")
    )
)

  
body <- 
    dashboardBody(
        fluidPage(
    
    # Application title
    # titlePanel("U-CIE"),
    titlePanel(div(HTML("<strong>U-CIE</strong> <em> [/juː 'siː/] </em> "))),
    
    
    
    tabItems(tabItem(tabName = "upload",
                     fileInput("file1", "Upload File",
                               multiple = FALSE,
                               accept = c("text/csv",
                                          "text/comma-separated-values,text/plain",
                                          ".csv")),
                     
                     selectInput(inputId = "dataset",
                                 label = "Example datasets:",
                                 choices = c("-" ,"GSE75748_time_course", "GSE75748_cell_type")),
                     
                     
                     
                     # Horizontal line ----
                     tags$hr(),
                     
                     # Input: Checkbox if file has header ----
                     selectInput("matrix", "Type of Matrix:",
                                  choices = c('-' = "-",
                                              'Single-cells' = 'Single-cells',
                                              'High Dimensional' = 'High Dimensional',
                                              'Distance matrix'=  'Distance matrix',
                                              '3D data' = '3D data'),
                                  selected = '-'),
                     
                     # Input: Select separator ----
                     checkboxInput("header", "Header", TRUE),
                     
                     radioButtons("sep", "Separator",
                                  choices = c(Comma = ",",
                                              Semicolon = ";",
                                              Tab = "\t"),
                                  selected = "\t"),
                     
                     # Input: Select quotes ----
                     radioButtons("quote", "Quote",
                                  choices = c(None = "",
                                              "Double Quote" = '"',
                                              "Single Quote" = "'"),
                                  selected = '"'),
                     
                     # Horizontal line ----
                     tags$hr(),
                     
                     # Input: Select number of rows to display ----
                     radioButtons("disp", "Display",
                                  choices = c(Head = "head",
                                              All = "all"),
                                  selected = "head"),
                     tags$hr(),
                     uiOutput("uploaded_dataset"),
                     tableOutput("contents")
                     
                     
    ),
    
    tabItem(tabName = "umap",
                     # h2("Parameters and UMAP"),
            column(2, sliderInput("weightL",
                                 "L* weight (Brightness):",
                                 min = 1,
                                 max = 3,
                                 value = 1,
                                 step = 0.1, 
                                 width = '1000px'),
                     sliderInput("weightA",
                                 "a* weight (Red - Green):",
                                 min = 1,
                                 max = 3,
                                 value = 1,
                                 step = 0.1, 
                                 width = '1000px'),
                     sliderInput("weightB",
                                 "b* weight (Yellow - Blue):",
                                 min = 1,
                                 max = 3,
                                 value = 1,
                                 step = 0.1, 
                                 width = '1000px'),
                     sliderInput("scaling",
                                 "Scaling factor multiplier:",
                                 min = 1,
                                 max = 2,
                                 value = 1,
                                 step = 0.1, 
                                 width = '1000px'),
                     actionButton("weightButton", "Weight", width = '350px'),
                     actionButton("scalingButton", "Scale", width = '350px'),
                   br(),
                   br(),
                   br(),
                   br(),
                   dataTableOutput("table")),
                    
            column(1, offset = 2, plotlyOutput("plotly_plot")),
                     uiOutput("list_of_parameters"),
                     
                     
                     
                     
    ),
    tabItem(tabName = "satellites",
            # h2("Satellites"),
            column(2, 
                   actionButton("remove_genes", "Remove Genes and Re-color!"),
                   actionButton("reset_genes", "Reset"),
                   uiOutput("reset")
                   ),
            column(1, plotlyOutput("satellite1"))
            
    ),
    tabItem(tabName = "Download",
            column(1,downloadButton('downloadData', 'Download')),
            hr(),
            hr(),
            hr(),
            dataTableOutput("download_table")

    )
    
    )

        
        # Show a plot of the generated distribution
        # mainPanel(
            
            # plotlyOutput("plotly_plot"),
            # uiOutput("list_of_parameters"),
            # dataTableOutput("table"),
            # 
            # plotOutput("satellite1"),
            # plotOutput("satellite2")
            # verbatimTextOutput("summary")
            
        
    # )
))

dashboardPage(skin = "purple",
    dashboardHeader(title = "U-CIE",
                    tags$li(actionLink("openModal", label = "", icon = icon("info")),
                            class = "dropdown")
                    
    ),
    sidebar,
    body
)

