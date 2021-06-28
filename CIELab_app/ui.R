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

# shinyUI(

sidebar <- dashboardSidebar(
    # sidebarMenu(id="tabs", sidebarMenuOutput("menu"))
    sidebarMenu(menuItem("Upload Files", tabName = "upload", icon = icon("dashboard")),
                menuItem("UMAP", tabName = "umap", icon = icon("dashboard")),
                menuItem("Satellites", icon = icon("th"), tabName = "satellites",
                         badgeLabel = "new", badgeColor = "green")
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
                     sliderInput("weightL",
                                 "L* weight (Brightness):",
                                 min = 1,
                                 max = 3,
                                 value = 1,
                                 step = 0.1),
                     sliderInput("weightA",
                                 "a* weight (Red - Green):",
                                 min = 1,
                                 max = 3,
                                 value = 1,
                                 step = 0.1),
                     sliderInput("weightB",
                                 "b* weight (Yellow - Blue):",
                                 min = 1,
                                 max = 3,
                                 value = 1,
                                 step = 0.1),
                     sliderInput("scaling",
                                 "Scaling factor multiplier:",
                                 min = 1,
                                 max = 2,
                                 value = 1,
                                 step = 0.1),
                     actionButton("weightButton", "Weight"),
                     actionButton("scalingButton", "Scale"),
                     br(),
                     br(),
                     br(),
                     br(),
                     plotlyOutput("plotly_plot"),
                     uiOutput("list_of_parameters"),
                     dataTableOutput("table")
                     
                     
                     
    ),
    tabItem(tabName = "satellites",
            # h2("Satellites"),
            plotlyOutput("satellite1"),
            # plotOutput("satellite2")
    )
    
    ),

        
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

