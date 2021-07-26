#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Define UI for application that draws a histogram
source('global.R', local = TRUE)
options(shiny.maxRequestSize = 100*1024^2)
# shinyUI(


sidebar <- dashboardSidebar(
    # sidebarMenu(id="tabs", sidebarMenuOutput("menu"))
    sidebarMenu(id="tabs",
                menuItem("Upload Files", tabName = "upload", icon = icon("file-upload", lib = "font-awesome")),
                menuItem("3D View", tabName = "umap", icon = icon("bar-chart-o")),
                menuItem("2D projections", icon = icon("bar-chart-o"), tabName = "satellites"),
                menuItem("Console output", icon = icon("terminal", lib = "font-awesome"), tabName = "console"),
                         # badgeLabel = "new", badgeColor = "green"),
                conditionalPanel(
                  condition = "input.btnanalysis != 0",
                hr(),
                # convertMenuItem(
                # menuItem("Adjust the colors", tabName = "sliders_axes", icon = icon("bar-chart-o"),
                chooseSliderSkin("Square"),
                # chooseSliderSkin("Flat", color = "#112446"),
                # setSliderColor(c("DeepPink ", "#FF4500", "", "Teal"), c(1)),
                sliderInput("weightL",
                            "L* weight (Brightness):",
                            min = 1,
                            max = 3,
                            value = 1,
                            step = 0.1,
                            width = '500px'),
                sliderInput("weightA",
                            "a* weight (Red - Green):",
                            min = 1,
                            max = 3,
                            value = 1,
                            step = 0.1,
                            width = '500px'),
                sliderInput("weightB",
                            "b* weight (Yellow - Blue):",
                            min = 1,
                            max = 3,
                            value = 1,
                            step = 0.1,
                            width = '500px'),
                actionButton("weightButton", "Re-optimize", width = '100px'),
                br(),
                sliderInput("scaling",
                            "Scaling factor:",
                            min = 1,
                            max = 2,
                            value = 1,
                            step = 0.1,
                            width = '500px'),
                
                actionButton("scalingButton", "Re-scale", width = '100px'),
                # ),
                # ,"umap"),
                hr()
                ), #conditionalPanel
                menuItem("Download", icon = icon("file-export", lib = "font-awesome"), tabName = "Download")
    )
)

  
body <- 
    dashboardBody(
        fluidPage(
          tags$script(HTML(
            'if (!document.cookie.includes("UCIEsession=")) document.cookie = "UCIEsession="+parseInt(Math.floor(Math.random()*1024))'
          )),
          
        # fixedPage(
          useShinyjs(),
    # Application title
    # titlePanel("U-CIE"),
    titlePanel(div(HTML("<strong>U-CIE</strong> <em> [/juː 'siː/] </em> "))),
    
    
    tabItems(tabItem(tabName = "upload",
                     # fileInput("file1", "Upload File",
                     #           multiple = FALSE,
                     #           accept = c("text/csv",
                     #                      "text/comma-separated-values,text/plain",
                     #                      ".csv")),

                     br(),
                     
                     conditionalPanel(
                       condition = "input.uiLoadGraphOptionsInput == 'oF'",
                       prettyRadioButtons("matrix", "1: Type of Data",
                                          choices = c('Single-cells' = 'Single-cells',
                                                      'High Dimensional' = 'High Dimensional',
                                                      'Distance matrix'=  'Distance matrix',
                                                      '3D data' = '3D data'), 
                                          selected = character(0), status = 'warning', inline = T),
                     ),
                     br(),
                     
                     # After ISMB
                     # conditionalPanel(
                     #   condition = "input.matrix != null",
                       
                     selectInput("uiLoadGraphOptionsInput",
                       "2: Choose File(s)",
                       c(
                         "File upload" = "oF",
                         "Single-cells example: GSE75748_time_course" = "oR_Example1",
                         # "Single-cells Big example: GSE75748_cell_type" = "oR_Example2",
                         "3D example" = "oR_Example3"
                         )
                     ),
                     
                     uiOutput("uiLoadGraphOptionsOutput"),
                     # ), # conditionPannel
                     conditionalPanel(
                       condition = "input.uiLoadGraphOptionsInput == 'oF'",
                     prettyCheckbox("header", "Header", TRUE, status = 'danger', bigger = T)),
                     
                     tags$head(
                       tags$style(HTML('#btnanalysis{border-color:red}'))
                     ),
                     actionButton("btnanalysis", "Analysis"),

                     # Horizontal line ----
                     tags$hr(),
                     
                     # Input: Checkbox if file has header ----
                     # selectInput("matrix", "Type of Matrix:",
                     #              choices = c('-' = "-",
                     #                          'Single-cells' = 'Single-cells',
                     #                          'High Dimensional' = 'High Dimensional',
                     #                          'Distance matrix'=  'Distance matrix',
                     #                          '3D data' = '3D data'),
                     #              selected = '-'),
                     
                     
                     # br(),
                     # Input: Select separator ----
                     
                     # prettyRadioButtons("sep", "Type of file",
                     #              choices = c("Comma-serepated" = ",",
                     #                          "Semicolon-seperated" = ";",
                     #                          "Tab-seperated" = "\t",
                     #                          "Excel" = "xlsx"),
                     #              selected = "\t"),
                     
                     
                     
                     # Horizontal line ----
                     # tags$hr(),
                     
                     # Input: Select number of rows to display ----
                     # prettyRadioButtons("disp", "Display",
                     #              choices = c(Head = "head",
                     #                          All = "all"),
                     #              selected = "head", status = 'info', inline = T),
                     # tags$hr(),
                     uiOutput("uploaded_dataset"),
                     # box(tableOutput("contents"), width=12,background ="purple"),
                     dataTableOutput("contents")
                     
                     
    ),
    
    tabItem(tabName = "umap",
                     # h2("Parameters and UMAP"),
            br(),
            br(),
            # fluidRow(box(plotlyOutput("plotly_plot"),width=12,title="UMAP with CIE L* a* b* colors", height = "800px", background ="black",collapsible = F),
            #          box(dataTableOutput("table"))),
            box(plotlyOutput("plotly_plot"), title = "UMAP with CIE L* a* b* colors", width = NULL, height = "800px"),
            box(dataTableOutput("table"), title = "Color-options based on optimization algorithm", width = NULL),
            uiOutput("list_of_parameters")), 
            
    
    tabItem(tabName = "console",
            pre(id = "console")),
    # tabItem(tabName = "sliders_axes",
           # fixedRow(
           #  chooseSliderSkin("Square"),
           #  # chooseSliderSkin("Flat", color = "#112446"),
           #  # setSliderColor(c("DeepPink ", "#FF4500", "", "Teal"), c(1)),
           # sliderInput("weightL",
           #                       "L* weight (Brightness):",
           #                       min = 1,
           #                       max = 3,
           #                       value = 1,
           #                       step = 0.1, 
           #                       width = '500px'),
           #           sliderInput("weightA",
           #                       "a* weight (Red - Green):",
           #                       min = 1,
           #                       max = 3,
           #                       value = 1,
           #                       step = 0.1, 
           #                       width = '500px'),
           #           sliderInput("weightB",
           #                       "b* weight (Yellow - Blue):",
           #                       min = 1,
           #                       max = 3,
           #                       value = 1,
           #                       step = 0.1, 
           #                       width = '500px'),
           #           sliderInput("scaling",
           #                       "Scaling factor multiplier:",
           #                       min = 1,
           #                       max = 2,
           #                       value = 1,
           #                       step = 0.1, 
           #                       width = '500px'),
           #           actionButton("weightButton", "Weight", width = '100px'),
           #           actionButton("scalingButton", "Scale", width = '100px'),
           # ),
                   
            # progressBar(id = "pb1", value = 0, size = "xs", display_pct = F, status = 'info', striped = T),
                     
                     

    # ),
    tabItem(tabName = "satellites",
            # h2("Satellites"),
            br(),
            br(),
            column(8,plotlyOutput("satellite"),
                   hidden(actionButton("remove_genes", "Remove Genes and Re-color!")), # Hidden up to ISMB
                   hidden(actionButton("reset_genes", "Reset")), # Hidden up to ISMB
                   # uiOutput("reset")
                   ), 
            
            column(4,dataTableOutput("legend"))
            # ,
            # column(4,dataTableOutput("table2")),


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

