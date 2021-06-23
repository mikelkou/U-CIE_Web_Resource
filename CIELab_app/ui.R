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
    sidebarMenu(menuItem("Umap", tabName = "umap", icon = icon("dashboard")),
                menuItem("Satellites", icon = icon("th"), tabName = "widgets",
                         badgeLabel = "new", badgeColor = "green")
    )
)

# dashboardPage(skin = "purple",
#     dashboardHeader(title = "U-CIE"),
#     
#     dashboardSidebar(
#         # sidebarMenu(id="tabs", sidebarMenuOutput("menu"))
#         sidebarMenu(menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
#                     menuItem("Widgets", icon = icon("th"), tabName = "widgets",
#                              badgeLabel = "new", badgeColor = "green")
#         )
#     )
#         
#     ,
#     
    
body <- 
    dashboardBody(
        fluidPage(
    
    # Application title
    # titlePanel("U-CIE"),
    titlePanel(div(HTML("<strong>U-CIE</strong> <em> [/juː 'siː/] </em> "))),

    tabItems(tabItem(tabName = "umap",
                     # h2("Parameters and UMAP"),
                     selectInput(inputId = "dataset",
                                 label = "Choose a dataset:",
                                 choices = c("GSE75748_time_course", "GSE75748_cell_type")),
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
    tabItem(tabName = "widgets",
            # h2("Satellites"),
            plotOutput("satellite1"),
            plotOutput("satellite2")
    )),

        
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
    dashboardHeader(title = "U-CIE"),
    sidebar,
    body
)

