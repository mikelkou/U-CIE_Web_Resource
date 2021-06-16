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
dashboardPage(
    dashboardHeader(title = "Content"),
    dashboardSidebar(
        sidebarMenu(id="tabs", sidebarMenuOutput("menu"))
    ),
    
    dashboardBody(
        fluidPage(
    
    # Application title
    # titlePanel("U-CIE"),
    titlePanel(div(HTML("<strong>U-CIE</strong> <em> [/ˈjuː siː/] </em> "))),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "dataset",
                        label = "Choose a dataset:",
                        choices = c("GSE75748_time_course", "GSE75748_cell_type")),
            
            sliderInput("weightL",
                        "L* weight:",
                        min = 1,
                        max = 3,
                        value = 1,
                        step = 0.1),
            sliderInput("weightA",
                        "a* weight:",
                        min = 1,
                        max = 3,
                        value = 1,
                        step = 0.1),
            sliderInput("weightB",
                        "b* weight:",
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
            actionButton("scalingButton", "Scale")
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            
            plotlyOutput("plotly_plot"),
            uiOutput("list_of_parameters"),
            dataTableOutput("table"),

            plotOutput("satellite1"),
            plotOutput("satellite2")
            # verbatimTextOutput("summary")
            
        )
    )
))
)
