#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
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

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("CIELab"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "dataset",
                        label = "Choose a dataset:",
                        choices = c("GSE75748_time_course", "GSE75748_cell_type")),
            
            sliderInput("weightL",
                        "How big should the L* be:",
                        min = 0.5,
                        max = 3,
                        value = 1,
                        step = 0.1),
            sliderInput("weightA",
                        "How big should the a* be:",
                        min = 0.5,
                        max = 3,
                        value = 1,
                        step = 0.1),
            sliderInput("weightB",
                        "How big should the b* be:",
                        min = 0.5,
                        max = 3,
                        value = 1,
                        step = 0.1),
            sliderInput("scaling",
                        "Increase scaling factor by x times:",
                        min = 0.4,
                        max = 5,
                        value = 1,
                        step = 0.2)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotlyOutput("plotly_plot"),
            plotOutput("satellite1"),
            plotOutput("satellite2")
            # verbatimTextOutput("summary")
            
        )
    )
))
