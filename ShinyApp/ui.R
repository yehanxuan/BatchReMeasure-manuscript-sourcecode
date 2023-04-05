rm(list = ls())
library(shiny)
library(shinythemes)
library(ggplot2)
library(reshape)
library(cowplot)
library(shinyscreenshot)

shinyUI(
  fluidPage(
    titlePanel("Power Calculation"),
    sidebarPanel(
      numericInput(inputId = "nc1", label = "Control sample size" , value = 50,
                   min = 10, max = 200, step = 1),
      numericInput(inputId = "nt2", label = "Case sample size" , value = 50,
                   min = 10, max = 200, step = 1),
      numericInput(inputId = "r2", "Inter-batch correlation", 0.3, min = -0.95, max = 0.95, step = 0.1 ),
      numericInput(inputId = "a0", "Cohen\'s d", 0.5, min = 0, max = 5, step = 0.1),
      numericInput(inputId = "alpha", "Significant level", 0.05, min = 0.02, max = 0.1, step = 0.01),
      actionButton("goButton", "Go!"),
      p("Click the button to update the plot displayed in the main panel.")
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Absolute Power", plotOutput("plot1")),
                  tabPanel("Relative Power", plotOutput("plot2"))
      )
    ),
    
    
    downloadButton("download", "Download the plot"),
    screenshotButton( label = "Capture entire page", scale = 2),
  )
)