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
    # put sidebarPanel and mainPanel in the fluidRow to make alignment 
    fluidRow( 
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
               ),
               downloadButton("download", "Download the plot"), 
               screenshotButton( label = "Capture entire page", scale = 2)
             )
    ),
    
    fluidRow(
      
      column(
        width = 12, 
        h3("Definitions"),
        p(strong("Control sample size:"), "the size of control samples in the batch 1."),
        p( strong("Case sample size:"), "the size of case samples in the batch 2."),
        p( strong("Inter-batch correlation:"), "the correlation between the remeasured samples in batch 1 and 2."),
        p(strong("Cohen's d:"), "the standardized effect size according to Cohen's criterion."),
        p(strong("Significant level:"), " the probability of rejecting the null hypothesis, given that it is true."), 
    #  ),
    #  column(
    #    width = 6,
    #    h3("Power Definitions"),
        p(strong("Absolute power:"), "the probability of correctly rejecting the null hypothesis when it is false."),
    #    p(strong("Relative power:"), "the ratio of the absolute power for two different sample sizes or effect sizes.")
        p(strong("Relative power:"), "The ratio of the absolute power calculated based on the remeasured sample size to the power calculated using all remeasured samples.")
      )
    )
    
  )
  
)