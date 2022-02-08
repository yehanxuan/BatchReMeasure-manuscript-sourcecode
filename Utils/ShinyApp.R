#### shiny app ###
library(shiny)
#runExample("01_hello")

ui <- fluidPage(
  sliderInput(inputId = "nc1", label = "Choose the control sample size" , value = 50,
              min = 10, max = 200),
  sliderInput(inputId = "nt2", label = "Choose the treatment sample size" , value = 50,
              min = 10, max = 200),
  sliderInput(inputId = "nc2", label = "Choose the re-measure sample size" , value = 5,
              min = 5, max = 100),
  verbatimTextOutput("stats"),
  textInput(inputId = "title", 
            label = "Write a title",
            value = "Histogram of Random Normal Value"),
  actionButton(inputId = "go", label = "Update"),
  plotOutput("hist")
)

server <- function(input, output){
  ## Builds a reactive object (reactive expression)
  # data <- reactive({
  #   rnorm(input$nc1)
  # })  
  # Also we can prevent the graph from updating until we hit the button 
  data <- eventReactive(input$go, {
    rnorm(input$nc1)
  })
  
  output$hist <- renderPlot({
   # hist(rnorm(input$nc1), 
    #     main = isolate({input$title}))
    hist(data())
  })
  # output$stats <- renderPrint({
  #   summary(data())
  # })
  #observeEvent(input$clicks, { print(as.numeric(input$clicks)) })
  observe({ print(as.numeric(input$clicks)) })
}

# eventReactive

shinyApp(ui = ui, server = server)




