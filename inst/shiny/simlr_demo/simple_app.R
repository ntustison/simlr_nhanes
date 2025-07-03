library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Toggle y = x Line"),
  
  sidebarLayout(
    sidebarPanel(
      actionButton("toggle", "Toggle y = x Line")
    ),
    
    mainPanel(
      plotOutput("line_plot")
    )
  )
)

server <- function(input, output, session) {
  # Reactive value to store the toggle state
  show_line <- reactiveVal(TRUE)
  
  observeEvent(input$toggle, {
    show_line(!show_line())  # Toggle the value
  })
  
  output$line_plot <- renderPlot({
    df <- data.frame(x = seq(-1, 1, length.out = 100))
    
    p <- ggplot(df, aes(x = x)) +
      ylim(-1, 1) +
      xlim(-1, 1) +
      theme_minimal() +
      labs(title = "y = x Line Toggle")
    
    if (show_line()) {
      p <- p + geom_line(aes(y = x), color = "blue", linewidth = 1)
    }
    
    p
  })
}

shinyApp(ui = ui, server = server)
