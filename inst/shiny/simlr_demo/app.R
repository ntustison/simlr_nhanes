# app.R - Shiny demo for NHANES SiMLR projection using full dataset

library(shiny)
library(ggplot2)
library(dplyr)
library(simlr.nhanes)  # Correct package name

# Load data
data("nhanescog_2011_2014")
data("nhanes_dict")

df <- nhanescog

# Ensure SiMLR projection exists in the data
if (!all(c("simlr_x", "simlr_y") %in% names(df))) {
  stop("Projection columns simlr_x and simlr_y are missing from nhanescog_2011_2014.")
}

# Filter dictionary to include only present variables
dict_vars <- nhanes_dict %>%
  filter(variable %in% names(df)) %>%
  arrange(label)

# UI
ui <- fluidPage(
  titlePanel("NHANES SiMLR Projection Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("color_var", "Color by:",
                  choices = setNames(dict_vars$variable, dict_vars$label),
                  selected = "age"),
      checkboxInput("show_labels", "Show SEQN (Subject IDs)", FALSE)
    ),
    
    mainPanel(
      plotOutput("simlr_plot", height = "600px")
    )
  )
)

# Server
server <- function(input, output, session) {
  output$simlr_plot <- renderPlot({
    req(input$color_var)
    
    p <- ggplot(df, aes(x = simlr_x, y = simlr_y)) +
      geom_point(aes(color = .data[[input$color_var]]), size = 2, alpha = 0.8) +
      labs(
        x = "SiMLR Dimension 1",
        y = "SiMLR Dimension 2",
        color = input$color_var
      ) +
      theme_minimal()
    
    if (input$show_labels && "SEQN" %in% names(df)) {
      p <- p + geom_text(aes(label = SEQN), size = 3, vjust = -0.7)
    }
    
    p
  })
}

# Run the app
shinyApp(ui = ui, server = server)
