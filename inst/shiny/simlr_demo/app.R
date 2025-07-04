library(shiny)
library(ggplot2)
library(dplyr)
library(fmsb)  # for radar plot
library(readr)
library(scales) # for alpha in radar plot

# ---- Function to compute embedding, distance, and percentile ----
computeViewScores <- function(indivDf, nhDf, projMatrix) {
  indivEmbed <- as.matrix(indivDf) %*% as.matrix(projMatrix)
  nhEmbed <- as.matrix(nhDf) %*% as.matrix(projMatrix)
  
  mu <- colMeans(nhEmbed)
  covMat <- cov(nhEmbed)

  covMat <- covMat + diag(1e-6, ncol(covMat))
  
  dM <- mahalanobis(indivEmbed, center = mu, cov = covMat)
  nhDists <- mahalanobis(nhEmbed, center = mu, cov = covMat)  
  
  percentile <- ecdf(nhDists)(dM)
  list(distance = dM, percentile = percentile, embedding = indivEmbed)
}

# ---- Load input data and SiMLR projection matrices ----

nh_list <- list()
for (i in seq.int(5)) {
  nh_list[[i]] <- read.csv(paste0("nh_list_", i, ".csv"))
  nh_list[[i]]$X <- NULL
}

nh_list_ages <- read.csv("nh_list_ages.csv")

categories <- c("demog", "diet", "exposures", "mentalhealth", "physical")
simlr_file_prefix <- "simlr_regression_ica_centerAndScale_np"
resultNHv <- read_simlr_data_frames(simlr_file_prefix, categories)

# ---- UI ----

ui <- fluidPage(
  titlePanel("NHANES SiMLR Atypicality Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("subject_index", "Select Subject Index:", 
                  min = 1, max = nrow(nh_list[[1]]), value = 1, step = 1),
      numericInput("subject_text", "Or enter Subject Index:", 
                   value = 1, min = 1, max = nrow(nh_list[[1]])),
      checkboxInput("show_radar", "Show Radar Plot", TRUE),
      tableOutput("interpretTable")
    ),
mainPanel(
  textOutput("overallScore"),
  fluidRow(
    column(
      width = 6,
      verbatimTextOutput("subject_info"),
    ),
    column(
      width = 6,
      verbatimTextOutput("categoryScores")
    )
  ),
  plotOutput("radarPlot")
)
  )
)

# ---- Server ----

server <- function(input, output, session) {

   observeEvent(input$subject_index, {
     updateNumericInput(session, "subject_text", value = input$subject_index)
   })

   observeEvent(input$subject_text, {
     idx <- input$subject_text
     if (!is.na(idx) && idx >= 1 && idx <= nrow(nh_list[[1]])) {
       updateSliderInput(session, "subject_index", value = idx)
     }
   })

   sampled_subject <- reactive({
     idx <- input$subject_index
     
     indiv_list <- lapply(nh_list, function(df) df[idx, , drop = FALSE])
     
     list(idx = idx, indiv_list = indiv_list)
   })
  
  indiv_scores <- reactive({
    req(sampled_subject())
    indiv_list <- sampled_subject()$indiv_list
    
    cat_scores <- list()
    indiv_embeddings <- c()
    percentiles <- c()
    
    for (i in seq_along(categories)) {
      nh_df <- nh_list[[i]]
      proj_matrix <- resultNHv[[i]]
      
      score <- computeViewScores(indiv_list[[i]], nh_df, proj_matrix)
      cat_scores[[categories[i]]] <- score
      indiv_embeddings <- c(indiv_embeddings, score$embedding)
      percentiles <- c(percentiles, score$percentile)
    }
    
    # Overall combined embedding
    nh_combined <- do.call(cbind, lapply(seq_along(categories), function(i) {
      as.matrix(nh_list[[i]]) %*% as.matrix(resultNHv[[i]])
    }))
    
    mu_comb <- colMeans(nh_combined)
    cov_comb <- cov(nh_combined)
    cov_comb <- cov_comb + diag(1e-6, ncol(cov_comb))
    
    d_M_overall <- mahalanobis(indiv_embeddings, center = mu_comb, cov = cov_comb)
    nh_comb_dists <- mahalanobis(nh_combined, center = mu_comb, cov = cov_comb)
    percentile_overall <- ecdf(nh_comb_dists)(d_M_overall)
    
    list(cat_scores = cat_scores, percentiles = percentiles, 
         overall_distance = d_M_overall, overall_percentile = percentile_overall)
  })
  
  output$subject_info <- renderPrint({
    req(sampled_subject())
    idx <- sampled_subject()$idx
    gender_code <- nh_list[[1]]$riagendr_1[idx]
    gender <- ifelse(gender_code == 1, "Male", "Female")

    
    cat("Sampled Subject Info\n")
    cat("----------------------------\n")
    cat("Index: ", idx, "\n")
    cat("Age: ", nh_list_ages$Age[idx], "\n")
    cat("Gender: ", gender, "\n")
  })

  output$interpretTable <- renderTable({
    data.frame(
      `Percentile Range` = c("0–20%", "40–60%", "80–100%"),
      Interpretation = c(
        "Very typical (close to population mean)",
        "Moderately typical",
        "Highly atypical (far from mean)"
      ),
      check.names = FALSE
    )
  })

  output$radarPlot <- renderPlot({
    req(indiv_scores())
    if (!input$show_radar) return(NULL)
    
    percentiles <- indiv_scores()$percentiles
    overall_score <- round( indiv_scores()$overall_percentile * 100, 2)

    output$categoryScores <- renderPrint({
      req(indiv_scores())
      
      percentiles <- indiv_scores()$percentiles
      names(percentiles) <- categories
      
      cat("Overall atypicality percentile: ", overall_score, "\n\n")
      cat("Per-category atypicality percentiles:\n")
      for (cat_name in categories) {
        score <- round(percentiles[cat_name] * 100, 2)
        cat(paste0(" - ", cat_name, ": ", score, "\n"))
      }
    })    

    # Create radar plot data frame
    max_row <- rep(1, length(categories))
    min_row <- rep(0, length(categories))
    user_row <- as.numeric(percentiles)
    names(user_row) <- categories

    radar_df <- rbind(max_row, min_row, user_row)
    colnames(radar_df) <- categories
    rownames(radar_df) <- c("Max", "Min", "User")

    fmsb::radarchart(as.data.frame(radar_df),
                     axistype = 1,
                     pcol = "blue", pfcol = scales::alpha("blue", 0.4),
                     plwd = 2, cglcol = "grey", cglty = 1, axislabcol = "grey",
                     vlcex = 0.9, title = "Atypicality across categories (percentile)")
  })
  
}

# ---- Run App ----
shinyApp(ui, server)
