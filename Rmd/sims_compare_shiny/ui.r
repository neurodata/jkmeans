library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("j-Clustering"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("K", "K:",
                list("2" = "2", 
                     "5" = "5")),
    selectInput("weight", "Weight:",
                list("Balanced" = "1", 
                     "Imbalanced" = "0")),
    selectInput("sigma", "Sigma:",
                list("0.5" = "0.5", 
                     "1" = "1")),
    selectInput("initial_truth", "Initialization from:",
                list("Truth" = "1", 
                     "K-means" = "0")),
    
    checkboxInput("fixW","Show fixed-weight models (k-means family)",  TRUE),
    checkboxInput("GMM","Show flexible-weight models (GMM family)",  TRUE)
    ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(

    plotOutput("ratePlotMCE"),
    plotOutput("ratePlotRMSE")
  )
))