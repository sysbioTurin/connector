#' Shiny app
#'
#' @description
#' function to lunch the shiny application.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @import shiny shinydashboard shinyFiles  shinybusy shinyjs tools ggplot2
#' @export
#' 
RunConnectorShiny <- function() {
  x = T
  shiny::runApp(
    appDir = system.file("Shiny", package = "connector")
                )
}

