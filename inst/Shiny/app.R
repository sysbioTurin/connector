
Appui <- system.file("Shiny","ui.R", package = "connector")
Appserver <- system.file("Shiny","server.R", package = "connector")

shinyApp(Appui, Appserver, options =  options(shiny.maxRequestSize=500*1024^2, shiny.launch.browser = .rs.invokeShinyWindowExternal) )




