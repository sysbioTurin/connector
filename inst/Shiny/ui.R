library("shiny")
library("shinydashboard")
library("shinyFiles")
library("shinybusy")
library("shinyjs")
library("shinyWidgets")
library("tools")
library("ggplot2")
library("dashboardthemes")

ui <- dashboardPage(
  
  dashboardHeader(
    title = "CONNECTOR",
    tags$li(class = "dropdown",
            tags$style("#mydropdown{cursor: pointer; background-color: #852fb4 ;
                              color: white;
                              border: 1px solid #852fb4; height: 50px;left:0;}"),
            dropdownButton(inputId = "mydropdown",label = "",
                           right = T,
                           circle = FALSE,
                           icon = icon("download"),
                           tags$li(shinyjs::useShinyjs(),
                                   downloadButton(outputId = 'download',href = "Connector.RData",
                                                  label = "R enviroment",title = "Save the R enviroment",
                                                  style = "cursor: pointer; width: 98%;
                                      text-align: center; vertical-align: middle;
                                      border: 1px solid #9809AF;",
                                                  download = "Connector.RData",
                                                  class="dlButton")),
                           tags$li(shinyjs::useShinyjs(),
                                   downloadButton(outputId = 'reportGen',
                                                  href = "report.html",
                                                  download = "report.html",
                                                  label = "Report",title = "Report generation",
                                                  style = "cursor: pointer; width: 98%;
                                      text-align: center; vertical-align: middle;
                                      border: 1px solid #9809AF;",
                                                  class="dlButton")
                                   )
                           ),
            block = TRUE),
    tags$li(a(onclick = "onclick =window.open('https://github.com/qBioTurin/connector')",
              href = NULL,
              icon("github"),
              title = "GitHub",
              style = "cursor: pointer;"),
            class = "dropdown")
            ),
### Sidebar content

  dashboardSidebar(div(style="overflow-y: scroll"),
    sidebarMenu(
      id = "tabs",
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Loading data", tabName = "Loading", icon = icon("upload")),
      menuItem("Pre processing", tabName = "PreProc", icon = icon("edit")),
      menuItem("Model selection", icon = icon("edit"),
               menuSubItem("Setting p", tabName = "pSelection", icon = icon("dashboard")),
               menuSubItem("Setting h", tabName = "hSelection", icon = icon("dashboard")),
               menuSubItem("Setting G", tabName = "FCM", icon = icon("dashboard") )
               ),
      menuItem("Cluster plot", icon = icon("edit"),
               tabName = "ClPlots"
               )
      )
    ),
## Body content
  dashboardBody(
    ### changing theme
    shinyDashboardThemes(
      #theme = "blue_gradient"
      theme = "onenote"
    ),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$script(src = "message-handler.js")
      ),
    
    tabItems(
      # Loading tab content
      tabItem(tabName = "home",
              h2("Connector "),
              h1(""),
              p("Presented by: ",strong("Quantitative Biology")," group of the University of Turin, ",
                tags$a(href= "onclick =window.open('http://beta.di.unito.it/index.php/english/research/groups/computational-and-mathematical-methodologies-modelling-multi-omics-systems/about')", "(q-Bio Turin)")
              ),
              p("Authors: S. Pernice, R. Sirovich, M. Beccuti and F. Cordero",
                style = "font-family: 'times'; font-si16pt"),
              img(src = "groups.png", height = 140, width = 400, style="display: block; margin-left: auto; margin-right: auto"),
              img(src = "founds.png", height = 70, width = 200, style="display: block; position: absolute; bottom: 5%; right: 5%")
      ),
      tabItem(tabName = "Loading",
              h2("Loading Data"),
              tabBox(
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "tabBoxLoading",width = 12 ,
                tabPanel("From files:", "",
                         fluidRow(
                           p("Two files are requested to run the data analysis exploiting the CONNECTOR package:"),
                           p("  - the excel file, namely ", strong("GrowDataFile", style = "font: bold"),", reporting the growth evolution data, "),
                           p("  - the csv file, namely", strong("AnnotationFile", style = "font: bold"),", containing the annotation information associated with the samples. "),
                           p("Hence, the growth data associated with an experiment must be stored into GrowDataFile as a table with two columns for each sample. The first column, labeled  Time, contains the time points of a sample. The second column, labeled by the sample name, contains  the data volume over the time. Instead, the second file (i.e. AnnotationFile)  stores  the annotated information associated with the samples as a table in csv format so that  number of rows is equals to the total number of samples."),
                           column(width=6,
                                  fileInput("fileGrowth", "GrowDataFile:", placeholder = "Choose a xlsx file",
                                            accept = c(
                                              "application/vnd.ms-excel",
                                              "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet") )
                                  ),
                           column(width=6,
                                  fileInput("fileFeat", "AnnotationFile",  placeholder = "Choose a csv file",
                                            accept = c(".csv") )
                                  ),
                           fluidRow(column(width = 10,offset = 1,verbatimTextOutput("LoadingError1"))),
                           fluidRow(column(width = 10,offset = 1,htmlOutput("SummaryLoading"))),
                           column(width=2,
                                  actionButton( label = "Load", inputId = "LoadingData", icon = icon("file-upload") )
                                 ),
                           column(width=1,
                                  actionButton( label = "Next", inputId = "FromS1toS2", icon = icon("arrow-circle-right") )
                           ),
                               # box(title = "Outputs",width = 12,
                               #     background = "black",collapsible = TRUE,
                               #     textOutput(outputId = "console")
                               # ),
                               add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80)
                         )
                ),
                tabPanel("From ConnectorList.RData:", "",
                         fluidRow(
                           column(width=12,
                                  fileInput("RDataImport1", "Connector List:",
                                            placeholder = "Select a Connector.RData",
                                            width = "100%"
                                            ),
                                  fluidRow(column(width = 10,offset = 1,verbatimTextOutput("LoadingError2"))),
                                  fluidRow(column(width = 10,offset = 1,htmlOutput("SummaryLoadingRdata")))
                               ),
                           column(width=2,
                                  actionButton( label = "Load", inputId = "LoadingConnectorList1", icon = icon("file-upload") )
                                  ),
                           add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80)
                           )
                         ),
                tabPanel("From RData:", "",
                         fluidRow(
                           column(width=12,
                                  p("Select an RDs storing the data.frame with three columns: ID for the identification number of each curve, V for the values and Time for the respective time points"),
                                  fileInput("RDataImportGrowthData", "Growth curves data:",
                                            placeholder = "Select an RDs storing the data.frame for the growth data",
                                            width = "100%"
                                            ),
                                  p("Select an RDs storing the data.frame where the first column must be named ID and it stores the identification numbers contained in the GrowDataFrame. Then it is possible to add a column per feature to consider and to associate to the respective sample (depending on the identification number in the first column. If NULL, then in the ConnectorList only the feature ID will be reported. " ),
                                  fileInput("RDataImportFeatureData", "Feature data:",
                                            placeholder = "Select an RDs storing the data.frame for the feature data",
                                            width = "100%"
                                  ),
                                  fluidRow(column(width = 10,offset = 1,verbatimTextOutput("LoadingError3"))),
                                  fluidRow(column(width = 10,offset = 1,htmlOutput("SummaryLoadingRdata2")))
                           ),
                           column(width=2,
                                  actionButton( label = "Load", inputId = "LoadingConnectorList2", icon = icon("file-upload") )
                           ),
                           add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80)
                         )
                )
                )
              ),
      tabItem(tabName = "PreProc",
              h2("Preprocessing"),
              box(width = 12, title = "Data truncation:",collapsible = TRUE,collapsed = TRUE,
                  sliderInput(inputId = "truncTime", label = h4("Truncation time:"),
                              min = 0, max = 0, value = c(0,0),step = 1),
                  column(width=2,
                         actionButton( label = "Truncate", inputId = "TruncateData", icon = icon("cut") )
                  ),
                  column(width=1,
                         actionButton( label = "Next", inputId = "FromS2toS3", icon = icon("arrow-circle-right") )
                  ),
                  add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80),
                  fluidRow(column(width = 10,offset = 1,htmlOutput("SummaryCutting")))
              ),
              box(width = 12, title = "Data visualization:",collapsible = TRUE,
                  selectInput(inputId = "feature", label = "Features:",choices = "" ),
                    column(width = 7,h4("Growth curves"),
                           plotOutput("visualGrowth")
                          ), 
                    column(width = 5, h4("Time grid density"),
                           plotOutput("visualTimes")
                    ),
                  box(width = 12, title = "Save the plots:",collapsible = TRUE,
                      column(width =12,
                             textInput(inputId = "titleGrowth", label = "Plot title",value = "" )
                             ),
                      column(width = 6,
                             textInput(inputId = "xlabGrowth", label = "x-axis label",value = ""  )
                             ),
                      column(width = 6,
                             textInput(inputId = "ylabGrowth", label = "y-axis label",value = ""  )
                             ),
                      column(width = 2,
                             actionButton(label="Change the plot", inputId = "ChangeGrowth", icon = icon("palette"))
                             ),
                      column(width = 1,
                             downloadButton(label="Download",outputId = "PlotDownloadGrowth",
                                            href = "GrowthCurves.pdf",
                                            download = "GrowthCurves.pdf" ,
                                            title = "Save the GrowthCurves plot",
                                            style = "cursor: pointer;",
                                            class="dlButton"))
                  )
               )
      ),
    tabItem(tabName = "pSelection",
              h2("Setting p:"),h4("Cross-log likelihood"),
              fluidRow(
                box(width = 12,h4("Spline dimension:"),
                    column(width = 6,
                       numericInput("Minp", "Min: ", 2, min = 2)
                    ),
                    column(width = 6,
                           numericInput("Maxp", "Max: ", 3, min = 3)
                    ),
                    actionButton( label = "Run", inputId = "RunForP" ),
                    fluidRow(column(width = 10,offset = 1,verbatimTextOutput("ErrorInP"))),
                    add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80)
                  ),
                box(width = 12, title = "CrossLog-likelihood visualization:",collapsible = TRUE,collapsed = T,
                   plotOutput("visualCLLik"),
                   sliderInput(inputId = "pValue", label = h4("Select p value:"),
                               min = 0,max = 0, value = "",step = 1),
                   column(width=2,
                          actionButton( label = "Select",inputId = "pSelection", icon = icon("check") )
                          ),
                   column(width=1,
                          actionButton( label = "Next", inputId = "FromS3toS4", icon = icon("arrow-circle-right") )
                          ),"",
                   box(width = 12, title = "Save the plots:",collapsible = TRUE,
                       column(width =12,
                              textInput(inputId = "titleCLL", label = "Plot title",value = "" )
                       ),
                       column(width = 6,
                              textInput(inputId = "xlabCLL", label = "x-axis label",value = ""  )
                       ),
                       column(width = 6,
                              textInput(inputId = "ylabCLL", label = "y-axis label",value = ""  )
                       ),
                       column(width = 2,
                              actionButton(label="Change the plot", inputId = "ChangeCLL", icon = icon("palette"))
                       ),
                       column(width = 1,
                              downloadButton(label="Download",outputId = "PlotDownloadCLL",
                                             href = "CLL.pdf",
                                             download = "CLL.pdf" ,
                                             title = "Save the CrossLog-likelihood plot",
                                             style = "cursor: pointer;",
                                             class="dlButton"))
                   )
                   )
            )
    ),
    tabItem(tabName = "hSelection",
            h2("Setting h:"),h4("PCA"),
            fluidRow(
              box(width = 6,
                  title = "Run the PCA:",
                  numericInput("pValueInHpanel", "Value of p: ", 1, min = 1),
                  actionButton( label = "Run", inputId = "RunForH" ),
                  add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80),
                  fluidRow(column(width = 10,offset = 1,verbatimTextOutput("ErrorInH")))
              ),
              box(width = 12, title = "PCA visualization:",collapsible = TRUE,collapsed = T,
                  column(width = 8,
                         plotOutput("visualPCA") 
                         ),
                  column(width = 4,
                         h3("Value of h to select:"),
                         numericInput("hValue", "Value of h: ", 1, min = 1),
                         column(width=2,
                                actionButton( label = "Select",inputId = "hSelection", icon = icon("check") )
                         ),
                         column(width=1,offset = 1,
                                actionButton( label = "Next", inputId = "FromS4toS5", icon = icon("arrow-circle-right") )
                         ),"",
                         box(width = 12, title = "Save the plots:",collapsible = TRUE,
                             column(width =12,
                                    textInput(inputId = "titlePCA", label = "Plot title",value = "" )
                             ),
                             column(width = 6,
                                    textInput(inputId = "xlabPCA", label = "x-axis label",value = ""  )
                             ),
                             column(width = 6,
                                    textInput(inputId = "ylabPCA", label = "y-axis label",value = ""  )
                             ),
                             column(width = 2,
                                    actionButton(label="Change the plot", inputId = "ChangePCA", icon = icon("palette"))
                             ),
                             column(width = 1,offset = 1,
                                    downloadButton(label="Download",outputId = "PlotDownloadPCA",
                                                   href = "PCA.pdf",
                                                   download = "PCA.pdf" ,
                                                   title = "Save the PCA plot",
                                                   style = "cursor: pointer;",
                                                   class="dlButton"))
                         )
                  )
              )
            )
    ),
    tabItem(tabName = "FCM",
            h2("Setting G:"),h4("FCM"),
            verbatimTextOutput("FCMallError"),
            fluidRow(
              box(width = 12,
                  title = "Run the FCM:",
                  column(width = 6,
                    numericInput("pValueInGpanel", "Value of p: ", 0, min = 0),
                    numericInput("hValueInGpanel", "Value of h: ", 0, min = 0),
                    numericInput("Nruns", "Number of runs: ", 100, min = 1)
                  ),
                  column(width = 6,
                    box(width = 12,
                      title = "Number of clusters:",
                      numericInput("MinG", "Min G: ", 1, min = 1),
                      numericInput("MaxG", "Max G: ", 2, min = 2)
                    )
                  ),
                  actionButton( label = "Run", inputId = "RunForG" ),
                  add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80)
                  )
              ),
              tabBox(
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "tabBoxSelectingG",width = 12 ,
                tabPanel("Boxplots:", "",
                         fluidRow(
                           column(width=11,
                                  plotOutput("visualIndexes"),
                                  column(width = 1,
                                         downloadButton(label="Download",outputId = "PlotDownloadBoxPlotIndexes",
                                                        href = "BoxPlotIndexes.pdf",
                                                        download = "BoxPlotIndexes.pdf" ,
                                                        title = "Save the BoxPlot Indexes plot",
                                                        style = "cursor: pointer;",
                                                        class="dlButton")
                                         )
                                  )
                           )
                         ),
                tabPanel("Consensus matrix:", "",
                         fluidRow(
                           column(width=11,
                                  sliderInput(inputId = "GConsMat", label = h4("Number of clusters:"),
                                              min = 0,max = 0, value = "",step = 1),
                                  plotOutput("visualConsMatrix"),
                                  column(width = 1,
                                         downloadButton(label="Download",outputId = "PlotDownloadConsMatrix",
                                                        href = "ConsMatrix.pdf",
                                                        download = "ConsMatrix.pdf" ,
                                                        title = "Save the Consensus Matrix plot",
                                                        style = "cursor: pointer;",
                                                        class="dlButton")
                                         )
                             )
                           )
                         ),
                tabPanel("All fDB indexes:", "",
                         fluidRow(
                           column(width =6,
                                  h4("Tightness indexes:"),
                                  tableOutput(outputId = "table_tight_all")," ",
                                  h4("fDB1 indexes:"),
                                  tableOutput(outputId = "table_fdb1_all")," "
                                  ),
                           column(width =6,
                                  h4("fDB indexes:"),
                                  tableOutput(outputId = "table_fdb_all")," ",
                                  h4("fDB2 indexes:"),
                                  tableOutput(outputId = "table_fdb2_all")," "
                                  )
                           )
                         ),
                tabPanel("G selection:", "",
                         fluidRow(
                           column(width= 11,
                                  numericInput("Gvalue", "Number of clusters: ", 1, min = 1)
                           ),
                           column(width= 2,offset = 9,
                                  actionButton("GExtrap", "Select"),
                                  add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80)
                           ) 
                         )
                ),
                tabPanel("fDB indexes of the selected cluster:", "",
                         fluidRow(
                           column(width= 11,
                                  verbatimTextOutput("IndexesClusterError")
                                  ),
                           box(width = 12, title = "Coefficents",collapsible = T,collapsed = T, 
                               column(width =6,
                                      h4("Erre indexes:"),
                                      tableOutput("table_fdb_erre"),
                                      h4("Esse indexes:"),
                                      tableOutput("table_fdb_esse")
                                   ),
                               column(width =6,
                                      h4("Emme indexes:"),
                                      tableOutput("table_fdb_emme")
                                      )
                               ),
                           box(width = 12, title ="First derivative coefficents",collapsible = T,collapsed = T,
                               column(width =6,
                                      h4("Erre indexes:"),
                                      tableOutput("table_fdb1_erre"),
                                      h4("Esse indexes:"),
                                      tableOutput("table_fdb1_esse")
                                   ),
                               column(width =6,
                                      h4("Emme indexes:"),
                                      tableOutput("table_fdb1_emme")
                                   )
                               ),
                           box(width = 12,
                               title= "Second derivative coefficents",collapsible = T,collapsed = T,
                               column(width =6,
                                      h4("Erre indexes:"),
                                      tableOutput("table_fdb2_erre"),
                                      h4("Esse indexes:"),
                                      tableOutput("table_fdb2_esse")
                                   ),
                               column(width =6,
                                      h4("Emme indexes:"),
                                      tableOutput("table_fdb2_emme")
                                   )
                               )
                           )
                         )
                )
    ),
    tabItem(tabName = "ClPlots",
            h2(""),
            verbatimTextOutput("IndexesClusterErrorCLplot"),
            tabBox(
              id = "tabBoxCLplot",width = 12,
              tabPanel("Clusters", "",
                       fluidRow(
                         column(width= 6,
                                selectInput(inputId = "featureCLplot",
                                            label = "Features:",choices = "ID",
                                            selected = "ID")
                                ),
                         column(width= 6,
                                selectInput(inputId = "typeCLplot",
                                            label = "Plot:",
                                            choices = c("All clusters","Cluster means") ,
                                            selected = "All clusters"  )
                                ),
                         plotOutput("visualClplot"),
                         column(width = 1,
                                downloadButton(label="Download",outputId = "PlotDownloadClplot",
                                               href = "Clplot.pdf",
                                               download = "Clplot.pdf" ,
                                               title = "Save the Cluster plot",
                                               style = "cursor: pointer;",
                                               class="dlButton")
                                )
                         ),
                       add_busy_spinner(spin = "fading-circle",color = "white",height = 80,width=80)
                       ),
              tabPanel("Discriminant", "",
                       fluidRow(
                         column(width = 4,
                                selectInput(inputId = "ColorDiscrPlot",
                                            label = "Curve colors:",choices = c("Clusters", "Features" ),
                                            selected = "Cluster")
                                ),
                         column(width = 4,
                                uiOutput("FeatSelection")
                                ),
                         plotOutput("DiscrPlot"),
                         column(width = 1,
                                downloadButton(label="Download",outputId = "PlotDownloadDiscrPlot",
                                               href = "DiscrPlot.pdf",
                                               download = "DiscrPlot.pdf" ,
                                               title = "Save the Discriminant plot",
                                               style = "cursor: pointer;",
                                               class="dlButton")
                                )
                         )
                       ),
              tabPanel("Spline", "",
                       fluidRow(
                         column(width= 11,
                              selectInput(inputId = "SplineID",
                                   label = "Sample:",
                                   choices = ""),
                              plotOutput("visualSplinePlot"),
                              column(width = 1,
                                     downloadButton(label="Download",outputId = "PlotDownloadSplinePlot",
                                                    href = "SplinePlot.pdf",
                                                    download = "SplinePlot.pdf" ,
                                                    title = "Save the Spline plot",
                                                    style = "cursor: pointer;",
                                                    class="dlButton")
                              )
                              )
                         )
                       ),
              tabPanel("Inspection", "",
                       fluidRow(
                         column(width= 3,
                                selectInput(inputId = "FeatInsp",
                                            label = "Feature:",choices = "ID",selected = "ID")
                                )
                              ),
                       fluidRow(
                         column(width = 6,
                                tableOutput(outputId = "tableClusterNames")
                                ),
                         column(width = 6,
                                tableOutput(outputId = "tableCounting")
                                )
                         )
                       )
              )
    )
    # ,
    # tabItem(tabName = "Plots",
    #         h2("Plots obtained through the anaylsys"),
    #         fluidRow(
    #           box(width = 12, title = "Data Visualization:",collapsible = TRUE,collapsed = T,
    #               column(width = 7,h4("Plots"),
    #                      plotOutput("visualGrowthEND")
    #               ),
    #               column(width = 5,h4("Plot features"),
    #                      #selectInput(inputId = "feature", label = "Features:",choices = "" ),
    #                      textInput(inputId = "title1", label = "Plot title",value = "" ),
    #                      column(width = 6,
    #                             textInput(inputId = "xlab1", label = "x-axis label",value = ""  ),
    #                             actionButton(label="Change", inputId = "Change1")),
    #                      column(width = 6,
    #                             textInput(inputId = "ylab1", label = "y-axis label",value = ""  ),
    #                             downloadButton(label="Download",outputId = "PlotDownload1",
    #                                            href = "GrowthCurves.pdf",
    #                                            download = "GrowthCurves.pdf" ,
    #                                            title = "Save the GrowthCurves plot",
    #                                            style = "cursor: pointer;",
    #                                            class="dlButton"))
    #               )
    #           ),
    #           box(width = 12,
    #               title = "CrossLog-likelihood visualization:",collapsible = TRUE,collapsed = T,
    #               column(width = 7,h4("Plots"),
    #                      plotOutput("visualCLLikEND")
    #               ),
    #               column(width = 5,h4("Plot features"),
    #                      textInput(inputId = "title2", label = "Plot title",value = "" ),
    #                      column(width = 6,
    #                             textInput(inputId = "xlab2", label = "x-axis label",value = ""  ),
    #                             actionButton(label="Change", inputId = "Change2")),
    #                      column(width = 6,
    #                             textInput(inputId = "ylab2", label = "y-axis label",value = ""  ),
    #                             downloadButton(label="Download",outputId = "PlotDownload2",
    #                                            href = "CrossLogLikelihood.pdf",
    #                                            title = "Save the CrossLogLikelihood plot",
    #                                            style = "cursor: pointer;",
    #                                            class="dlButton")
    #                      )
    #               )
    #           )
    #        )
    #)# tabitem
  )# tabitems
  )
)




