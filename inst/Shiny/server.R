
UpdatingFCMall<-function(S.cl,input, output, session,rResult,Flags)
{
  output$table_tight_all<- renderTable({ 
    Tight = S.cl$BoxPlots[[paste("h=",input$hValueInGpanel)]]$Data$Tight
    OutTable<-data.frame(Min = apply(Tight,1,"min"),
                         Mean = apply(Tight,1,"mean"),
                         Median = apply(Tight,1,"median"),
                         Max = apply(Tight,1,"max") )
    rownames(OutTable) <- rownames(Tight)
    OutTable
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  output$table_fdb_all <- renderTable({ 
    fDB = S.cl$BoxPlots[[paste("h=",input$hValueInGpanel)]]$Data$fDB
    OutTable<-data.frame(Min = apply(fDB,1,"min"),
                         Mean = apply(fDB,1,"mean"),
                         Median = apply(fDB,1,"median"),
                         Max = apply(fDB,1,"max") )
    rownames(OutTable) <- rownames(fDB)
    OutTable
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  output$table_fdb1_all<- renderTable({ 
    fDB.1 = S.cl$BoxPlots[[paste("h=",input$hValueInGpanel)]]$Data$fDB.1
    OutTable<-data.frame(Min = apply(fDB.1,1,"min"),
                         Mean = apply(fDB.1,1,"mean"),
                         Median = apply(fDB.1,1,"median"),
                         Max = apply(fDB.1,1,"max") )
    rownames(OutTable) <- rownames(fDB.1)
    OutTable
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  output$table_fdb2_all<- renderTable({ 
    fDB.2 = S.cl$BoxPlots[[paste("h=",input$hValueInGpanel)]]$Data$fDB.2
    OutTable<-data.frame(Min = apply(fDB.2,1,"min"),
                         Mean = apply(fDB.2,1,"mean"),
                         Median = apply(fDB.2,1,"median"),
                         Max = apply(fDB.2,1,"max") )
    rownames(OutTable) <- rownames(fDB.2)
    OutTable
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  
  IndexBoxPlot<- BoxPlot.Extrapolation(stability.list = S.cl, h = input$hValueInGpanel)
  
  rResult$IndexBoxPlot <- IndexBoxPlot
  rResult$CONNECTORListFCM_all <- S.cl
  output$visualIndexes <- renderPlot(IndexBoxPlot)
  
  
  
  ## Selecting a number of clusters:
  observeEvent(input$GConsMat,{
    ConsMatrix <- ConsMatrix.Extrapolation(stability.list = S.cl,
                                          h = input$hValueInGpanel,
                                          G = input$GConsMat)
    output$visualConsMatrix <- renderPlot(ConsMatrix)
    rResult$ConsMatrix <- ConsMatrix
  })
  
  observeEvent(input$GExtrap,{
    if(!is.null(rResult$CONNECTORList.FCM))
    {
      session$sendCustomMessage(type = 'testmessage', 
                                message = 'The number of clusters choise has been updated!') 
    }
    
    CONNECTORList.FCM <- MostProbableClustering.Extrapolation(stability.list = S.cl,
                                                              h = input$hValueInGpanel,
                                                              G = input$Gvalue )
    rResult$CONNECTORList.FCM  <- CONNECTORList.FCM 
    rResult$G <- input$Gvalue
  })
  
  Flags$Continue <- FALSE
  
  output$IndexesClusterError <- renderText({
      validate(need(!is.null(rResult$CONNECTORList.FCM), "Please select a number of clusters!!!") )
 
    Flags$Continue = TRUE
    "Done!"
  })
  observeEvent(Flags$Continue, {
    if(Flags$Continue)
    {
      ConnFCM <- rResult$CONNECTORList.FCM
      output$table_fdb_erre <- renderTable({
        ConnFCM$Cl.Info$Coefficents$erre
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb_esse <- renderTable({
        ConnFCM$Cl.Info$Coefficents$esse
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb_emme <- renderTable({
        ConnFCM$Cl.Info$Coefficents$emme
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb1_erre <- renderTable({
        ConnFCM$Cl.Info$Deriv.Coefficents$erre
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb1_esse <- renderTable({
        ConnFCM$Cl.Info$Deriv.Coefficents$esse
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb1_emme <- renderTable({
        ConnFCM$Cl.Info$Deriv.Coefficents$emme
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb2_erre <- renderTable({
        ConnFCM$Cl.Info$Deriv2.Coefficents$erre
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb2_esse <- renderTable({
        ConnFCM$Cl.Info$Deriv2.Coefficents$esse
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb2_emme <- renderTable({
        ConnFCM$Cl.Info$Deriv2.Coefficents$emme
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
    }
  })
  
}

PlotModify <- function(plot,title,x,y,feature = NULL)
{
    if(is.null(feature)){
      plot = plot + labs( x = x, y= y, title = title)
    }else{
      plot = plot + labs( x = x, y= y, title = title)
    }
    
    return(plot)
}

withConsoleRedirect <- function(expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr )
  if (length(txt) > 0) {
  #   # insertUI(paste0("#", containerId), where = "beforeEnd",
  #   #          ui = paste0(txt, "\n", collapse = "") )
     return(txt)
  }
  return("")
}

Visual_1Step <- function(CONNECTORList,input, output, session,rResult,Flags)
{
  GrowPlot<-GrowthCurve(CONNECTORList, feature = "ID")
  TimeGrid <- TimeGridDensity(CONNECTORList)
  OutGrowPlot <- GrowPlot$GrowthCurve_plot
  
  output$visualGrowth <- renderPlot({OutGrowPlot})
  output$visualGrowthEND <- renderPlot({OutGrowPlot})
  rResult$GrowPlot <- OutGrowPlot
  
  output$visualTimes <- renderPlot({TimeGrid})
  # updating the feature to select:
  Feat = colnames(CONNECTORList$LabCurv)
  updateSelectInput(session, "feature",
                    choices = Feat )
  updateSliderInput(session,"truncTime",
                    min = min(CONNECTORList$TimeGrid),
                    max = max(CONNECTORList$TimeGrid),
                    value = c(min(CONNECTORList$TimeGrid), max(CONNECTORList$TimeGrid) ) )
  
  ## Growth plot:  
  observeEvent(c(input$feature, input$TruncateData),{
    MinTrunc<-input$truncTime[1]
    MaxTrunc<-input$truncTime[2]
    if(MaxTrunc != 0 &
       MinTrunc != 0 &
       (MaxTrunc != max(CONNECTORList$TimeGrid) || MinTrunc != min(CONNECTORList$TimeGrid)) )
    {
      if(input$feature == "" ){
        txt <- withConsoleRedirect(
          CONNECTORList_trunc<- DataTruncation(CONNECTORList, feature="ID",
                                               truncTime = input$truncTime )
        )
        output$SummaryCutting <- renderUI({HTML(paste(txt[-1], collapse = "<br/>"))} ) 
        
      }else{
        txt <- withConsoleRedirect(
          CONNECTORList_trunc<- DataTruncation(CONNECTORList,
                                               feature= input$feature,
                                               truncTime = input$truncTime)
        )
        output$SummaryCutting <- renderUI({HTML(paste(txt[-1], collapse = "<br/>"))} ) 
        
      }
      
      rResult$CONNECTORList_trunc <- CONNECTORList_trunc
      
      OutGrowPlot <- CONNECTORList_trunc$GrowthCurve_plot
      output$visualGrowth <- renderPlot({OutGrowPlot})
      rResult$GrowPlot <- OutGrowPlot
      
      TimeGridOut <- TimeGrid +
        geom_vline(xintercept = input$truncTime) +
        geom_hline(yintercept = input$truncTime)
      
      output$visualTimes <- renderPlot({TimeGridOut})
    }else{
      if(input$feature != "" ){
        GrowPlot <-GrowthCurve(CONNECTORList, feature = input$feature)
        OutGrowPlot <- GrowPlot$GrowthCurve_plot
        rResult$GrowPlot <- OutGrowPlot
        output$visualGrowth <- renderPlot({OutGrowPlot})
      }
    } 
    
    output$visualGrowth <- renderPlot({OutGrowPlot })
  })
  
  ## Now it is possible to save something:
  Flags$save1 = 0
  Flags$GoS1toS2 = 0
}

server <- function(input, output, session) {
  options(shiny.maxRequestSize=500*1024^2)
  
  CONNECTORList<-NULL
  CONNECTORList_trunc<-NULL
  
  # DECLARE REACTIVEVALUES FUNCTION HERE
  rResult <- reactiveValues(CONNECTORList_trunc = CONNECTORList_trunc,
                            CONNECTORList = CONNECTORList,
                            p = NULL,
                            h = NULL,
                            G = NULL,
                            GrowPlot=NULL,
                            TimeGrid=NULL,
                            CrossLL = NULL,
                            PCAplot = NULL,
                            IndexBoxPlot = NULL,
                            FCMplots = NULL,
                            DiscriminantPlot = NULL,
                            DiscriminantPlotF = NULL,
                            ConsMatrix = NULL
                            )
  
  Flags  <- reactiveValues(save1 = 1, GoS1toS2 = 1, GoS3toS4 = 1, GoS4toS5 = 1)
## Data import:

  Flags$ContinueLoadingData <- FALSE
  Flags$ContinueEstimP<- FALSE
  Flags$ContinueEstimH<- FALSE
  Flags$ContinueEstimG<- FALSE
  
  observeEvent(input$LoadingData, {
      output$LoadingError1 <- renderText({
    validate(
      need(!is.null(input$fileGrowth) & !is.null(input$fileFeat) , 
           if(is.null(input$fileGrowth)) "Please select a xlsx file!!"
           else "Please select a csv/txt file!!" )
    )
    
    Flags$ContinueLoadingData <- TRUE
    ""
  })
    
    if(Flags$ContinueLoadingData)
    {
      txt <- withConsoleRedirect(
        CONNECTORList <- DataImport(input$fileGrowth$datapath,input$fileFeat$datapath)
      )
      rResult$CONNECTORList <- CONNECTORList
      output$SummaryLoading <- renderUI({HTML(paste(txt[-1], collapse = "<br/>"))} ) 
      
      Visual_1Step(CONNECTORList,input, output, session,rResult,Flags)
      }
  })

## Run CLLik:
  observeEvent(input$Minp, {
    if(input$Minp >= input$Maxp) updateNumericInput(session,"Maxp",value = input$Minp+1, min = input$Minp+1 )
    
  }) 
  observeEvent(input$RunForP, {
    output$ErrorInP <- renderText({
      validate(
        need(!is.null(rResult$CONNECTORList_trunc) || !is.null(rResult$CONNECTORList) , 
             "Please upload the RData with a ConnectorList object or run the preprocessing step!" )
      )
      Flags$ContinueEstimP<- TRUE
      ""
    })
        
    if(Flags$ContinueEstimP)
    { 
      if(!is.null(rResult$CONNECTORList_trunc)) CONNECTORList_now <- rResult$CONNECTORList_trunc
      else if(!is.null(rResult$CONNECTORList)) CONNECTORList_now <- rResult$CONNECTORList
      
      CrossLogLike<-BasisDimension.Choice(CONNECTORList_now,input$Minp:input$Maxp)
      CrossLogLikePlot<-CrossLogLike$CrossLogLikePlot
      rResult$CrossLL <- CrossLogLikePlot
      
      updateSliderInput(session,"pValue",
                        min = input$Minp,
                        max = input$Maxp,
                        value = input$Minp)
      
      observeEvent(input$pValue,{
        if(input$pSelection == input$Minp)
          CrossLogLikePlot_out <- CrossLogLikePlot
        else CrossLogLikePlot_out <- CrossLogLikePlot + geom_vline(xintercept = input$pValue)
        
        output$visualCLLik <- renderPlot({CrossLogLikePlot_out})
        
        observeEvent(input$pSelection,{
          rResult$p <- input$pValue
          Flags$GoS3toS4 <- 0
          
          mess <- paste("The value of p selected is",input$pValue)
          output$ErrorInP <- renderText({mess})
        })
        
        updateNumericInput(session,"pValueInHpanel",value = input$pValue, min = input$Minp, max = input$Maxp)
        updateNumericInput(session,"pValueInGpanel",value = input$pValue, min = input$Minp, max = input$Maxp)
      })
    }
    
  })

## Run PCA:
  observeEvent(input$pValueInHpanel,{
    if(!is.null(rResult$p) & input$pValueInHpanel != 0)
      if(input$pValueInHpanel != rResult$p)
        session$sendCustomMessage(type = 'testmessage', message = "The value of p selected is different from the one selected in the CrossLog Likelihood step!") 
  })
  observeEvent(input$RunForH, {
    output$ErrorInH <- renderText({
      validate(
        need(!is.null(rResult$CONNECTORList_trunc) || !is.null(rResult$CONNECTORList) , 
             "Please upload the RData with a ConnectorList object or run the preprocessing step!" )
      )
      Flags$ContinueEstimH <- TRUE
      ""
    })
    
    if(Flags$ContinueEstimH)
    {
      if(!is.null(rResult$CONNECTORList_trunc)) CONNECTORList_now <- rResult$CONNECTORList_trunc
      else if(!is.null(rResult$CONNECTORList)) CONNECTORList_now <- rResult$CONNECTORList
      
      pca <- PCA.Analysis(CONNECTORList_now ,p = input$pValueInHpanel)
      
      PCAplot<- pca$plot
      rResult$PCAplot <- PCAplot
      
      output$visualPCA<- renderPlot(PCAplot)
      
      observeEvent(input$hSelection,{
        rResult$h <- input$hValue
        mess <- paste("The value of h selected is",input$hValue)
        output$ErrorInP <- renderText({mess})
        Flags$GoS4toS5 <- 0
        updateNumericInput(session,"hValueInGpanel",value = input$hValue, min = 1)
      })
    }

  })
  
## Run FCM:
  observeEvent(c(input$pValueInGpanel,input$hValueInGpanel),{
    if(!is.null(rResult$p) & input$pValueInGpanel != 0 )
      if(input$pValueInGpanel != rResult$p)
        session$sendCustomMessage(type = 'testmessage', message = "The value of p selected is different from the one selected in the CrossLog Likelihood step!") 
    
    if(!is.null(rResult$h) & input$hValueInGpanel != 0 )
      if(input$hValueInGpanel != rResult$h)
        session$sendCustomMessage(type = 'testmessage', message = "The value of h selected is different from the one selected in the PCA step!") 
    
  })
  observeEvent(input$MinG, {
    if(input$MinG >= input$MaxG)
      updateNumericInput(session,"MaxG",value = input$MinG, min = input$MinG+1 )
  }) 
  observeEvent(input$RunForG, {
    output$FCMallError <- renderText({
      validate(need(!(is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList_trunc)) ,
                    "Please upload the RData with a ConnectorList object or run the preprocessing step!") )
      Flags$ContinueEstimG <- T
      ""
      })
    
    if(Flags$ContinueEstimG)
    {
      if(!is.null(rResult$CONNECTORList_trunc)) CONNECTORList_now <- rResult$CONNECTORList_trunc
      else if(!is.null(rResult$CONNECTORList)) CONNECTORList_now <- rResult$CONNECTORList
      
      updateSliderInput(session, "GConsMat",
                        min =  input$MinG,max = input$MaxG)
      updateNumericInput(session, "Gvalue",
                         value =  input$MinG, min = input$MinG,max = input$MaxG)
      
      S.cl <-StabilityAnalysis(CONNECTORList_now,
                               G = input$MinG:input$MaxG ,
                               h = input$hValueInGpanel,
                               p = input$pValueInGpanel,
                               runs=input$Nruns)
      
      UpdatingFCMall(S.cl,input, output, session,rResult,Flags)
    }
  })

## Run CLuster Visualization
  Flags$ContinueCLplot = FALSE
  
  observeEvent(rResult$CONNECTORList.FCM,{
    
    if(!is.null(rResult$CONNECTORList_trunc)) {
      CONNECTORList_now <- rResult$CONNECTORList_trunc
    }else{
      CONNECTORList_now <- rResult$CONNECTORList
    }
          output$IndexesClusterErrorCLplot <- renderText({
            validate(need(!is.null(rResult$CONNECTORList.FCM) ,
                          if((is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList_trunc)) ){
                            "Please upload the RData with a ConnectorList object or run the preprocessing step!"
                      }else{
                        "Please run the FCM and select a number of clusters!!!"
                        }
                      ) )
            Flags$ContinueCLplot = TRUE

            clNames <- rResult$CONNECTORList.FCM$FCM$cluster$cluster.names
            Feat <- names(CONNECTORList_now$LabCurv)
    
    updateSelectInput(session, "featureCLplot",
                      choices = Feat,
                      selected = "ID" )
    updateSelectInput(session, "FeatInsp",
                      choices = Feat,
                      selected = "ID" ) 
    updateSelectInput(session, "typeCLplot",
                      choices = c("All clusters","Cluster means",paste(clNames,"Cluster") ),
                      selected = "All clusters" )
    
    ""
  })
})
  observeEvent(c(input$ColorDiscrPlot,Flags$ContinueCLplot),{
    if(!is.null(rResult$CONNECTORList_trunc)) {
      CONNECTORList_now <- rResult$CONNECTORList_trunc
    }else{
      CONNECTORList_now <- rResult$CONNECTORList
    }
    
    if(input$ColorDiscrPlot == "Features" & Flags$ContinueCLplot){
      output$FeatSelection <- renderUI({
        selectInput("FeatDiscrCol", "Features:",
                    choices =  names(CONNECTORList_now$LabCurv) )
      })
    }
  })
  observeEvent(c(input$featureCLplot,input$typeCLplot,Flags$ContinueCLplot),{
    if(Flags$ContinueCLplot)
    {
      
      if(!is.null(rResult$CONNECTORList_trunc)) {
        CONNECTORList_now <- rResult$CONNECTORList_trunc
      }else{
        CONNECTORList_now <- rResult$CONNECTORList
      }
      FCMplots<- ClusterWithMeanCurve(clusterdata = rResult$CONNECTORList.FCM,
                                        data= CONNECTORList_now,
                                        feature = input$featureCLplot)
      rResult$FCMplots <- FCMplots
      
      updateSelectInput(session, "SplineID",
                        choices = names(FCMplots$spline.plots),
                        selected = names(FCMplots$spline.plots)[1]
                        )

      if(input$typeCLplot == "All clusters")
      {
        visualClplot <- FCMplots$plotsCluster$ALL
      }else if(input$typeCLplot == "Cluster means"){
        visualClplot <- FCMplots$plotMeanCurve
      }else{
        visualClplot <- FCMplots$plotsCluster[[input$typeCLplot]]
      }
      output$visualClplot <- renderPlot({visualClplot})
      rResult$visualClplot<- visualClplot
      
      observeEvent(input$SplineID,{
        output$visualSplinePlot <- renderPlot({FCMplots$spline.plots[[input$SplineID]]})
        rResult$visualSplinePlot <- FCMplots$spline.plots[[input$SplineID]]
      })
        
    }
    
    })
  observeEvent(c(input$ColorDiscrPlot,input$FeatDiscrCol,Flags$ContinueCLplot),{
    if(Flags$ContinueCLplot)
    {
      if(!is.null(rResult$CONNECTORList_trunc)) {
        CONNECTORList_now <- rResult$CONNECTORList_trunc
      }else{
        CONNECTORList_now <- rResult$CONNECTORList
      }
      output$DiscrPlot <- renderPlot({
        validate(need( rResult$h <= 2,   "h value must be lower or equal to 2."   ) )
        if(input$ColorDiscrPlot == "Clusters"){
          DiscrPlot <- DiscriminantPlot(clusterdata = rResult$CONNECTORList.FCM,
                                        data = CONNECTORList_now,
                                        h = rResult$h,
                                        feature = "ID")
          visualDiscrPlot <- DiscrPlot$ColCluster
          
        }else{
          DiscrPlot <- DiscriminantPlot(clusterdata = rResult$CONNECTORList.FCM,
                                        data = CONNECTORList_now,
                                        h = rResult$h,
                                        feature = input$FeatDiscrCol)
          visualDiscrPlot <- DiscrPlot$ColFeature

        }
        output$DiscrPlot <- renderPlot({visualDiscrPlot})
        rResult$visualDiscrPlot<- visualDiscrPlot
        rResult$DiscriminantPlotF <- DiscrPlot$ColFeature
        rResult$DiscriminantPlot <- DiscrPlot$ColCluster
      })
      
    }
  })
  observeEvent(c(input$FeatInsp,Flags$ContinueCLplot),{
    if(Flags$ContinueCLplot){
      
      if(!is.null(rResult$CONNECTORList_trunc)) {
        CONNECTORList_now <- rResult$CONNECTORList_trunc
      }else{
        CONNECTORList_now <- rResult$CONNECTORList
      }
      
      NumberSamples<-CountingSamples(clusterdata= rResult$CONNECTORList.FCM,
                                     data = CONNECTORList_now,
                                     feature = input$FeatInsp )
      output$tableCounting<- renderTable({ 
        NumberSamples$Counting
        },bordered = TRUE, colnames = TRUE)
      output$tableClusterNames<- renderTable({ 
        NumberSamples$ClusterNames
      },bordered = TRUE, colnames = TRUE) 
      
    }

  })
  
## Costume Plots
# Data Visual
  observeEvent(input$ChangeGrowth,{
    pl<-rResult$GrowPlot
    output$visualGrowthEND <- renderPlot( PlotModify(pl,input$titleGrowth,input$xlabGrowth,input$ylabGrowth) )
  }) 
  output$PlotDownloadGrowth <- downloadHandler(
    filename = function() { paste("GrowthCurve.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$GrowPlot , device = "pdf",width = 10,height = 10,units = "cm")
    }
  )
# Cross Log Likelihood
  observeEvent(input$ChangeCLL,{
    pl<-rResult$CrossLL
    output$visualCLLikEND <- renderPlot( PlotModify(pl,input$titleCLL,input$xlabCLL,input$ylabCLL) )
  }) 
  output$PlotDownloadCLL <- downloadHandler(
    filename = function() { paste("CrossLogLike.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$CrossLL , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# PCA
  observeEvent(input$ChangePCA,{
    pl<-rResult$PCAplot
    output$visualCLLikEND <- renderPlot( PlotModify(pl,input$titlePCA,input$xlabPCA,input$ylabPCA) )
  }) 
  output$PlotDownloadPCA <- downloadHandler(
    filename = function() { paste("PCA.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$PCAplot , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# Cons Matrix
  output$PlotDownloadBoxPlotIndexes <- downloadHandler(
    filename = function() { paste("BoxPlotIndexes.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$BoxPlotIndexes , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# Cons Matrix
  output$PlotDownloadConsMatrix <- downloadHandler(
    filename = function() { paste("ConsMatrix.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$ConsMatrix , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# Cluster plot
  output$PlotDownloadClplot <- downloadHandler(
    filename = function() { paste("Clplot.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$visualClplot , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# DiscrPlot
  output$PlotDownloadDiscrPlot <- downloadHandler(
    filename = function() { paste("DiscrPlot.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$visualDiscrPlot , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# spline plot
  output$PlotDownloadSplinePlot <- downloadHandler(
    filename = function() { paste("SplinePlot.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$visualSplinePlot , device = "pdf",width = 10,height = 8,units = "cm")
    }
  ) 
## Change panel
  observe({
    if(Flags$GoS1toS2 == 0){
      shinyjs::enable("FromS1toS2")
    }else{
      shinyjs::disable("FromS1toS2")
    }
  })
  observe({
    if(Flags$GoS3toS4 == 0){
      shinyjs::enable("FromS3toS4")
    }else{
      shinyjs::disable("FromS3toS4")
    }
  })
  observe({
    if(Flags$GoS4toS5 == 0){
      shinyjs::enable("FromS4toS5")
    }else{
      shinyjs::disable("FromS4toS5")
    }
  })
  
  observeEvent(input$FromS1toS2, {
    updateTabsetPanel(session, "tabs",
                      selected = "PreProc"
    )
  })
  observeEvent(input$FromS2toS3, {
    updateTabsetPanel(session, "tabs",
                      selected = "pSelection"
    )
  })
  observeEvent(input$FromS3toS4, {
    updateTabsetPanel(session, "tabs",
                      selected = "hSelection"
    )
  })
  observeEvent(input$FromS4toS5, {
    updateTabsetPanel(session, "tabs",
                      selected = "FCM"
    )
  })
  
## Upload/ Save / Download:
  observeEvent(input$LoadingConnectorList1, {
    
    Flags$ContinueUploading <- FALSE
    
    output$LoadingError2 <- renderText({
      validate(
        need(!is.null(input$RDataImport1$datapath), "Please select a Connector.RData!!!")
      )
      Flags$ContinueUploading = TRUE
      ""
    })
    
    observeEvent(Flags$ContinueUploading, {

      if(Flags$ContinueUploading)
      {
        load(input$RDataImport1$datapath)
        if(is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList) )
        {
          load  = TRUE
        }else{ ### alert!!! if it is aleready present! 
          showModal(modalDialog(
            title = "Important message",
            "Do you want to update the ConnectorList already present?",
            easyClose = TRUE,
            footer= tagList(actionButton("confirmUpload", "Update"),
                            modalButton("Cancel")
                            )
          )) 
          
        }
        
        observeEvent(input$confirmUpload, {
          removeModal()
          mess = "Value uploaded:  "
          
          rResult$CONNECTORList <- Info$CONNECTORList
          rResult$CONNECTORList_trunc <- Info$CONNECTORList_trunc

          if(is.null(Info$p)){  
            mess <- c(mess,paste(" -No p value; "))
          }else{ 
            mess <- c(mess,paste(" -p = ",Info$p,"; "))
            rResult$p <- Info$p
            updateNumericInput(session,"pValueInHpanel",value = Info$p, min =1 )
            updateNumericInput(session,"pValueInGpanel",value = Info$p, min =1 )
          }
          
          if(is.null(Info$h)){
            mess <- c(mess,paste(" -No h value; "))
          }else{
            mess <- c(mess,paste(" -h = ",Info$h,"; "))
            rResult$h <- Info$h
            updateNumericInput(session,"hValueInGpanel",value = Info$h, min =1 )
          }
          
          if(is.null(Info$G)){
            mess <- c(mess,paste(" -No G value. "))
          }else{
            mess <- c(mess,paste(" -G = ",Info$G,". "))
            rResult$G <- Info$G
          }
          if(!is.null(Info$CONNECTORListFCM_all)){
            rResult$CONNECTORListFCM_all <- Info$CONNECTORListFCM_all
            UpdatingFCMall(rResult$CONNECTORListFCM_all,input, output, session,rResult,Flags)
            Cl<-names(Info$CONNECTORListFCM_all$ConsensusInfo[[1]])
            Cl.number<-as.numeric(gsub("G= ", "", Cl))
            
            updateSliderInput(session, "GConsMat",
                              min =  min(Cl.number),max = max(Cl.number))
            updateNumericInput(session, "Gvalue",
                               value =  min(Cl.number), min =  min(Cl.number),max = max(Cl.number))
            updateNumericInput(session,"MaxG",value =  max(Cl.number) )
            updateNumericInput(session,"MinG",value =  min(Cl.number) )
          }
          if(!is.null(Info$CONNECTORList.FCM)){
            rResult$CONNECTORList.FCM <- Info$CONNECTORList.FCM
          }
          
          if(!is.null(Info$GrowPlot))
          {
            output$visualGrowth <- renderPlot({Info$GrowPlot})
            output$visualGrowthEND <- renderPlot({Info$GrowPlot})
          }
          if(!is.null(Info$TimeGrid))
          {
            output$visualTimes <- renderPlot({Info$TimeGrid})
          }

          output$SummaryLoadingRdata <- renderUI({HTML(paste(mess, collapse = "<br/>"))} ) 
        })
        
      }
    })
    
    ## Now it is possible to save something:
    Flags$save1 = 0
  }) 
  observeEvent(input$LoadingConnectorList2, {
    ### è sbagliato, sto passando dei rdata, invece che dataframe!!!
    
    Flags$ContinueUploading <- FALSE
    
    output$LoadingError3 <- renderText({
      validate(
        need(!is.null(input$RDataImportGrowthData$datapath), "Please select a Connector.RData!!!")
      )
      Flags$ContinueUploading = TRUE
      ""
    })
    
    observeEvent(Flags$ContinueUploading, {
      
      if(Flags$ContinueUploading)
      {
        if(is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList) )
        {
          load  = TRUE
        }else{ ### alert!!! if it is aleready present! 
          showModal(modalDialog(
            title = "Important message",
            "Do you want to update the ConnectorList already present?",
            easyClose = TRUE,
            footer= tagList(actionButton("confirmUpload2", "Update"),
                            modalButton("Cancel")
            )
          )) 
        }
        
        GrowthData <- readRDS(input$RDataImportGrowthData$datapath)
        
        observeEvent(input$confirmUpload2, {
          removeModal()
          if(!is.null(input$RDataImportFeatureData$datapath)){
            FeatureData <- readRDS(input$RDataImportFeatureData$datapath)
            txt <- withConsoleRedirect(
              CONNECTORList <- DataFrameImport(GrowthData,FeatureData)
            )
          }else{
            txt <- withConsoleRedirect(
              CONNECTORList <- DataFrameImport(GrowthData)
            )
          }
          rResult$CONNECTORList <- CONNECTORList
          output$SummaryLoading <- renderUI({HTML(paste(txt[-1], collapse = "<br/>"))} ) 
          Visual_1Step(CONNECTORList,input, output, session,rResult,Flags)
          
        })
        
      }
    })
    
    ## Now it is possible to save something:
    Flags$save1 = 0
  }) 
  observeEvent(input$saveRData1, {
    cat("save.image 1\n") 
    
    if( is.null(rResult$CONNECTORList) & is.null(rResult$CONNECTORList_trunc))
    {
      session$sendCustomMessage(type = 'testmessage',
                                message = 'There are no Connector objects to save!') 
    }else{
      CONNECTORList<-rResult$CONNECTORList
      CONNECTORList_trunc<-rResult$CONNECTORList_trunc
      Info <- reactiveValuesToList(rResult)
      
      save(CONNECTORList, CONNECTORList_trunc,Info, file = "Connector.RData")
      
      Flags$save1 = 0
    }
  })
  
  observe({
    if(Flags$save1 == 0){
      shinyjs::enable("download")
    }else{
      shinyjs::disable("download")
    }
  })
  
  output$download <- downloadHandler(
    filename = function(){
      paste("Connector.RData")
    },
    content = function(file) {
      #if(file_ext(file)!="RData") paste0(file_path_sans_ext(file),".RData")
      cat("file.copy\n")
      file.copy(from = "Connector.RData", to = file)
    },
    contentType = NULL
  )
  
## Report
  output$reportGen <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      report.path <- system.file("Shiny","report.Rmd", package = "connector")
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy(report.path, tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      infoReport <-  reactiveValuesToList(rResult)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = list(infoReport = infoReport),
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
## end
}


