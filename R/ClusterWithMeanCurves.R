#' Mean Cluster Curves
#'
#'@description
#' Visualizes the plots regarding the fitted and clustered data with respect to one model among FCM, Malthus, Gompertz and Logistic.
#'
#' @param clusterdata Object belonging to the class funcyOutList if the model in study is the Functional Clustering Model (see \code{\link[funcy]{funcyOutList-class}}). Otherwise a list derived from fitting and clustering the data using Malthus, Gompertz or Logistic model storing the parameters and the cluster membership for each sample, the  parameters of the center and the mean curve values for each cluster (see \code{\link{FittingAndClustering}}).
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param feature he column name reported in the AnnotationFile containing the feature  to be investigated.
#' @param title The  string containing  the plot title. 
#' @param labels  The vector containing the axis names. 
#' @param save If TRUE then the plots of the fitted and clustered growth curves are saved into two  pdf files.
#' @param path  The folder path where the plots will be saved. If it is missing, the plots are saved in the current working  directory.
#' 
#' @return ClusterWithMeanCurve returns a list with two objects, (i) a list storing the growth curves plots partitioned according to cluster membership and (ii) the cluster mean curves plot.
#'
#' @seealso \code{\link{ClusterChoice}},  \code{\link{FittingAndClustering}}.
#'
#' @examples
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'### Merge curves and target file
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'### Truncation
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",60,labels = c("time","volume","Tumor Growth"))
#'
#'
#' ###  FCM
#'
#'CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=2)
#'
#'CONNECTORList.FCM.k4.h2<- CONNECTORList.FCM$FCM_all$`k= 4`$`h= 2`
#'
#'FCMplots <- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.k4.h2,data= CONNECTORList,feature = "Progeny",labels = c("Time","Volume"),title= " FCM model ")
#'
#' ###  Malthus
#'
#' lower<-c(10^(-5),0)
#' upper<-c(10^2,10^3)
#' init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=1)
#'
#'Malthus1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",feature="Progeny",fitting.method="optimr",lower=lower,upper=upper,init=init)
#'MalthusPlots1<-ClusterWithMeanCurve(clusterdata=Malthus1,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "Optimr Malthus model")
#'
#'MalthusPlots1$plotsCluster$ALL
#'
#' @import ggplot2 
#' @importFrom cowplot plot_grid
#' @export
ClusterWithMeanCurve<-function(clusterdata, data, feature ,title="", labels=c("","") ,save=FALSE,path=NULL )
{
axis.x<-labels[1]
axis.y<-labels[2]
  
    if(isS4(clusterdata))
    {
      k<-clusterdata@k
      
      Cluster(clusterdata)->classes
      out.fit<-clusterdata@models$fitfclust@fit
      # Check if it is regular
      if(clusterdata@reg==1)
      {
        fitfclust.curvepred(out.fit)$meancurves->meancurves
      } else{
        fitfclust.curvepredIrreg(out.fit)$meancurves->meancurves
      }
    }else{
      
      clusterdata -> classification
      classification$cluster ->classes
      classification$meancurves->meancurves
      length(meancurves[1,])->k
    }
    
    symbols<-cluster.symbol(k)
    
    ##################
    
    classificate <- rep(classes,data$LenCurv)
    curves <- data.frame(Times=data$Dataset$Time,Vol=data$Dataset$Vol,ID=data$Dataset$ID,Cluster=classificate, Info=rep(t(data$LabCurv[feature]),data$LenCurv))
  
    # cut the meancurves at the growth curves' maximum time
    time1<-sort(unique(data$Dataset$Time))
    meancurves_truncated<-c()
    time3<-c()
    cluster<-c()
    for(clust in 1:k)
    {
      time2<-sort(unique(curves[curves$Cluster==clust,]$Times))
      m<-meancurves[,clust][time1<=max(time2)]
      time3<-c(time3,time1[time1<=max(time2)])
      meancurves_truncated<-c(meancurves_truncated,m)
      cluster<-c(cluster,rep(clust,length(time1[time1<=max(time2)])))
    }
  
  
  
    plot_data<-data.frame(time=time3,means=meancurves_truncated,clusters=cluster)
    PlotMeanCurve<-ggplot()+
                      geom_line(data=plot_data, aes(x=time,y=means,group=clusters) )+
                      labs(title=title, x=axis.x, y = axis.y)+
                      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
  
      col<-as.character(unique(curves$Info))
      col1<-rainbow(length(col))
      plots<-list()
      ymax<-max(curves$Vol)
      for(i in 1:k)
      {
        plots[[paste(symbols[i],"Cluster")]]<-ggplot()+
          geom_line(data=plot_data[plot_data$clusters==i,], aes(x=time,y=means),size = 1.2 )+
          labs(title=paste(title,"",symbols[i],"Cluster"), x=axis.x, y = axis.y)+
          geom_line(data = curves[curves$Cluster==i,],aes(x=Times,y=Vol,group=ID,color=factor(Info)))+
          scale_colour_manual(values = col1,limits=col,breaks=col,name=feature)+
          theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())+
          ylim(0,ymax)
      }
       plots[["ALL"]]<-plot_grid(plotlist = plots)
       
       out<-list(plotMeanCurve=PlotMeanCurve,plotsCluster=plots)

  if(save)
  {
    if(is.null(path)) path <- getwd()
    for(i in 1:k)
    {
    ggsave(filename = paste(symbols[i],"Cluster.pdf",sep=""),plot=out$plotsCluster[[paste(symbols[i],"Cluster")]],width=29, height = 20, units = "cm",scale = 1,path = path)
    }
    ggsave(filename = "ALLCluster.pdf",plot=out$plotsCluster$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
    ggsave(filename = "MeanCurve.pdf",plot=out$plotMeanCurve,width=29, height = 20, units = "cm",scale = 1,path = path)
}

 return(out)
}


