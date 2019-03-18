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
#' @details   ADD THE COEFFS MEANING, to check the distance used in kmeans to calculate withtot 
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

    if(!is.null(clusterdata$fit))
    {
      k<-length(clusterdata$prediction$meancurves[1,])
      
      clusterdata$cluster$cluster.member -> classes
      clusterdata$prediction$meancurves -> meancurves
      clusterdata$fit$grid -> grid
      
    }else{
      clusterdata -> classification
      classification$cluster ->classes
      classification$meancurves->meancurves
      length(meancurves[1,])->k
      grid <- unique(data$Dataset$Time)
    }
    
    symbols<-cluster.symbol(k)
    

    classificate <- rep(classes,data$LenCurv)
    
    curves <- data.frame(Times=data$Dataset$Time,Vol=data$Dataset$Vol,ID=data$Dataset$ID,Cluster=classificate, Info=rep(t(data$LabCurv[feature]),data$LenCurv))
    

    
    # cut the meancurves at the growth curves' maximum time
    
    meancurves_truncated<-c()
    time3<-c()
    cluster<-c()
    
    for(clust in 1:k)
    {
      time2<-sort(unique(curves[curves$Cluster==clust,]$Times))
      m<-meancurves[,clust][grid%in%time2]
      time3<-c(time3,grid[grid%in%time2])
      meancurves_truncated<-c(meancurves_truncated,m)
      cluster<-c(cluster,rep(symbols[clust],length(grid[grid%in%time2])))
    }
    
    #cluster<-rep(symbols,each=length(grid))
    #plot_data<-data.frame(time=rep(grid,k),means=c(meancurves),clusters=cluster)
    
    plot_data<-data.frame(time=time3,means=meancurves_truncated,clusters=cluster)
    
    PlotMeanCurve<-ggplot()+
      geom_line(data=plot_data, aes(x=time,y=means,group=clusters,col= as.factor(clusters)) )+
      labs(title=title, x=axis.x, y = axis.y,colour="Cluster")+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    #+labs(subtitle = paste("Tight.E\ =\ ",as.integer(tightness$EucTight),"\ \ \ Tight.H\ =\ ",as.integer(tightness$HausTight),"\ \ \ coeff\ =\ ",signif(tightness$coeff , digits = 2)) )
    
      col<- as.character(unique(curves$Info))
      col1<- data$ColFeature
      plots<-list()
      plots_spline<-list()
      ymax<-max(plot_data$means,curves$Vol)
      xmax<-max(grid)
      
      for(i in 1:k)
      {
        plots[[paste(symbols[i],"Cluster")]]<-ggplot()+
          geom_line(data=plot_data[plot_data$clusters==symbols[i],], aes(x=time,y=means),size = 1.2 )+
          labs(title=paste(title,"",symbols[i],"Cluster"), x=axis.x, y = axis.y)+
          geom_line(data = curves[curves$Cluster==i,],aes(x=Times,y=Vol,group=ID,color=factor(Info)))+
          scale_colour_manual(values = col1,limits=col,breaks=col,name=feature)+
          theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())+
          ylim(0,ymax)+xlim(0,xmax)
      }
      
       plots[["ALL"]]<-plot_grid(plotlist = plots)
       
       
       if(save)
       {
         if(is.null(path)) path <- getwd()
         
         for(i in 1:k)
         {
           ggsave(filename = paste(symbols[i],"Cluster.pdf",sep=""),plot=plots[[paste(symbols[i],"Cluster")]],width=29, height = 20, units = "cm",scale = 1,path = path)
         }
           ggsave(filename = "ALLCluster.pdf",plot=plots$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
           ggsave(filename = "MeanCurve.pdf",plot=PlotMeanCurve,width=29, height = 20, units = "cm",scale = 1,path = path)
       }
         
######## Let's plot the spline fitting for each sample 
         
    if(!is.null(clusterdata$fit))
     {
      
      allsplineEstimation<-curve_prediction(cluster=classes ,object=clusterdata$prediction)

      out<-list(plotMeanCurve=PlotMeanCurve,plotsCluster=plots,spline.plots=allsplineEstimation)
    }
    else{
      out<-list(plotMeanCurve=PlotMeanCurve,plotsCluster=plots)
    }

  return(out)
}

curve_prediction<-function(cluster,object)
{
  index <- 1:length(table(object$data$curve)) #pdx
  
  timeIndx <- object$data$timeindex
  curveIndx <- object$data$curve
  
 plotCreation<-function(i,cluster,object)
 {
    cl<-cluster[i]
    grid <- object$grid
    upci <- object$upci[i,]
    uppi <- object$uppi[i,]
    lowci <- object$lowci[i,]
    lowpi <- object$lowpi[i,]
    gpred <- object$gpred[i,]
    meancurves <- (object$mean)[,cl]

        yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
    
    data.ggplot<-data.frame(grid=grid,upci=upci,uppi=uppi,lowci=lowci,lowpi=lowpi,gpred=gpred,meancurves=meancurves)
    data.real<-data.frame(time=grid[timeIndx[curveIndx==i]],vol=object$data$x[curveIndx==i] )
    
    gpl<-
      
      ggplot()+
      geom_ribbon(data=data.ggplot,aes(x=grid,ymin=lowci, ymax=upci), alpha=0.1)+
      geom_line(data=data.ggplot,aes(x=grid,y=gpred,linetype="Spline estimated",col="Spline estimated"))+
      geom_line(data=data.ggplot,aes(x=grid,y=meancurves,linetype="Cluster mean",col="Cluster mean"))+
      #geom_line(data=data.real,aes(x=time,y=gpredF,linetype="SplineFiltering",col="SplineFiltering"))+
      geom_line(data=data.real,aes(x=time,y=vol,col="Real points",linetype="Real points"))+
      geom_point(data=data.real,aes(x=time,y=vol),col="blue")+
      labs(title=paste("Sample ",i), x="Time", y="Growth value")+
      scale_colour_manual(values = c("black","red","blue") ,limits =c("Cluster mean","Spline estimated","Real points"),breaks= c("Cluster mean","Spline estimated","Real points") , name=" ")+
      guides(
             linetype = FALSE,
             colour = guide_legend(override.aes = list(linetype = c("solid","dashed","dashed"))))+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.key.width = unit(1, "cm"))
    
    
    return(gpl)
 }
 
 plot_list<-lapply(index,plotCreation,cluster,object)
 names(plot_list)<-paste(paste("Sample ",index))

 return(plot_list)
}

