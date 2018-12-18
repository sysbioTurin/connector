#' Separation and Tightness Visualization
#'
#' Visualizes the plots reporting the Separation and Tightness clustering measures.
#'
#' @param clusterdata Object belonging to the class funcyOutList if the model in study is the Functional Clustering Model (see \code{\link[funcy]{funcyOutList-class}}). Otherwise a list derived from fitting and clustering the data using Malthus, Gompertz or Logistic model storing the parameters and the cluster membership for each sample, the  parameters of the center and the mean curve values for each cluster (see \code{\link{FittingAndClustering}}).
#' @param Title Title of the Separation and Tightness plot.
#' @param save If TRUE then the Separation and Tightness plot is saved into a pdf file.
#' @param path The folder path where the plots will be saved. If it is missing, the plots are saved in the current working  directory.
#' 
#' @return  Returns the Separation and Tightness plot.
#' 
#' @details To visualize the Separation (how much the cluster mean curves are different from each other)  and the Tightness (how much the curves belonging to the same cluster are different from the corresponding mean curve) of the clusters three circles are drawn according to the mean and variance of the Euclidian distance of fitted curves within a cluster with respect to the cluster mean curve. In each circle the number represent the cluster mean curves while the symbols represent the objects in the cluster, one symbol for each cluster.
#' 
#' 
#' @examples
#' ### Data files
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#' 
#' ### Merge curves and target file
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#' 
#'### Truncation
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",60,labels = c("time","volume","Tumor Growth"))
#'
#' ############################## FCM ##############################
#'
#' CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=2)
#' CONNECTORList.FCM.k4.h2<- CONNECTORList.FCM$FCM_all$`k= 4`$`h= 2`
#'
#' PlotSeparationTightness(CONNECTORList.FCM.k4.h2,Title = "FCM Cluster betweenness and withinness")
#'
#' ############################## MALTHUS ##############################
#' lower<-c(10^(-5),0)
#' upper<-c(10^2,10^3)
#' init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=1)
#' 
#' 
#' Malthus1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",feature="Progeny",fitting.method="optimr",lower=lower,upper=upper,init=init)
#'
#'
#' PlotSeparationTightness(Malthus1,Title = "Malthus Cluster betweenness and withinness")
#' 
#' @import ggplot2 ggforce Matrix
#' @export
PlotSeparationTightness <- function(clusterdata,Title=NULL,save=FALSE,path=NULL)
{
  
  
  if(isS4(clusterdata))
  {
    Cluster(clusterdata)->cluster
    out.fit<-clusterdata@models$fitfclust@fit
    # Check if it is regular
    if(clusterdata@reg==1)
    {
      fitfclust.curvepred(out.fit)$meancurves-> MeanCurves
    } else{
      fitfclust.curvepredIrreg(out.fit)$meancurves-> MeanCurves
    }
    ClustCurve <- out.fit$ClustCurves
    
  }else{
    
    MeanCurves <- clusterdata$meancurves
    ClustCurve <- clusterdata$Summary
    clusterdata$cluster -> cluster
  }
  
  # calculating the area to sort the cluster
  
  distFromZero<-apply(MeanCurves,2,sum)/length(MeanCurves[,1])
  
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol<-cluster.symbol(K)
  cluster.palette <- rainbow(K)
  
  names(distFromZero)<-ClustSymbol
  
  distFromZero.sorted<-sort(distFromZero)
  ClustSymbol.sorted <- names(distFromZero.sorted)
  i<-match(ClustSymbol.sorted, ClustSymbol)
  MeanCurves.sorted<-MeanCurves[,i]
  
  nc <- dim(MeanCurves.sorted)[2]
  nt <- dim(MeanCurves.sorted)[1]
  res <- t(sapply(2:nc, function(x) sum(sqrt((MeanCurves.sorted[,x]-MeanCurves.sorted[,1])^2))) /nt)
  
  
  
  shift<-sapply(2:K , function(x){
    sum( sqrt((MeanCurves.sorted[,1]-MeanCurves.sorted[,x])^2 ) )
  } )/length(MeanCurves.sorted[,1]) 
  
  shift<- c(0,shift) # taking the first cluster as zero!
  
  cluster.magnitudo <- apply(table(unique(ClustCurve[,c(4,1)])),1,sum)[i]
  circles <- data.frame(
    x0 = numeric(3*K)     ,
    y0 = numeric(3*K)     ,
    r = numeric(3*K)      ,
    distance = numeric(3*K),
    Cluster = numeric(3*K)
  )
  colnames(circles) <- c("x0","y0","r","distance","Cluster")
  
  n <- length(unique(ClustCurve[,1]))
  WithDist <- data.frame(
    x1 = numeric(n),
    y1 = numeric(n),
    Cluster = numeric(n)
  )
  colnames(WithDist) <- c("x1","y1","Cluster")
  index <- matrix(seq(1,3*K),nrow=K,byrow=TRUE)
  
  dataplotALL <- DataFrameWithinness.2(clusterdata,MeanCurves,cluster,shift,ClustSymbol=ClustSymbol.sorted,K)
  
  counter <- 1
  for(k in 1:K)
  {
    ### Data frames
    dataplot <- dataplotALL[[i[k]]]
    circles[index[k,],] <- dataplot$circles
    WithDist[counter:(cumsum(cluster.magnitudo)[k]),] <- dataplot$WithDist
    counter <- cumsum(cluster.magnitudo)[k] + 1
  }
  
  circles$distance <- factor(circles$distance)
  circles$Cluster <- factor(circles$Cluster)
  
  WithDist$Cluster <- factor(WithDist$Cluster)
  
  plots <- ggplot() + geom_circle(aes(x0=x0, y0=y0, r=r,linetype=distance,color=Cluster), data=circles,size=.6)+ylim(-max(shift)/2,max(shift)/2)
  plots <- plots + scale_shape_manual(values=c(0:(K-1)))+
    geom_point(data=WithDist,aes(x=x1,y=y1,shape=Cluster),size=4)+
    scale_linetype_manual("",values = c( "2"= "dashed","1"="solid"),
                          labels=c("Mean ","Mean+- sd"))
  
  if(is.null(Title)) Title<-"Cluster Separation and Tightness"
  
  plots <- plots  + geom_text(aes(x=x0, y=y0,label=Cluster), data=circles) +
    #scale_colour_manual(values = cluster.palette,name="Cluster Colors") +
    labs(title=Title,x="distance",y="distance")+
    theme(plot.title = element_text(hjust = 0.5),                                                                      axis.line = element_line(colour = "black"))
  
  
  if(save==TRUE)
  {
    if(is.null(path)) path <- getwd()
    ggsave(filename="Separation&Tightness.pdf",plot =plots,width=29, height = 20, units = "cm",scale = 1,path=path)
  }
  
  return(plots)
}