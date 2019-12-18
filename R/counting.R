#' Counting samples
#'
#' @description
#' Counts for each combination of cluster and feature the number of samples and their name are returned.
#'
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{MostProbableClustering.Extrapolation}}). 
#' @param data CONNECTORList.  (see \code{\link{DataImport}} or \code{\link{DataTruncation}})
#' @param feature The column name reported in the AnnotationFile containing the feature  to be investigated.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return CountingSample returns a list containing  1) a matrix composed by three columns: (i) cluster, (ii) the feature name, and (iii) the number of samples, and 2) a dataframe storing the sample names reported in the AnnotationFile and the corresponding  cluster membership.
#' 
#' 
#' @examples
#'
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
#' NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.k4.h2,CONNECTORList,feature = "Progeny")
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
#' NumberSamples<-CountingSamples(clusterdata=Malthus1,CONNECTORList,feature = "Progeny")
#'
#' @import plyr
#' @export
#' 
CountingSamples<-function(clusterdata,data,feature="ID")
{
  if( is.null(clusterdata$FCM) &  is.null(clusterdata$cluster) )
  {
    warning("In input is needed the FCM or StabilityAnalysis file. ",immediate. = T)
  }else{
    if(!is.null(clusterdata$FCM))    clusterdata <- clusterdata$FCM
  }
  
  ClustCurve <- clusterdata$cluster$ClustCurve
  
  ClustCurve <- data.frame(merge(ClustCurve,data$LabCurv,by="ID"))
  
  ClustCurve$Cluster <- clusterdata$cluster$cluster.names[ClustCurve$Cluster]
  
  countRes<-count(ClustCurve,
                  c("ID", "Cluster",feature))
  countRes<-countRes[,-length(count(ClustCurve, c("ID", "Cluster",feature)))]
  
  
  Counting<-count(countRes,c( "Cluster",feature))

  return(list(Counting=Counting,ClusterNames=data.frame(Cluster=ClustCurve$Cluster[cumsum(CONNECTORList$LenCurv)],SampleName=unique(ClustCurve$SampleName))))
}

