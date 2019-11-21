#' Discriminant Plot
#'
#'@description
#'
#'
#' @param clusterdata Object belonging to the class funcyOutList if the model in study is the Functional Clustering Model (see \code{\link[funcy]{funcyOutList-class}}). 
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param h The  number between 1 or 2 representing the dimension of the cluster mean space(see \code{\link{PCA.Analysis}}).
#' @param save If TRUE then the growth curves plot truncated at the ``truncTime'' is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#'  @return
#' 
#' @seealso \code{\link[funcy]{funcit}}.
#' 
#' @references
#' Gareth M. James and Catherine A. Sugar, (2003). Clustering for Sparsely Sampled Functional Data. Journal of the American Statistical Association.
#' 
#' @examples
#'
#'
#' @import  ggplot2 
#' @export
#' 

DiscriminantPlot<-function(clusterdata,h,data,feature,save=FALSE,path=NULL)
{ 
  G<-length(clusterdata$FCM$prediction$meancurves[1,])
  discrplot<-list()
  
  outpred<-clusterdata$FCM$cluster$cl.info
  
  projectedcurve=outpred$alpha.hat
  projectedclustcenters=clusterdata$FCM$fit$parameters$alpha
  outpred$Calpha -> stdevalpha
  outpred$class.pred -> classes
  outpred$distance -> distanze
  Feature <- data$LabCurv[paste(feature)]
  
  col<- as.character(unlist(unique(Feature)) )
  col1<- data$ColFeature
  
  symbols<-clusterdata$FCM$cluster$cluster.names
  
  if(h==1)
  {     
    DataFrame<-data.frame(PrjCurv=projectedcurve,stdaplh=stdevalpha,Feature=unlist(Feature),Cluster=symbols[classes] )
    
    cl.names<-1:G-1
    names(cl.names)<-symbols
    
    discrplot[["ColCluster"]]<-ggplot(data = DataFrame)+
      geom_point(aes(x=PrjCurv,y=stdaplh,colour=Cluster,pch=Cluster),size=4)+
      geom_vline(xintercept = projectedclustcenters)+
      scale_shape_manual("Cluster",values=cl.names)+ 
      labs(title="Discriminant plot")+
      xlab("Alpha")+ylab('Standard Deviation')+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    

    discrplot[["ColFeature"]]<-ggplot(data = DataFrame)+
      geom_point(aes(x=PrjCurv,y=stdaplh,colour=as.factor(Feature),shape=Cluster),size=4)+
      geom_vline(xintercept = projectedclustcenters)+
      xlab("Alpha")+ylab('Standard Deviation')+
      labs(title="Discriminant plot")+
      scale_colour_manual(values = col1,limits=col,breaks=col,name=feature )+
      scale_shape_manual("Cluster",values=cl.names) +
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
  }
  else{
    DataFrameSamples<-data.frame(PrjCurv1=projectedcurve[,1],PrjCurv2=projectedcurve[,2],Feature=unlist(Feature),Cluster=symbols[classes] )
    DataFrameCluster<-data.frame(projectedclustcenters1=projectedclustcenters[,1],projectedclustcenters2=projectedclustcenters[,2],Cluster=symbols,Center=paste("c.",symbols,sep=""))
    
    cl.names<-1:G-1
    names(cl.names)<-symbols
    
    discrplot[["ColCluster"]]<-ggplot()+
      geom_point(data = DataFrameSamples,aes(x=PrjCurv1,y=PrjCurv2,colour=Cluster,pch=Cluster),size=4)+
      geom_text(data = DataFrameCluster,aes(x=projectedclustcenters1,y=projectedclustcenters2,label=paste("c",symbols,sep=""),colour=Cluster),size=5,show.legend = F)+
      scale_shape_manual("Cluster",values=cl.names,labels=paste(names(cl.names)," (c",symbols,")",sep=""),breaks=names(cl.names)) +
      scale_colour_discrete("Cluster",limits=names(cl.names),breaks=names(cl.names),labels=paste(names(cl.names)," (c",symbols,")",sep=""))+
      xlab("Alpha 1")+ylab('Alpha 2')+
      labs(title="Discriminant plot")+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    
    discrplot[["ColFeature"]]<-ggplot()+geom_point(data = DataFrameSamples,aes(x=PrjCurv1,y=PrjCurv2,colour=as.factor(Feature),pch=Cluster),size=4)+
      geom_text(data = DataFrameCluster,aes(x=projectedclustcenters1,y=projectedclustcenters2,label=Center),size=5,show.legend = F)+
      xlab("Alpha 1")+ylab('Alpha 2')+
      labs(title="Discriminant plot")+
      scale_colour_manual(values = col1,limits=col,breaks=col,name=feature )+
      scale_shape_manual("Cluster",values=cl.names,labels=paste(names(cl.names),", c.",symbols,sep=""),breaks=names(cl.names)) +
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
  }
  
  if(save==TRUE)
  {
    if(is.null(path))
    {
      path <- getwd()
    }
    ggsave(filename="DiscrPlot_ColFeature.pdf",plot =discrplot[["ColFeature"]],width=29, height = 20, units = "cm",scale = 1,path=path )
    ggsave(filename="DiscrPlot_ColClust.pdf",plot =discrplot[["ColCluster"]],width=29, height = 20, units = "cm",scale = 1,path=path )
  }
  
  return(discrplot)
}
