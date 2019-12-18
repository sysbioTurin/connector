#' Maximum Discrimination Function
#'
#'@description
#' Visualizes the h curve(s) representing the optimal weights to apply to each dimension for determining the cluster membership.
#'
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{MostProbableClustering.Extrapolation}}).
#' @param absvalue If TRUE, the absolute values of the weights are plotted.
#'  @return
#' MaximumDiscriminationFunction generates h plots as ggplot objects, showing the      the weights, applied to each dimension for determining the cluster membership, plotted as a single curve. 
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @references
#' Gareth M. James and Catherine A. Sugar, (2003). Clustering for Sparsely Sampled Functional Data. Journal of the American Statistical Association.
#' 
#' @examples
#'
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
#' 
MaximumDiscriminationFunction<-function(clusterdata,absvalue=TRUE)
{ 
  if( is.null(clusterdata$FCM) &  is.null(clusterdata$fit) )
  {
    warning("In input is needed the FCM or StabilityAnalysis file. ",immediate. = T)
  }else{
    if(!is.null(clusterdata$FCM))    fit<-clusterdata$FCM$fit
    if(!is.null(clusterdata$fit))    fit<-clusterdata$fit
  }
    
  
  base <- fit$FullS
  sigma <- fit$par$sigma
  nt <- nrow(base)
  Gamma <- fit$par$Gamma
  Sigma <- base%*%Gamma%*%t(base)+sigma*diag(nt)
  Lambda <- fit$par$Lambda
  discrim <- solve(Sigma)%*%base%*%Lambda
  
  if (absvalue)
    discrim <- abs(discrim)
  n <- ncol(discrim)
  
  plots<-list()
  
  DiscrPlot<-function(i){
    DiscrFrame<-data.frame(Time=fit$grid,Discrim=discrim[,i])
    plots[[i]]<-ggplot(data=DiscrFrame,aes(x=Time,y=Discrim))+geom_line()+
                geom_hline(yintercept = 0, linetype="dashed")+
                xlab("Time")+ylab(paste("Discriminant Function ",i))
  }
  DiscriminantFunctions<-lapply(1:n,DiscrPlot)
  
  if(n>1)  
  {
    DiscriminantFunctionsALL<-plot_grid(plotlist = DiscriminantFunctions)
    return(list(DiscriminantFunctionsALL,DiscriminantFunctions))
  }
  else return(DiscriminantFunctions)
}
