#' Maximum Discrimination Function
#'
#'@description
#'
#'
#' @param clusterdata Object belonging to the class funcyOutList if the model in study is the Functional Clustering Model (see \code{\link[funcy]{funcyOutList-class}}). Otherwise a list derived from fitting and clustering the data using Malthus, Gompertz or Logistic model storing the parameters and the cluster membership for each sample, the  parameters of the center and the mean curve values for each cluster (see \code{\link{FittingAndClustering}}).
#' @param k The number of clusters.
#' @param h The  vector/number representing the dimension of the cluster mean space. If NULL, ClusterChoice set the $h$ value equals to the number of PCA components needed to explain the 95\% of variability of the natural cubic spline coefficients, but the \emph{PCAperc} is needed (see \code{\link{PCA.Analysis}}).
#' @param absvalue 
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
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
#' 
MaximumDiscriminationFunction<-function(clusterdata,k=NULL,h=NULL,absvalue=TRUE)
{ 
  if( !isS4(clusterdata) )
  {
    warning("In input is needed a funcyOutList-class file. ")
  }else{
    fit<-clusterdata@models$fitfclust@fit
  }
  
  base <- fit$base
  sigma <- fit$par$sigma
  nt <- nrow(base)
  Gamma <- fit$par$Gamma
  Sigma <- base%*%Gamma%*%t(base)+sigma*diag(nt)
  Lambda <- fit$par$Lambda
  discrim <- solve(Sigma)%*%base%*%Lambda
  
  if (absvalue)
    discrim <- abs(discrim)
  n <- ncol(discrim)
  nrows <- ceiling(sqrt(n))
  
  par(mfrow=squareGrid(n))
  plots<-list()
  
  DiscrPlot<-function(i){
    DiscrFrame<-data.frame(Time=fit$grid,Discrim=discrim[,i])
    plots[[i]]<-ggplot(data=DiscrFrame,aes(x=Time,y=Discrim))+geom_line()+
                geom_hline(yintercept = 0, linetype="dashed")+
                xlab("Time")+ylab(paste("Discriminant Function ",i))
  }
  DiscriminantFunctions<-lapply(1:n,DiscrPlot)
  DiscriminantFunctionsALL<-plot_grid(plotlist = DiscriminantFunctions)
  return(list(DiscriminantFunctionsALL,DiscriminantFunctions))
}
