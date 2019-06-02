#' Extrapolate objects from StabilityAnalysis
#'
#' @param stability.list CONNECTORList obtained from the StabilityAnalysis function. (see \code{\link{StabilityAnalysis}})
#' @param k  The number of clusters.
#' @param h  The  number representing the dimension of the cluster mean space (see \code{\link{PCA.Analysis}}).
#' 
#' @name ExtrapolationFuncs
NULL
#> NULL

#' @rdname ExtrapolationFuncs
#' @export
BoxPlot.Extrapolation <- function(stability.list, h){
  
  pl<-stability.list$BoxPlots[[paste("h= ",h)]]
    
  return(pl$Plot)
}

#' @rdname ExtrapolationFuncs
#' @export
ConsMatrix.Extrapolation <- function(stability.list, h, k){
  
  pl<-stability.list$ConsensusInfo[[paste("h= ",h)]][[paste("k= ",k)]]
  
  return(pl$ConsensusPlot)
}

