#' Extrapolate objects from StabilityAnalysis
#' 
#' @description Extrapolate objects from StabilityAnalysis output list (see \code{\link{StabilityAnalysis}}).
#'
#' @param stability.list The list obtained from the StabilityAnalysis function. (see \code{\link{StabilityAnalysis}})
#' @param G  The number of clusters.
#' @param h  The number representing the dimension of the cluster mean space (see \code{\link{PCA.Analysis}}).
#' 
#' @details 
#' \itemize{
#' \item{BoxPlot.Extrapolation}{ extrapolates from StabilityAnalysis output list the box plot fixing the h value.}
#' \item{ConsMatrix.Extrapolation}{extrapolates from StabilityAnalysis output list the Consensus Matrix for G and h fixed.}
#' \item{MostProbableClustering.Extrapolation}{extrapolates from StabilityAnalysis output list the most probable cluster among the several runs obtained  for G and h fixed.}
#' }
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @name ExtrapolationFuncs
NULL
#> NULL

#' @rdname ExtrapolationFuncs
#' @export
BoxPlot.Extrapolation <- function(stability.list, h){
  
  pl<-stability.list$BoxPlots[[paste("h=",h)]]
    
  return(pl$Plot)
}

#' @rdname ExtrapolationFuncs
#' @export
ConsMatrix.Extrapolation <- function(stability.list, h, G){
  
  pl<-stability.list$ConsensusInfo[[paste("h=",h)]][[paste("G=",G)]]
  
  return(pl$ConsensusPlot)
}

#' @rdname ExtrapolationFuncs
#' @export
MostProbableClustering.Extrapolation <- function(stability.list, h, G){
  
  pl<-stability.list$ConsensusInfo[[paste("h=",h)]][[paste("G=",G)]]
  
  return(pl$MostProbableClustering)
}

