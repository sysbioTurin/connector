#' Save the plots of the spline exploited to fit and cluster the samples.
#' 
#' @description 
#' Saves the Cubic Spline plots from the ClusterWithMeanCurve output list (see \code{\link{ClusterWithMeanCurve}}). Each plot shows (i) in blue the sample curve, (ii) in red the cubic spline estimated from the FCM, (iii) in black the correspondive cluster mean curve, and finally (iv) the grey area represents the confidence interval.
#'
#' @param FCM.plots CONNECTORList obtained from the ClusterWithMeanCurve function. (see \code{\link{ClusterWithMeanCurve}})
#' @param All If TRUE then all the plots, for each sample curve, are saved.
#' @param SampleNumber The number/vector of the sample ID(s) referring to the curve(s) that will be saved.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#'  
#'  
#' @import ggplot2
#' 
Spline.plots <- function(FCM.plots, All=TRUE, SampleNumber = NULL, path = NULL){
  
  n.samples<-length(FCM.plots$spline.plots)
  
  if(ALL){ 
    samplesIndexes<- 1:n.samples
  }else if(!is.null(SampleNumber)) {
    # select only the sample number between 1 and n.samples
    samplesIndexes<-SampleNumber[which( SampleNumber %in% 1:n.samples == TRUE)]
    if(length(samplesIndexes)==0) {
      warning(paste0("The number of the sample must be between 0 and ",n.samples,", that is the number of samples. It will be saved the first plot."),immediate. = TRUE)
      samplesIndexes<-1
    }
  }else{
    warning(paste0("Select a number of the sample between 0 and ",n.samples,", or put the parameter ALL = True. Since neither the parameter All or SampleNumber were selected, it will be saved only the first curve."),immediate. = TRUE)
    samplesIndexes<-1
  }
  
  if(is.null(path)) path <- getwd()
      
      for(i in samplesIndexes)
      {
        ggsave(filename = paste("Spline",i,"sample.pdf",sep="_"),plot=FCM.plots$spline.plots[[paste("Sample",i)]],width=29, height = 20, units = "cm",scale = 1,path = path)
      }
      
}