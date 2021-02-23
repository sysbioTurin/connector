#' Data Frame Import
#'
#'@description
#' Reads the data.frame and creates a list storing all the information.
#'
#'
#' @param GrowDataFrame Dataframe of three columns storing the time series. The first columns, called "\emph{ID}", stores the identification number of each sample. The second column, labeled "\emph{Vol}", contains  the data volume over the time and the respective time point is reported in the third column, labeled  "\emph{Time}"
#' @param AnnotationFrame If NULL, then in the ConnectorList only the feature ID will be reported. Otherwise a dataframe must be passed, with a number of rows equal to the number of samples in the GrowDataFrame. The first column must be named "ID" and it stores the identification numbers exploited for the matching with the samples saved in the GrowDataFrame. Then it is possible to add a column per feature to consider and to associate to the respective sample (depending on the identification number in the first column)
#'
#' @return   DataFrameImport returns a list, called ConnectorList, with four arguments: (i) Dataset: data frame with three columns (i.e. ID, data and time values) encoding the growth data for each sample, (ii) LenCurv: the vector reporting the number of observations collected per sample  (iii) LabCurv: the data frame matching the samples with their annotations, (iv) TimeGrid: the vector storing all the sample time points (i.e. time grid). Furthermore, it prints a brief summary of the input data, i.e. the total number of curves (samples), the minimum and the maximum curve length.
#' 
#' @details One data frame is requested to run the data analysis exploiting the CONNECTOR package:
#' \itemize{
#' \item the \emph{GrowDataFrame}, reporting the growth evolution data. 
#' }
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#' 
#' @examples
#'
#' @import readxl
#' @export
DataFrameImport <- function(GrowDataFrame,AnnotationFrame=NULL) {
  alldata<-list()
  
  ###Read growth curves dataframe
  colnames(GrowDataFrame) <-c("ID","Vol","Time")
  
  ID.sample <- unique(GrowDataFrame$ID)
  samplesize <- length(ID.sample)
  alldata$LenCurv <- table(GrowDataFrame$ID)
  alldata$LenCurv <- alldata$LenCurv[paste(unique(GrowDataFrame$ID))]
  ID <- 1: samplesize
  if(!all(ID==ID.sample))
  { # the id passed in input is correct, from 1 to number of samples
    GrowDataFrame$ID <- rep(ID, alldata$LenCurv )
    switchIndex <- 1:samplesize
    names(switchIndex) <- ID.sample

    changeAnnotation<-T
  }else{
    changeAnnotation<-F
  }
  
  if(is.null(AnnotationFrame))
  {
    if(changeAnnotation){
      annotations <- data.frame(ID = ID,ID.sample = ID.sample)
    }else{
      annotations <- data.frame(ID = ID)
    }
  }else{
    if(colnames(AnnotationFrame)[1]!= "ID") 
      {
        warning("The first column name of the AnnotationFrame must be 'ID' \n ")
        return()
      }else{
        if(unique(ID.sample %in% AnnotationFrame[,1]) %in% FALSE)
          {
            warning("The AnnotationFrame 'ID's do not correspond to the one stored in the GrowDataFrame!\n ")
            return()
          }
        }
    
    if(changeAnnotation){
      colnames(AnnotationFrame)[1] <- "ID.Sample"
      annotations <- data.frame(ID = switchIndex[paste(AnnotationFrame$ID.Sample)], AnnotationFrame )
    }else{
      annotations <- AnnotationFrame
    }
 
  }
  
  if(! "SampleName" %in% colnames(annotations) )
  {
    annotations$SampleName <- paste("Sample", ID)
  }
  ###Check the column names
  
  if(samplesize!=length(annotations[,1]) )
  {
    warning("Number of samples in the GrowDataFrame is different from the number of samples stored in the target file.")
    return()
  }
  
  ### Inizialize :
  alldata$Dataset <- GrowDataFrame
  alldata$TimeGrid <- sort(unique(GrowDataFrame$Time))
  alldata$LabCurv <- annotations[order(annotations$ID),]
    
  cat("############################### \n######## Summary ##############\n")
  cat("\n Number of curves:",samplesize,";\n Min curve length: ",min(alldata$LenCurv),"; Max curve length: ",max(alldata$LenCurv),".\n")
  cat("############################### \n")
  
  return(alldata)
}
