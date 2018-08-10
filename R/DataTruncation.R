#' DataTruncation
#'
#' @description
#'
#' Truncates the growth data at a specific time point chosen by the user.
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param feature The column name reported in the AnnotationFile containing the feature  to be investigated.
#' @param truncTime  An integer number corresponding to the time where  the curves will be truncated.
#' @param labels  Vector containing the text for the title of axis names and plot title.
#' @param save If TRUE then the growth curves plot truncated at the ``truncTime'' is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @return  DataTruncation returns the line plot of growth curves with a vertical line at the truncation time and the CONNECTORList updated with the following elements: (i)a data frame with three variables (ID curves, volume and time values truncated at the chosen time), (ii) a vector collecting the number of truncated observations collected per sample, a data frame matching curses  with the chosen feature, (iv) the vector storing all the truncated time points of the samples (i.e. truncated time grid).
#'
#' @examples
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))
#'
#' @import ggplot2
#' @export
DataTruncation <- function(data,feature,truncTime=NULL,labels=NULL,save=FALSE,path=NULL)
{

### Variables initialization
growth.curve.ls <- GrowthCurve(data,feature,labels = labels)
### Plot growth curves with truncation time
growth.curve.tr <- growth.curve.ls$GrowthCurve_plot
if(! is.null(truncTime)) growth.curve.tr <- growth.curve.tr + geom_vline(xintercept=truncTime, color="black", size=1)

if(save==TRUE)
{
 if(is.null(path))
 {
  path <- getwd()
 }
 ggsave(filename="DataTruncation.pdf",plot =growth.curve.tr,width=29, height = 20, units = "cm",scale = 1,path=path )
}

### Truncated dataset
dataset <- data$Dataset
sample.size <- max(unique(dataset[,1]))
lencurv.tr <- numeric(sample.size)

# Data truncation
if(!is.null(truncTime)) data.tr <- DataTrunc(data,truncTime=truncTime)
else data.tr <- data



plot(growth.curve.tr)

return(data.tr)
}

#######################à

DataTrunc <- function(data,truncTime=NULL)
{
#### Truncates the growth data at a specific time point chosen by the user.
  # Variables inizialization
  dataset <- data$Dataset
  sample.size <- max(unique(dataset[,1]))
  lencurv.tr <- numeric(sample.size)
  
  # Data truncation
  
  if(!is.null(truncTime))
  {
    max.time<-max(dataset[,3])
    if(max.time<truncTime)  warning("Truncation time greater than maximum time in the dataset.")
    dataset.tr <- dataset[dataset[,3]<=truncTime,]
    for (i in 1:sample.size)  lencurv.tr[i]<-length(dataset.tr[,1][dataset.tr[,1]==i])
    timegrid.tr <- data$TimeGrid[data$TimeGrid<=truncTime]
    data.tr=list(Dataset=dataset.tr,LenCurv=lencurv.tr,LabCurv=data$LabCurv,TimeGrid=timegrid.tr)
  }
  else data.tr <- data
  
  return(data.tr)
##### The CONNECTORList updated with the following arguments: a data frame with three variables (ID curves, volume and time values truncated at the chosen time), a vector collecting the number of truncated observations collected per sample, a data frame with curves labeled according to target file feature chosen and a vector for overall truncated time grid.
}
