#' Time Grid Density
#'
#'@description
#'Plots the grid of the time points at which the growth data are collected. Each point of the grid is colored with respect to its frequency in the input data.
#'
#' @param  data CONNECTORList. (see \code{\link{DataImport}})
#' @param  save If TRUE then the density time grid plot is saved into a pdf file.
#' @param  path	 The folder path where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @return Time Grid Density returns the time grid density plot as a ggplot object.
#' 
#' @examples
#' GrowDataFile<-"data/745dataset.xls"
#' AnnotationFile <-"data/745info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#' TimeGridDensity(CONNECTORList)
#'
#' @import ggplot2 plyr
#' @export
TimeGridDensity <- function(data,save=FALSE,path=NULL)
{

 ### Variables initialization
 TimeMeasure <-data$Dataset[,c(1,3)]
 SampleSize <- length(unique(data$Dataset[,1]))
 LenCurve <- data$LenCurv
 PointsCoord<-matrix(0,nrow =LenCurve%*%LenCurve,ncol = 2)
 k<-0

 # Generate time grid
  for(i in 1:SampleSize)
  {
    l<-length(TimeMeasure[TimeMeasure[,1]==i,1])
    PointsCoord[(k+1):(k+l^2),1]<-rep(TimeMeasure[TimeMeasure[,1]==i,2],l)
    PointsCoord[(k+1):(k+l^2),2]<-rep(TimeMeasure[TimeMeasure[,1]==i,2],each=l)
    k<-k+l^2
  }

  ### Generate data frame for ggplot
  Time1 <- PointsCoord[,1]
  Time2 <- PointsCoord[,2]
  
  df <- data.frame(Time1=Time1,Time2=Time2)
  df<-ddply(df,.(Time1,Time2),nrow)
  

  ### Plot density grid
  TimeGrid_plot <- ggplot(df) +
                   geom_point(aes(Time1, Time2,col=V1), size = 4) +
	               coord_fixed(ratio = 1) +
                scale_color_gradientn(colours = c("#baffc9","#FF0000"),name="Number of \nobservations")+
                 theme_bw() +
	               labs(title="Time grid",x="Time", y = "Time")+
                   theme(plot.title = element_text(hjust = 0.5),title =element_text(size=12, face='bold'))

  if(save==TRUE)
  {
    if(is.null(path))
	{
	path <- getwd()
    ggsave(filename="TimeGrid.pdf",plot =TimeGrid_plot,width=29, height = 20, units = "cm",scale = 1,path=path )
    }
  }
  return(TimeGrid_plot)
}
