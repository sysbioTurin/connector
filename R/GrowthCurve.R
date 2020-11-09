#' Growth curves
#'
#' Generates the line plot of the growth data. The curves are colored with respect to the feature chosen by the user and/or reported in the AnnotationFile.
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param feature The column name reported in the AnnotationFile containing the feature to be investigated.
#' @param labels   The vector containing the text for the axis names and plot title.
#' @param save If TRUE then the plot is saved into a pdf file.
#' @param path The folder path  where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @return GrowthCurve returns a list containing the line plot as ggplot object, namely GrowthCurve_plot, and data, i.e. the CONNECTORList.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @examples
#' GrowDataFile<-"data/745dataset.xls"
#' AnnotationFile <-"data/745info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#' GrowthCurve(CONNECTORLis,"Progeny",labels=c("Time","Volume","Tumor Growth"))
#'
#' @import ggplot2
#' @export

GrowthCurve <- function(data,feature,labels=NULL,save=FALSE,path=NULL)
{
  if(is.null(labels))
  {
    axes.x<-""
    axes.y<-""
    title<-""

  }else{
    axes.x<-labels[1]
    axes.y<-labels[2]
    title<-labels[3]
  }


  ### dataframe for ggplot
  dataplot <- data$Dataset
  dataplot <- data.frame(merge(data$Dataset,data$LabCurv[,c("ID",feature)],by="ID"))
  dataplot[,feature]<-factor(as.matrix(dataplot[feature]))
  
  col <- as.character(unique(curves$Info))
  colFetaure <- rainbow(dim(unique(data$Lab[feature]))[1])
  
  ### Set growth curve plot with ggplot
  GrowthCurve <- ggplot(data=dataplot, aes(x=Time, y=Vol,group=ID,col=dataplot[,feature])) +
    geom_line() +
    geom_point() +
    labs(title=title,x=axes.x, y = axes.y,col=feature)+
    theme(plot.title = element_text(hjust = 0.5),title =element_text(size=10, face='bold'))+
    scale_colour_manual(values = colFetaure,limits = col, breaks = sort(col), name = feature)

  ### Enrich data with colour palette
  data$FeatureColour <- feature.palette

  GrowthCurve.ls <- list(GrowthCurve_plot=GrowthCurve,data=data,PlotData = dataplot)
  
  if(save==TRUE)
  {
    ggsave(filename="GrowthCurves.pdf",plot = GrowthCurve,width=29, height = 20, units = "cm",scale = 1,path=path )
  }

  return( GrowthCurve.ls )
}
