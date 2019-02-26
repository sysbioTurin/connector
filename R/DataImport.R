#' Data Import
#'
#'@description
#' Reads the files and creates a list storing all the information.
#'
#'
#' @param GrowDataFile The name of the excel file storing the  growth evolution data. See \sQuote{Details}.
#'
#' @param AnnotationFile The name of a csv file storing  the annotation data. See \sQuote{Details}.
#'
#' @return   DataImport returns a list called CONNECTORList, with four arguments: (i) data frame with three columns (i.e. ID, data and time values) encoding the growth data for each sample, (ii) the vector reporting the number of observations collected per sample  (iii) the data frame matching the samples with their annotations, (iv) the vector storing all the sample time points (i.e. time grid).
#' 
#' @details Two files are requested to run the data analysis exploiting the CONNECTOR package:
#' \itemize{
#' \item the excel file, namely \emph{GrowDataFile}, reporting the growth evolution data,
#' \item the csv file, namely \emph{AnnotationFile}, containing the annotation information associated with the samples.
#' }
#' Hence, the growth data associated with an experiment must be stored into GrowDataFile as a table with two columns for each sample.
#' The first column, labeled  \emph{Time}, contains the time points of a sample. The second column, labeled by the sample name, contains  the data volume over the time.
#' 
#' Instead, the second file (i.e. AnnotationFile)  stores  the annotated information associated with the samples as a table in csv format so that  number of rows is equals to the total number of samples.
#'
#' @examples
#'
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#' @import readxl
#' @export
DataImport <- function(GrowDataFile,AnnotationFile) {
 ###Read Data File
  dataset <- read_excel(GrowDataFile,col_names=T)
 ### Read Target File
  labcurv  <- read.csv(file=AnnotationFile,header=TRUE)

 ###Check the column names
  c_names<-colnames(dataset[2*(1:(length(dataset[1,])/2))])
  if(length(c_names)!=(length(labcurv$ID)))
  {
    warning("Number of columns in the excel file is different from the number of samples stored in the target file.")

  }else{

    if(all(c_names==labcurv$SampleName)==FALSE)
    {
      warning("SampleNames in the target file do not correspond to the names in the excel file.")
    }
  }



  ### Inizialize :
  ### vector for curves lenghts
  nvar       <- dim(dataset)[2]
  nobs       <- dim(dataset)[1]
  tot        <- nobs*nvar
  samplesize <- dim(dataset)[2]/2
  lencurv    <- numeric(samplesize)
  ### vectors for times, volume and ID curves values
  tot        <- length(dataset)
  TimeIndex  <- seq(1,nvar,2)
  VolIndex   <- seq(2,nvar,2)
  TimeValue  <- as.matrix(dataset[TimeIndex])
  VolValue   <- as.matrix(dataset[VolIndex])

  times      <- numeric(tot)
  vol        <- numeric(tot)
  ID         <- numeric(tot)

  ### Organize times, volume and ID curves values, removing NA
  for (cappa in 1:samplesize)
  {
    tempv <- as.double(VolValue[,cappa])
    tempv <- tempv[!is.na(tempv)]
    lencurv[cappa] <- length(tempv)
    tempt <-as.double(TimeValue[,cappa])
    tempt <-tempt[!is.na(tempt)]

    if (cappa == 1)
    {
      vol[1:lencurv[cappa]]    <- tempv
      times[1:lencurv[cappa]]  <- tempt

    }

    else
    {
      lcum <- cumsum(lencurv)
      vol[(lcum[cappa-1]+1):lcum[cappa]] <- tempv
      times[(lcum[cappa-1]+1):lcum[cappa]] <- tempt
    }
  }

  ndata    <- sum(lencurv)
  vol      <- vol[1:ndata]
  times    <- times[1:ndata]
  ID       <- rep(labcurv$ID,times=lencurv)
  timegrid <- sort(unique(times))

  ### ID, volume and time data frame

  dataset <- data.frame(ID=ID,Vol=vol,Time=times)
  alldata <- list(Dataset=dataset,LenCurv=lencurv,LabCurv=labcurv,TimeGrid=timegrid)


  return(alldata)
}
