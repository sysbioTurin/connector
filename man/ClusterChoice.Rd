% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_choice.R
\name{ClusterChoice}
\alias{ClusterChoice}
\title{Cluster Choice}
\usage{
ClusterChoice(data, G, h = 1, p = 5, PCAperc = NULL, seed = 2404,
  tol = 0.001, maxit = 20, save = FALSE, path = NULL)
}
\arguments{
\item{data}{CONNECTORList. (see \code{\link{DataImport}})}

\item{G}{The vector/number of clusters.}

\item{h}{The  vector/number representing the dimension of the cluster mean space. If NULL, ClusterChoice set the $h$ value equals to the number of PCA components needed to explain the 95\% of variability of the natural cubic spline coefficients, but the \emph{PCAperc} is needed (see \code{\link{PCA.Analysis}}).}

\item{p}{The dimension of the natural cubic spline basis.}

\item{PCAperc}{The PCA percentages applied to the natural cubic spline coefficients, if  NULL then $h$ is needed (see \code{\link{PCA.Analysis}}).}

\item{seed}{Seed for the kmeans function.}

\item{save}{If TRUE then the growth curves plot truncated at the ``truncTime'' is saved into a pdf file.}

\item{path}{The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
 @return
ClusterChoice returns the matrices of the AIC and BIC values, a list of FCMList objects, the Elbow Method plot and the matrix containing the total withinness measures. The distance used to calculate the two last objects is the L2 distance.}
}
\description{
Fits and clusters the data with respect to the Functional Clustering Model [Sugar and James]. The BIC and AIC values considering G number of clusters and h dimension of the cluster mean space are calculated, and the plot based on the Elbow Method is generated. As explained in [Sugar and James], to have a simple low-dimensional representation of the individual curves and to reduce the number of parameters to be estimated, h value must be equals or lower than \eqn{\min(p,G-1)}.
}
\examples{

GrowDataFile<-"data/1864dataset.xls"
AnnotationFile <-"data/1864info.txt"

CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)

CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))


### Calculation of k and fitting using FCM
# Specifying the h value

CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=c(1:2))

# Using the PCA percentaes instead of the h value.
pca <- PCA.Analysis(CONNECTORList)
CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),PCAperc = pca$perc)



}
\references{
Gareth M. James and Catherine A. Sugar, (2003). Clustering for Sparsely Sampled Functional Data. Journal of the American Statistical Association.
}
\author{
Cordero Francesca, Pernice Simone, Sirovich Roberta
}
