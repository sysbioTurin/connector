library(readr)
library(readxl)
library(tibble)
library(xlsx)

setwd("~/Dropbox/Universita/Progetto CONNECTOR/newdata/Mirna")

Mean<-FALSE
MEAN<-20
Rescale<-FALSE

ALL<-FALSE
G3<-FALSE
G1.G2<-TRUE

source("LetturaTabelle.R")

library(connector)

GrowthCurve(CONNECTORList,"Patient")
DataVisualization(CONNECTORList,feature="Patient",labels = c("time","volume","Mirna"))

### Calculation of p

CrossLogLike<-DimensionBasis.Choice(CONNECTORList,2:6)

CrossLogLike$CrossLogLikePlot

p<-4
### Calculation of h

pca <- PCA.Analysis(CONNECTORList,p = 5)

pca$plot

h<-2

### Calculation of k and fitting using FCM

CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:10),h=1:3,p = 4)
CONNECTORList.FCM$ElbowMethod

k<-4


CONNECTORList.FCM.k<- CONNECTORList.FCM$FCM_all[[paste("k=",5)]][[paste("h=",2)]]

DiscriminantPlot(clusterdata = CONNECTORList.FCM.k, data= CONNECTORList,h=2,feature="Patient")

FCMplots<- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.k, data= CONNECTORList,feature = "Patient",labels = c("Time","Volume"),title= " FCM model h=2 ")


FCMplots$plots$plotMeanCurve
FCMplots$plots$plotsCluster$ALL

save()
