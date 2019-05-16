library(connector)

### Data files
GrowDataFile<-system.file("data", "1864dataset.xls", package = "connector")
AnnotationFile <-system.file("data", "1864info.txt", package = "connector")

### Merge curves and target file
CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)

### Visualization

gr<-GrowthCurve(CONNECTORList,"Progeny")

datavisual<-DataVisualization(CONNECTORList,feature="Progeny", labels = c("time","volume","Tumor Growth"))

### Truncation

CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",60,labels = c("time","volume","Tumor Growth"))

### Calculation of p

CrossLogLike<-BasisDimension.Choice(CONNECTORList,2:10)

CrossLogLike$CrossLogLikePlot

p<-3
### Calculation of h

pca <- PCA.Analysis(CONNECTORList,p = p)

pca$plot

h<-1

### Calculation of k and fitting using FCM

CONNECTORList.FCM <- ClusterChoice(CONNECTORList, k = 2:7, h = 1, p = 3,seed=2404)

CONNECTORList.FCM <- ClusterChoice(CONNECTORList, k = 2:6, PCAperc = pca$perc, p = p)

CONNECTORList.FCM$ElbowMethod
CONNECTORList.FCM$DBplot
k<-4

####### Stability Analysis

S.cl <-StabilityAnalysis(CONNECTORList, k =2:6, h = 1, p = 3 ,runs=50)

S.cl$ConsensusInfo$`h= 1`$`k= 3`$ConsensusPlot
S.cl$BoxPlots$`h=  1`$Boxplot$Tight

CONNECTORList.FCM.p3.k4.h1<-S.cl$ConsensusInfo$`h= 1`$`k= 4`$MostProbabilyClustering

#######

CONNECTORList.FCM.p3.k4.h1<- CONNECTORList.FCM$FCM_all[[paste("k=",4)]][[paste("h=",1)]]

### Plotting discriminant functions

MaxDiscrPlots<-MaximumDiscriminationFunction(clusterdata = CONNECTORList.FCM.p3.k4.h1$FCM)

### Plotting Mean Curves and Sample Curves depending on the cluster

FCMplots<- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.p3.k4.h1, data= CONNECTORList,feature = "Progeny",labels = c("Time","Volume"),title= " FCM model h=1 ")

FCMplots$plotsCluster$ALL
FCMplots$plotMeanCurve
FCMplots$plotsCluster$ALL

### Disriminant Plot (goodness of the cluster) just for h = 1 or 2

DiscriminantPlot(clusterdata = CONNECTORList.FCM.p3.k4.h1, data= CONNECTORList,h= 1,feature="Progeny")

### Counting samples distribution into the clusters

NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.k4.h1$FCM,CONNECTORList,feature = "Progeny")

