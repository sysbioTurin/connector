library(connector)

### Data files
GrowDataFile<-system.file("data", "1864dataset.xls", package = "connector")
AnnotationFile <-system.file("data", "1864info.txt", package = "connector")

### Merge curves and target file
CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)

### Visualization

GrowthCurve(CONNECTORList,"Progeny")

DataVisualization(CONNECTORList,feature="Progeny",labels = c("time","volume","Tumor Growth"))

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

CONNECTORList.FCM <- ClusterChoice(CONNECTORList, k = 2:7, h = 1:2, p = p,seed=2404)

CONNECTORList.FCM <- ClusterChoice(CONNECTORList, k = 2:6, PCAperc = pca$perc, p = p)

CONNECTORList.FCM$ElbowMethod

k<-4

####### Stability Analysis

S.cl <-StabilityAnalysis(CONNECTORList, k = 2:7, h = 1:2, p = p ,runs=50)

#######

CONNECTORList.FCM.p3.k4.h1<- CONNECTORList.FCM$FCM_all[[paste("k=",4)]][[paste("h=",1)]][["FCM"]]

### Plotting discriminant functions

MaxDiscrPlots<-MaximumDiscriminationFunction(clusterdata = CONNECTORList.FCM.p3.k4.h1)

### Plotting Mean Curves and Sample Curves depending on the cluster

FCMplots<- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.p3.k4.h1, data= CONNECTORList,feature = "Progeny",labels = c("Time","Volume"),title= " FCM model h=1 ")


FCMplots$plotMeanCurve
FCMplots$plotsCluster$ALL

### Disriminant Plot (goodness of the cluster) just for h = 1 or 2

DiscriminantPlot(clusterdata = CONNECTORList.FCM.p3.k4.h1, data= CONNECTORList,h= 1,feature="Progeny")

### Counting samples distribution into the clusters

NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.k4.h1,CONNECTORList,feature = "Progeny")


############################## MALTHUS ##############################
lower<-c(10^(-5),0)
upper<-c(10^2,10^3)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=1)

#The following function can require some  minutes to finish
Malthus1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",fitting.method="optimr",lower=lower,upper=upper,init=init)
#The following function can require some  minutes to finish
Malthus2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",fitting.method="GenSA",lower=lower,upper=upper,init=init)
#The following function can require some  minutes to finish
Malthus3<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",fitting.method="DEoptim",lower=lower,upper=upper)

CONNECTORList.Malthus1<-Malthus1$clusterList$`K= 4`
CONNECTORList.Malthus2<-Malthus2$clusterList$`K= 4`
CONNECTORList.Malthus3<-Malthus3$clusterList$`K= 4`

MalthusPlots1<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Malthus1,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "Optimr Malthus model")
MalthusPlots2<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Malthus2,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "GenSA Malthus model")
MalthusPlots3<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Malthus3,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "DEoptim Malthus model")

MalthusPlots1$plotsCluster$ALL
MalthusPlots2$plotsCluster$ALL
MalthusPlots3$plotsCluster$ALL


############################## LOGISTIC ##############################

lower<-c(10^(-5),0,0)
upper<-c(10^2,10^5,1)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=.5, b=.5)

#The following function can require some  minutes to finish
Logistic1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",fitting.method="optimr",lower=lower,upper=upper,init=init)
#The following function can require some  minutes to finish
Logistic2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",fitting.method="GenSA",lower=lower,upper=upper,init=init)
#The following function can require some  minutes to finish
Logistic3<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",fitting.method="DEoptim",lower=lower,upper=upper)

CONNECTORList.Logistic1<-Logistic1$clusterList$`K= 4`
CONNECTORList.Logistic2<-Logistic2$clusterList$`K= 4`
CONNECTORList.Logistic3<-Logistic3$clusterList$`K= 4`

LogisticPlots1<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Logistic1,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "Optimr Logistic model")
LogisticPlots2<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Logistic2,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "GenSA Logistic model")
LogisticPlots3<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Logistic3,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "DEoptim Logistic model")

LogisticPlots1$plotsCluster$ALL
LogisticPlots2$plotsCluster$ALL
LogisticPlots3$plotsCluster$ALL


############################## GOMPERTZ ##############################

lower<-c(10,0,10^(-4))
upper<-c(10^2,2,2)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=.5, b=.5)

#The following function can require some  minutes to finish
Gompertz1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",fitting.method="optimr",lower=lower,upper=upper,init=init)
#The following function can require some  minutes to finish
Gompertz2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",fitting.method="GenSA",lower=lower,upper=upper,init=init)
#The following function can require some  minutes to finish
Gompertz3<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",fitting.method="DEoptim",lower=lower,upper=upper)

CONNECTORList.Gompertz1<-Gompertz1$clusterList$`K= 4`
CONNECTORList.Gompertz2<-Gompertz2$clusterList$`K= 4`
CONNECTORList.Gompertz3<-Gompertz3$clusterList$`K= 4`

GompertzPlots1<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Gompertz1,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "Optimr Gompertz model")
GompertzPlots2<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Gompertz2,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "GenSA Gompertz model")
GompertzPlots3<-ClusterWithMeanCurve(clusterdata=CONNECTORList.Gompertz3,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "DEoptim Gompertz model")

GompertzPlots1$plotsCluster$ALL
GompertzPlots2$plotsCluster$ALL
GompertzPlots3$plotsCluster$ALL

#################################################

###### All meancurves together

pl1<-list(FCMplots$plotMeanCurve,MalthusPlots1$plotMeanCurve,LogisticPlots1$plotMeanCurve,GompertzPlots1$plotMeanCurve)

pl2<-list(FCMplots$plotMeanCurve,MalthusPlots2$plotMeanCurve,LogisticPlots2$plotMeanCurve,GompertzPlots2$plotMeanCurve)

pl3<-list(FCMplots$plotMeanCurve,MalthusPlots3$plotMeanCurve,LogisticPlots3$plotMeanCurve,GompertzPlots3$plotMeanCurve)



############## Counting the samples

NumberSamples<-CountingSamples(clusterdata=Malthus1,CONNECTORList,feature = "Progeny")

NumberSamples<-CountingSamples(clusterdata=Logistic1,CONNECTORList,feature = "Progeny")

NumberSamples<-CountingSamples(clusterdata=Gompertz1,CONNECTORList,feature = "Progeny")
