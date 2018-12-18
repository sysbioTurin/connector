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

CrossLogLike<-DimensionBasis.Choice(CONNECTORList,2:10)

CrossLogLike$CrossLogLikePlot

p<-3
### Calculation of h

pca <- PCA.Analysis(CONNECTORList,p = p)

pca$plot

h<-c(1,2)

### Calculation of k and fitting using FCM

CONNECTORList.FCM <- ClusterChoice(CONNECTORList, k = c(2:6), h = c(1,2) , p = 3)

CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),PCAperc = pca$perc , p = p)

CONNECTORList.FCM$ElbowMethod

k<-4


CONNECTORList.FCM.p3.k4.h2<- CONNECTORList.FCM$FCM_all[[paste("k=",4)]][[paste("h=",1)]]


#discriminantplot(CONNECTORList.FCM.p3.k4.h2)

a<-MaximumDiscriminationFunction(clusterdata = CONNECTORList.FCM.p3.k4.h2)

FCMplots<- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.p3.k4.h2, data= CONNECTORList,feature = "Progeny",labels = c("Time","Volume"),title= " FCM model h=1 ")


FCMplots$plots$plotMeanCurve
FCMplots$plots$plotsCluster$ALL

DiscriminantPlot(clusterdata = CONNECTORList.FCM.p3.k4.h2, data= CONNECTORList,h=2,feature="Progeny")

NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.k4.h2,CONNECTORList,feature = "Progeny")

######### separation and tightness plot considering the FCM

PlotSeparationTightness(clusterdata=CONNECTORList.FCM.p3.k4.h2,Title = "FCM Cluster betweenness and withinness h=2 ",save = TRUE,path="../../Dropbox/Universita/Progetto CONNECTOR/newdata/dati1864/FCMp3h2k3/")
PlotSeparationTightness(clusterdata=CONNECTORList.FCM.p3.k4.h2,Title = "FCM Cluster betweenness and withinness h=1 ")


############################## MALTHUS ##############################
lower<-c(10^(-5),0)
upper<-c(10^2,10^3)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=1)


Malthus1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",fitting.method="optimr",lower=lower,upper=upper,init=init)
Malthus2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",fitting.method="GenSA",lower=lower,upper=upper,init=init)
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

##### separation and tightness plot considering the Malthus model

PlotSeparationTightness(CONNECTORList.Malthus1,Title = "Malthus Optimr Cluster betweenness and withinness")

PlotSeparationTightness(CONNECTORList.Malthus2,Title = "Malthus Gensa Cluster betweenness and withinness")

PlotSeparationTightness(CONNECTORList.Malthus3,Title = "Malthus Deoptim Cluster betweenness and withinness")

############################## LOGISTIC ##############################

lower<-c(10^(-5),0,0)
upper<-c(10^2,10^5,1)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=.5, b=.5)

Logistic1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",fitting.method="optimr",lower=lower,upper=upper,init=init)
Logistic2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",fitting.method="GenSA",lower=lower,upper=upper,init=init)
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


##### separation and tightness plot considering the Logistic model

PlotSeparationTightness(CONNECTORList.Logistic1,Title = "Logistic Optimr Cluster betweenness and withinness")

PlotSeparationTightness(CONNECTORList.Logistic2,Title = "Logistic Gensa Cluster betweenness and withinness")

PlotSeparationTightness(CONNECTORList.Logistic3,Title = "Logistic Deoptim Cluster betweenness and withinness")

############################## GOMPERTZ ##############################

lower<-c(10,0,10^(-4))
upper<-c(10^2,2,2)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=.5, b=.5)

Gompertz1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",fitting.method="optimr",lower=lower,upper=upper,init=init)
Gompertz2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",fitting.method="GenSA",lower=lower,upper=upper,init=init)
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

###### separation and tightness plot considering the Gompertz model

a2<-PlotSeparationTightness(CONNECTORList.Gompertz1,Title = "Gompertz Optimr Cluster betweenness and withinness")
b2<-PlotSeparationTightness(CONNECTORList.Gompertz2,Title = "Gompertz Gensa Cluster betweenness and withinness")
c2<-PlotSeparationTightness(CONNECTORList.Gompertz3,Title = "Gompertz Deoptim Cluster betweenness and withinness")
#################################################

###### All meancurves together

pl1<-list(FCMplots$plotMeanCurve,MalthusPlots1$plotMeanCurve,LogisticPlots1$plotMeanCurve,GompertzPlots1$plotMeanCurve)

pl2<-list(FCMplots$plotMeanCurve,MalthusPlots2$plotMeanCurve,LogisticPlots2$plotMeanCurve,GompertzPlots2$plotMeanCurve)

pl3<-list(FCMplots$plotMeanCurve,MalthusPlots3$plotMeanCurve,LogisticPlots3$plotMeanCurve,GompertzPlots3$plotMeanCurve)



############## Counting the samples
NumberSamples1<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.k4.h2_1,CONNECTORList,feature = "Progeny")

NumberSamples2<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.k4.h2_2,CONNECTORList,feature = "Progeny")

NumberSamples3<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.k4.h2_3,CONNECTORList,feature = "Progeny")



NumberSamples<-CountingSamples(clusterdata=Malthus1,CONNECTORList,feature = "Progeny")

NumberSamples<-CountingSamples(clusterdata=Logistic1,CONNECTORList,feature = "Progeny")

NumberSamples<-CountingSamples(clusterdata=Gompertz1,CONNECTORList,feature = "Progeny")
