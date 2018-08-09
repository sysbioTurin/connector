
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

### Calculation of h

pca <- PCA.Analysis(CONNECTORList,p = 5)

pca$plot

### Calculation of k and fitting using FCM

CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=2)

CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),PCAperc = pca$perc)

CONNECTORList.FCM.k4.h2<- CONNECTORList.FCM$FCM_all$`k= 4`$`h= 2`


FCMplots <- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.k4.h2,data= CONNECTORList,feature = "Progeny",labels = c("Time","Volume"),title= " FCM model ")


FCMplots$plotsCluster$ALL

######### separation and tightness plot considering the FCM

PlotSeparationTightness(CONNECTORList.FCM.k4.h2,Title = "FCM Cluster betweenness and withinness")

############################## MALTHUS ##############################
lower<-c(10^(-5),0)
upper<-c(10^2,10^3)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=1)


Malthus1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",feature="Progeny",fitting.method="optimr",lower=lower,upper=upper,init=init)
Malthus2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",feature="Progeny",fitting.method="GenSA",lower=lower,upper=upper,init=init)
Malthus3<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",feature="Progeny",fitting.method="DEoptim",lower=lower,upper=upper)

MalthusPlots1<-ClusterWithMeanCurve(clusterdata=Malthus1,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "Optimr Malthus model")
MalthusPlots2<-ClusterWithMeanCurve(clusterdata=Malthus2,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "GenSA Malthus model")
MalthusPlots3<-ClusterWithMeanCurve(clusterdata=Malthus3,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "DEoptim Malthus model")

MalthusPlots1$plotsCluster$ALL
MalthusPlots2$plotsCluster$ALL
MalthusPlots3$plotsCluster$ALL

##### separation and tightness plot considering the Malthus model

PlotSeparationTightness(Malthus1,Title = "Malthus Optimr Cluster betweenness and withinness")

PlotSeparationTightness(Malthus2,Title = "Malthus Gensa Cluster betweenness and withinness")

PlotSeparationTightness(Malthus3,Title = "Malthus Deoptim Cluster betweenness and withinness")

############################## LOGISTIC ##############################

lower<-c(10^(-5),0,0)
upper<-c(10^2,10^5,1)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=.5, b=.5)

Logistic1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",feature="Progeny",fitting.method="optimr",lower=lower,upper=upper,init=init)
Logistic2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",feature="Progeny",fitting.method="GenSA",lower=lower,upper=upper,init=init)
Logistic3<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",feature="Progeny",fitting.method="DEoptim",lower=lower,upper=upper)

LogisticPlots1<-ClusterWithMeanCurve(clusterdata=Logistic1,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "Optimr Logistic model")
LogisticPlots2<-ClusterWithMeanCurve(clusterdata=Logistic2,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "GenSA Logistic model")
LogisticPlots3<-ClusterWithMeanCurve(clusterdata=Logistic3,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "DEoptim Logistic model")

LogisticPlots1$plotsCluster$ALL
LogisticPlots2$plotsCluster$ALL
LogisticPlots3$plotsCluster$ALL


##### separation and tightness plot considering the Logistic model

PlotSeparationTightness(Logistic1,Title = "Logistic Optimr Cluster betweenness and withinness")

PlotSeparationTightness(Logistic2,Title = "Logistic Gensa Cluster betweenness and withinness")

PlotSeparationTightness(Logistic3,Title = "Logistic Deoptim Cluster betweenness and withinness")

############################## GOMPERTZ ##############################

lower<-c(10,0,10^(-4))
upper<-c(10^2,2,2)
init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=.5, b=.5)

Gompertz1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",feature="Progeny",fitting.method="optimr",lower=lower,upper=upper,init=init)
Gompertz2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",feature="Progeny",fitting.method="GenSA",lower=lower,upper=upper,init=init)
Gompertz3<- FittingAndClustering(data = CONNECTORList, k = 4, model="Gompertz",feature="Progeny",fitting.method="DEoptim",lower=lower,upper=upper)

GompertzPlots1<-ClusterWithMeanCurve(clusterdata=Gompertz1,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "Optimr Gompertz model")
GompertzPlots2<-ClusterWithMeanCurve(clusterdata=Gompertz2,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "GenSA Gompertz model")
GompertzPlots3<-ClusterWithMeanCurve(clusterdata=Gompertz3,data = CONNECTORList, feature = "Progeny",labels = c("Time","Volume"),title= "DEoptim Gompertz model")

GompertzPlots1$plotsCluster$ALL
GompertzPlots2$plotsCluster$ALL
GompertzPlots3$plotsCluster$ALL

###### separation and tightness plot considering the Gompertz model

a2<-PlotSeparationTightness(Gompertz1,Title = "Gompertz Optimr Cluster betweenness and withinness")
b2<-PlotSeparationTightness(Gompertz2,Title = "Gompertz Gensa Cluster betweenness and withinness")
c2<-PlotSeparationTightness(Gompertz3,Title = "Gompertz Deoptim Cluster betweenness and withinness")
#################################################

###### All meancurves together

pl1<-list(FCMplots$plotMeanCurve,MalthusPlots1$plotMeanCurve,LogisticPlots1$plotMeanCurve,GompertzPlots1$plotMeanCurve)

pl2<-list(FCMplots$plotMeanCurve,MalthusPlots2$plotMeanCurve,LogisticPlots2$plotMeanCurve,GompertzPlots2$plotMeanCurve)

pl3<-list(FCMplots$plotMeanCurve,MalthusPlots3$plotMeanCurve,LogisticPlots3$plotMeanCurve,GompertzPlots3$plotMeanCurve)



############## Counting the samples

NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.k4.h2,CONNECTORList,feature = "Progeny")
