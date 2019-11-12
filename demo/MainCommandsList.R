library(connector)

### Data files
GrowDataFile<-system.file("data", "1864dataset.xls", package = "connector")
AnnotationFile <-system.file("data", "1864info.txt", package = "connector")

### Merge curves and target file
CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)

### Visualization

GrowthCurve(CONNECTORList,"Progeny")

datavisual<-DataVisualization(CONNECTORList,feature="Progeny", labels = c("time","volume","Tumor Growth"))

### Truncation
CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",50,labels = c("time","volume","Tumor Growth"))

### Calculation of p
CrossLogLike<-BasisDimension.Choice(CONNECTORList,2:10)

CrossLogLike$CrossLogLikePlot
p<-3

### Calculation of h
pca <- PCA.Analysis(CONNECTORList,p = p)

pca$plot
h<-2

####### Calculation of k and fitting using FCM

### Stability Analysis
S.cl <-StabilityAnalysis(CONNECTORList, k =2:6, h = h, p = p ,runs=10)

### Using the Box Plots you can understand the optimal number of cluster, k.
BoxPlot.Extrapolation(stability.list = S.cl, h = h) 

# Looking at the Consensus Matrix is possible to understand how much the clustering obtained is steable.
ConsMatrix.Extrapolation(stability.list = S.cl, h = h, k = 4)

# Fixed the h and k values, here we are able to extrapolate the most probable clustering.
CONNECTORList.FCM.p3.h1.k4<-MostProbableClustering.Extrapolation(stability.list = S.cl, h = h, k = 4)

######
### Plotting discriminant functions
MaxDiscrPlots<-MaximumDiscriminationFunction(clusterdata = CONNECTORList.FCM.p3.h1.k4)

### Plotting Mean Curves and Sample Curves depending on the cluster
FCMplots<- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.p3.h1.k4, data= CONNECTORList,feature = "Progeny",labels = c("Time","Volume"),title= paste(" FCM model h=",h) )

### Disriminant Plot (goodness of the cluster) just for h = 1 or 2
DiscriminantPlot(clusterdata = CONNECTORList.FCM.p3.h1.k4, data = CONNECTORList,h = h,feature = "Progeny")

### Counting samples distribution into the clusters
NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.k4.h1,CONNECTORList,feature = "Progeny")

