library(connector)

### Data files
GrowDataFile<-system.file("data", "475dataset.xlsx", package = "connector")
AnnotationFile <-system.file("data", "475info.txt", package = "connector")

### Merge curves and target file
CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)

### Visualization

GrowPlot<-GrowthCurve(CONNECTORList,"Progeny")

Datavisual<-DataVisualization(CONNECTORList,feature="Progeny", labels = c("Time","Volume","Tumor Growth"))

Datavisual

# trunc = 50

### Truncation
CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",50,labels = c("Time","Volume","Tumor Growth"))

### Calculation of p
CrossLogLike<-BasisDimension.Choice(CONNECTORList,2:6)

CrossLogLike$CrossLogLikePlot

# p is 
p<-3

### Calculation of h
pca <- PCA.Analysis(CONNECTORList,p = p)

pca$plot

# h is 
h<-1

####### Calculation of G and fitting using FCM

### Stability Analysis
S.cl <-StabilityAnalysis(CONNECTORList, G = 2:6, h = h, p = p, runs=10)

### Using the Box Plots you can understand the optimal number of cluster, G.
BoxPlot.Extrapolation(stability.list = S.cl, h = h)

# Both G = 4 or 5 are characterized by a low fDB index and a small variation in the Elbow plot

# Looking at the Consensus Matrix is possible to understand how much the clustering obtained is steable.
ConsMatrix.Extrapolation(stability.list = S.cl, h = h, G = 4)

ConsMatrix.Extrapolation(stability.list = S.cl, h = h, G = 5)

G=4

# Fixed the h and k values, here we are able to extrapolate the most probable clustering.
CONNECTORList.FCM.p3.h1.G4 <-MostProbableClustering.Extrapolation(stability.list = S.cl, h = h, G = G)

### Plotting Mean Curves and Sample Curves depending on the cluster

FCMplots<- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.p3.h1.G4, data= CONNECTORList,feature = "Progeny",labels = c("Time","Volume"),title= paste(" FCM model h=",h))

Spline.plots(FCMplots,All=T,path="~/Desktop/")

### Disriminant Plot (goodness of the cluster) just for h = 1 or 2
DiscriminantPlot(clusterdata = CONNECTORList.FCM.p3.h1.G4, data = CONNECTORList,h = h,feature = "Progeny")


### Counting samples distribution into the clusters
NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.p3.h1.G4,CONNECTORList,feature = "Progeny")


######
# Advanced Analysis
#######
### Plotting discriminant functions
######
MaxDiscrPlots<-MaximumDiscriminationFunction(clusterdata = CONNECTORList.FCM.p3.h1.G4)

MaxDiscrPlots
