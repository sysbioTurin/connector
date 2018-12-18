#' Cluster Choice
#'
#'@description
#'
#'  Fits and clusters the data with respect to the Functional Clustering Model [Sugar and James]. The BIC and AIC values considering k number of clusters and h dimension of the cluster mean space are calculated, and the plot based on the Elbow Method is generated. As explained in [Sugar and James], to have a simple low-dimensional representation of the individual curves and to reduce the number of parameters to be estimated, h value must be equals or lower than \eqn{\min(p,k-1)}.
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param k  The vector/number of clusters.
#' @param h The  vector/number representing the dimension of the cluster mean space. If NULL, ClusterChoice set the $h$ value equals to the number of PCA components needed to explain the 95\% of variability of the natural cubic spline coefficients, but the \emph{PCAperc} is needed (see \code{\link{PCA.Analysis}}).
#' @param p The dimension of the natural cubic spline basis.
#' @param PCAperc The PCA percentages applied to the natural cubic spline coefficients, if  NULL then $h$ is needed (see \code{\link{PCA.Analysis}}).
#'  @return
#' ClusterChoice returns the matrices of the AIC and BIC values, a list of FCMList objects belonging to class funcyOutList (see \code{\link[funcy]{funcyOutList-class}}) for each \emph{h} and \emph{k}, the Elbow Method plot and the matrix containing the total withinness measures. The distance used to calculate the two last objects is the Euclidian distance.
#' 
#' @seealso \code{\link[funcy]{funcit}}.
#' 
#' @references
#' Gareth M. James and Catherine A. Sugar, (2003). Clustering for Sparsely Sampled Functional Data. Journal of the American Statistical Association.
#' 
#' @examples
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))
#'
#'
#'### Calculation of k and fitting using FCM
#' # Specifying the h value
#' 
#'CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=c(1:2))
#'
#' # Using the PCA percentaes instead of the h value.
#' pca <- PCA.Analysis(CONNECTORList)
#' CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),PCAperc = pca$perc)
#'
#'
#'
#' @import  ggplot2 flexclust
#' @export
#' 
ClusterChoice<-function(data,k,h=1,p=5,PCAperc=NULL)
{
  K<-k
 
  database<-data$Dataset
  

  if(!is.null(PCAperc)){ h<-min(which(cumsum(PCAperc)>=95)) }
 
  H<-h
  
  output_k<-list()
  row_names <-c(paste("k=",K))
  col_names<-c(paste("h=",H))
  matrix_AIC<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))
  matrix_BIC<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))
  
  
  data.funcit <-matrix(c(database$ID,database$Vol,database$Time),ncol=3,byrow=F)
########## 
  Tot.within.Eucl<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))

  # return a list of K lists, in which is is stored the output for all h
  # We also create two matrixes with the BIC and AIC values
  for(k in K)
  {
    output_h<-list()
    for(h in H)
    {
      mycontfclust = new("funcyCtrl",baseType="splines",dimBase=p,init="kmeans",nrep=10,redDim=h)
      
      out.funcit<- funcit.simo(data.funcit,seed=2404,k,methods="fitfclust",funcyCtrl=mycontfclust,save.data=TRUE)
      
      Cluster(out.funcit)->cluster
      if(out.funcit@reg){
      ClustCurve <- data.frame(ID=database[,1],Times=database[,3],Vol=database[,2],Cluster= rep(cluster,out.funcit@timeNr[,1]))
      }else{
        ClustCurve <- data.frame(ID=database[,1],Times=database[,3],Vol=database[,2],Cluster= rep(cluster,out.funcit@timeNr[,2]))
      }
      out.funcit@models$fitfclust@fit[["ClustCurves"]] <- ClustCurve
      
      output_h[[paste("h=",h)]] <- out.funcit
 
      ##################     
      if(out.funcit@reg==1)
      {
        fitfclust.curvepred(out.funcit@models$fitfclust@fit)$meancurves->meancurves
      } else{
        fitfclust.curvepredIrreg(out.funcit@models$fitfclust@fit)$meancurves->meancurves
      }
      dist<-dist2centers.new(out.funcit@data, meancurves)
      minDist<-makeClMat(dist)[,1]
      Tot.within.Eucl[which(K==k),which(H==h)]<-sum(minDist)
      ##################     
      
      matrix_BIC[which(K==k),which(H==h)]<-output_h[[paste("h=",h)]]@models$fitfclust@BIC
      matrix_AIC[which(K==k),which(H==h)]<-output_h[[paste("h=",h)]]@models$fitfclust@AIC
    }
    output_k[[paste("k=",k)]]<-output_h
  }
  
  ####### Elbow method with Hausdorff distance
  ############
  Tot.within.Eucl<-data.frame(dist=c(Tot.within.Eucl),K=rep(K,length(H)),H=factor(rep(H,each=length(K))))
 
ElbowMethod.Eucl<-ggplot(data=Tot.within.Eucl,aes(x=K))+ geom_point(aes(y=dist,col=H))+
    geom_line(aes(y=dist,col=H))+
    labs(title="Elbow method",x="Cluster",y="total within-cluster")+
    theme(text = element_text(size=20))
  
  return(list(FCM_all=output_k,matrix_BIC=matrix_BIC,matrix_AIC=matrix_AIC,ElbowMethod=ElbowMethod.Eucl,Tot.within=Tot.within.Eucl))
}
