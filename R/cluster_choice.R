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
#' @param seed Seed for the kmeans function.
#' @param save If TRUE then the growth curves plot truncated at the ``truncTime'' is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#'  @return
#' ClusterChoice returns the matrices of the AIC and BIC values, a list of FCMList objects belonging to class funcyOutList (see \code{\link[funcy]{funcyOutList-class}}) for each \emph{h} and \emph{k}, the Elbow Method plot and the matrix containing the total withinness measures. The distance used to calculate the two last objects is the L2 distance.
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
#' @import  ggplot2 flexclust Matrix splines statmod
#' @export
#' 
ClusterChoice<-function(data,k,h=1,p=5,PCAperc=NULL,seed=NULL,tol = 0.001, maxit = 20,save=FALSE,path=NULL)
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
  DB.indexes<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))
  Tight.indexes<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))
  
  # return a list of K lists, in which is is stored the output for all h
  # We also create two matrixes with the BIC and AIC values
  
  grid <- data$TimeGrid
  
  ############### Calculation and integration of the Gauss points into the timegrid 
  ############### for calculating the distances between curves
  ## we need the splines values calculated in the gauss$nodes in the time interval [a,b]
  
  gauss.quad(10) -> gauss

  a <- min(grid)
  b <- max(grid)
  itempi <- (a+b)/2 + (b-a)/2*gauss$nodes
  
  grid <- sort(c(grid,itempi))
  match(itempi,grid) -> itimeindex 
  
  gauss.info<-list(gauss=gauss,itimeindex=itimeindex,a=a,b=b)
  ##########################################
  
  ## Let define the input parameters for the FCM function
  
  points<-database$Vol
  ID<-database$ID
  times<-database$Time
  match(times,grid) -> timeindex
 
  for(k in K)
  {
    output_h<-list()
    for(h in H)
    {
      out.funcit<-list()
      
      fcm.fit<- fitfclust(x=points,curve=ID,timeindex=timeindex,q=p,h=h,K=k,p=p,grid=grid,seed=seed,tol = tol, maxit = maxit)
      
      fclust.curvepred(fcm.fit) -> fcm.prediction
      fclust.pred(fcm.fit) -> fcm.PRED
      
      fcm.PRED$class.pred -> cluster
      
      ClustCurve <- data.frame(ID=database[,1],Times=database[,3],Vol=database[,2],Cluster= rep(cluster,data$LenCurv))
      
      
     ##################  Goodness coefficents calculation ############
      

      fcm.prediction$meancurves->meancurves
      
      distances <- L2dist.curve2mu(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 0)
      distances.zero<-L2dist.mu20(fcm.prediction,gauss.info,fcm.fit,deriv=0)
      
      Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 0)
      Deriv.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 1)
      Deriv2.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 2)
      
      ######### Let name the cluster with A->Z from the lower mean curve to the higher.

      if( k !=1 )
      {
        Cl.order<-order(distances.zero)
      }else{
        Cl.order<-1
      }
      
      cl.names<-LETTERS[order(Cl.order)]
      
      names(cluster)<-cl.names[cluster]
      
      ################
      
      out.funcit$cluster <- list(ClustCurve=ClustCurve,cl.info=fcm.PRED,cluster.member=cluster,cluster.names=cl.names)
      out.funcit$fit <- fcm.fit
      out.funcit$prediction <- fcm.prediction
      out.funcit$TimeGrid<-data$TimeGrid
      
      output_h[[paste("h=",h)]]$FCM <- out.funcit
      output_h[[paste("h=",h)]]$Cl.Info<- list(Coefficents=Coefficents,Deriv.Coefficents=Deriv.Coefficents,Deriv2.Coefficents=Deriv2.Coefficents)
        
      DB.indexes[which(K==k),which(H==h)]<-Coefficents$DB.index
      Tight.indexes[which(K==k),which(H==h)]<-sum(distances)
      
    }

      output_k[[paste("k=",k)]]<-output_h
  }
  
  ####### Elbow method with L2 distances
  ############
  
 L2dist <- data.frame(dist=c(Tight.indexes),k=rep(K,length(H)),h=factor(rep(H,each=length(K))))
 
  ElbowMethod<-ggplot(data=L2dist,aes(x=k))+ geom_point(aes(y=dist,col=h))+
    geom_line(aes(y=dist,col=h))+
    labs(title="Elbow method ",x="Cluster",y="Tightness")+
    theme(text = element_text(size=20))
  
  ####### DB indexes plot
  ############
  
  DB.indexes.4plot <- data.frame(dist=c(DB.indexes),k=rep(K,length(H)),h=factor(rep(H,each=length(K))))
  
  DBplot <-ggplot(data=DB.indexes.4plot,aes(x=k))+ geom_point(aes(y=dist,col=h))+
    geom_line(aes(y=dist,col=h))+
    labs(title="Elbow method ",x="Cluster",y="Tightness")+
    theme(text = element_text(size=20))
  
  
  if (h==1)
  {
    DBplot<-DBplot+ theme(legend.position="none")
    ElbowMethod<-ElbowMethod+ theme(legend.position="none")
  }
  ####################
  
  if(save==TRUE)
  {
    if(is.null(path))
    {
      path <- getwd()
    }
    ggsave(filename="ElbowMethod.pdf",plot =ElbowMethod,width=29, height = 20, units = "cm",scale = 1,path=path )
    ggsave(filename="DBplot.pdf",plot = DBplot, width=29, height = 20, units = "cm",scale = 1,path=path )
  }
  
  return(list(FCM_all=output_k,matrix_BIC=matrix_BIC,matrix_AIC=matrix_AIC,ElbowMethod=ElbowMethod,DBplot=DBplot,DB.indexes=DB.indexes,Tight.indexes=Tight.indexes))
}
