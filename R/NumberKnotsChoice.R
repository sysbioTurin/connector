#' Choice of the B-spline dimension 
#'
#'@description
#'
#'  
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param p The vector of the dimension of the natural cubic spline basis.
#' 
#' @return
#' DimensionBasis.Choice returns a plot with ...as a ggplot object.
#'
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
#' .................
#'
#'
#' @import ggplot2 MASS fda plyr
#' @export
#' 
DimensionBasis.Choice<-function(data,p.values)
{
  
  crossvalid<-list()
  n_sample<-length(data$LenCurv)
  perc<- as.integer(n_sample*0.1)
  
  for(step in 1:10)
  {
  
  SampleTestSet<-sample(1:n_sample,perc)
  SampleTestSet<-SampleTestSet[order(SampleTestSet)]
  
  TestSet<-CONNECTORList$Dataset[CONNECTORList$Dataset$ID%in%SampleTestSet,]
  TestSet <-data.frame(ID=rep(1:(length(SampleTestSet)),CONNECTORList$LenCurv[SampleTestSet]),Vol=TestSet$Vol[],Time=TestSet$Time)
  
  TrainingSet<-CONNECTORList$Dataset[-which(CONNECTORList$Dataset$ID%in%SampleTestSet),]
  
  data.funcit <-matrix(c(rep(1:(n_sample-length(SampleTestSet)),CONNECTORList$LenCurv[-SampleTestSet]),TrainingSet$Vol,TrainingSet$Time),ncol=3,byrow=F)
  
  CrossLikelihood<-sapply(p.values, CalcLikelihood, data.funcit,TestSet)
  
  crossvalid[[step]]<-data.frame(p=p.values,CrossLikelihood=CrossLikelihood,sim=step)
  
  }
  
  ALLcrossvalid<-ldply(crossvalid, rbind)
  mean<-sapply(p.values, function(x){ mean(ALLcrossvalid$CrossLikelihood[ALLcrossvalid$p==x]) } )
  meandata<-data.frame(p=p.values,mean=mean)
  
  ValidationPlot<-ggplot()+
                  geom_line(data=meandata,aes(x=p,y=mean,linetype="mean",col="mean"),size=1.2)+
                  geom_line(data=ALLcrossvalid,aes(x=p,y=CrossLikelihood,group=sim,linetype="test",col="test"),size=1.1)+
                  geom_point(data=meandata,aes(x=p,y=mean),size=2 )+
                  ylab(" Cross-LogLikelihood ")+xlab("number of knots")+
                  scale_color_manual("",breaks=c("mean","test"),values = c("black","grey"))+
                  scale_linetype_manual("",breaks=c("mean","test"),values = c("solid","dashed"))+
                  theme(legend.title=element_blank())
  
  return(list(CrossLogLikePlot=ValidationPlot,Meanvalues=meandata,CrossLogLike=ALLcrossvalid))
}


CalcLikelihood<-function(p,data.funcit,TestSet){
 
  perc<-max(TestSet[,1])
  
  mycontfclust = new("funcyCtrl",baseType="splines",dimBase=p,init="kmeans",nrep=10,redDim=1)
  
  out.funcit<- funcit.simo(data.funcit,seed=2404,1,methods="fitfclust",funcyCtrl=mycontfclust,save.data=TRUE)
  
  Gamma<-as.matrix(out.funcit@models$fitfclust@prms$Gamma)
  sigma<-out.funcit@models$fitfclust@prms$sigma
  Lambda<- out.funcit@models$fitfclust@prms$Lambda
  alpha<- out.funcit@models$fitfclust@prms$alpha
  mu<- out.funcit@models$fitfclust@prms$lambda.zero+Lambda*c(alpha)
  
  ##### make basis considering the test set
  
  TimeGrid<- unique(sort(TestSet[,3]))
  
  
  m <- length(TimeGrid)
  
  bObj<-create.bspline.irregular(c(TimeGrid[1],TimeGrid[m]),
                           nbasis=p,
                           norder=min(p, 4))
  base <- eval.basis(TimeGrid, bObj)
  

  
  Likelihood<-function(x)
  {
    data.temp<-TestSet[TestSet$ID==x,]
    time.temp<-data.temp$Time
    base.i<-base[TimeGrid%in%time.temp,]
    Yi<-data.temp$Vol
    n.i<-length(time.temp)
    Vi<-base.i%*%Gamma%*%t(base.i)+sigma^2*diag(n.i)
    mui<-base.i%*%mu
    invVi<-ginv(Vi)
    -n.i/2*log(det(Vi)) - n.i/2*log(2*pi) - 1/2*t(Yi-mui)%*%invVi%*%(Yi-mui)
  }
  
  Li<-sapply(1:perc,Likelihood)
  Likelihood<-sum(Li)
  return( Likelihood )
}
