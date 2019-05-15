#' Choice of the B-spline dimension 
#'
#'@description
#'	Generates a line plot reporting the cross-validated loglikelihood value for each number of knots. In details, for each number of knots 10\% of the curves from the whole data are removed and treated as a test set, then the remaing curves are fitted using the FCM and the loglikelihood on the test set is calculated. The process is then repeated nine more times.
#'  
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param p The vector of the dimension of the natural cubic spline basis.
#' @param save If TRUE then the plot is saved into a pdf file.
#' @param path The folder path  where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @return
#' DimensionBasis.Choice returns line plot of the cross-validated loglikelihood for each value of p, in grey the result of all ten repetitions of the likelihood calculation and in black the mean of them.
#'
#' 
#'  
#' @references
#' Gareth M. James and Catherine A. Sugar, (2000). Principal component models for sparse functional data.
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
#'CrossLogLike<-BasisDimension.Choice(CONNECTORList,2:10)
#'CrossLogLike$CrossLogLikePlot
#' 
#'
#'
#' @import ggplot2 fda plyr
#' @importFrom MASS ginv
#' @export
#' 
BasisDimension.Choice<-function(data,p,save=FALSE,path=NULL)
{
  p.values<-p
  crossvalid<-list()
  n_sample<-length(data$LenCurv)
  perc<- as.integer(n_sample*0.1)
  
  grid<-data$TimeGrid
  
  for(step in 1:10)
  {
  
  SampleTestSet<-sample(1:n_sample,perc)
  SampleTestSet<-SampleTestSet[order(SampleTestSet)]
  
  TestSet<-CONNECTORList$Dataset[CONNECTORList$Dataset$ID%in%SampleTestSet,]
  TestSet <-data.frame(ID=rep(1:(length(SampleTestSet)),CONNECTORList$LenCurv[SampleTestSet]),Vol=TestSet$Vol[],Time=TestSet$Time)
  
  TrainingSet<-CONNECTORList$Dataset[-which(CONNECTORList$Dataset$ID%in%SampleTestSet),]

  #### Let define the FCM input matrix of ID sample, Points and time indexes per sample into the grid.
  
  data.funcit <-matrix(c(rep(1:(n_sample-length(SampleTestSet)),CONNECTORList$LenCurv[-SampleTestSet]),TrainingSet$Vol,match(TrainingSet$Time,grid)),ncol=3,byrow=F)
 
  CrossLikelihood<-sapply(p.values, CalcLikelihood, data.funcit,TestSet,grid)
  
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
  
  if(save==TRUE)
  {
    ggsave(filename="CrossLogLikePlot_Estimation_p.pdf",plot = ValidationPlot,width=29, height = 20, units = "cm",scale = 1,path=path )
  }
  
  
  return(list(CrossLogLikePlot=ValidationPlot))
}


CalcLikelihood<-function(p,data.funcit,TestSet,grid){
 
  perc<-max(TestSet[,1])
  
  #grid<-sort(unique(data.funcit[,3] ) )
  points<-data.funcit[,2]
  ID<-data.funcit[,1]
  timeindex<-data.funcit[,3]
  
  fcm.fit <- fitfclust(x=points,curve=ID,timeindex=timeindex,q=p,h=1,K=1,p=p,grid=grid,seed=2404)
  
  Gamma<-as.matrix(fcm.fit$par$Gamma)
  sigma<-fcm.fit$par$sigma
  Lambda<- fcm.fit$par$Lambda
  alpha<- fcm.fit$par$alpha
  mu<- fcm.fit$par$lambda.zero+Lambda*c(alpha)
  
  ##### make basis considering the test set
  
  # TimeGrid<- unique(sort(TestSet[,3]))
  # 
  # tempBase <-cbind(1, ns(TimeGrid, df = (p - 1)))
  # base <- svd(tempBase)$u
  
  base<-fcm.fit$FullS
  
  
  Likelihood<-function(x,base,TestSet)
  {
    data.temp<-TestSet[TestSet$ID==x,]
    time.temp<-data.temp$Time
    base.i<-base[grid%in%time.temp,]
    Yi<-data.temp$Vol
    n.i<-length(time.temp)
    Vi<-base.i%*%Gamma%*%t(base.i)+sigma^2*diag(n.i)
    mui<-base.i%*%mu
    invVi<-ginv(Vi)
    -n.i/2*log(det(Vi)) - n.i/2*log(2*pi) - 1/2*t(Yi-mui)%*%invVi%*%(Yi-mui)
  }
  
  Li<-sapply(1:perc,Likelihood,base,TestSet)
  Likelihood<-sum(Li)
  return( Likelihood )
}
