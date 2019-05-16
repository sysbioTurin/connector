#' Cluster Stability
#'@description
#'
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param k  The vector/number of clusters.
#' @param h The  vector/number representing the dimension of the cluster mean space. If NULL, ClusterChoice set the $h$ value equals to the number of PCA components needed to explain the 95\% of variability of the natural cubic spline coefficients, but the \emph{PCAperc} is needed (see \code{\link{PCA.Analysis}}).
#' @param p The dimension of the natural cubic spline basis.
#' @param runs Number of runs.
#' @param seed Seed for the kmeans function.
#' @param save If TRUE then the growth curves plot truncated at the ``truncTime'' is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#'  @return
#' 
#' @examples
#'
#' @import ggplot2 reshape2 RColorBrewer statmod
#' @export

StabilityAnalysis<-function(data,k,h,p,runs=50,seed=NULL,save=FALSE,path=NULL)
{
  ALL.runs<-list()
  ConsensusInfo<-list()
  
  if( !is.null(seed)) set.seed(seed)

  ALL.runs<-lapply(1:runs, function(i) ClusterChoice(CONNECTORList, k = k, h = h, p = p,seed=NULL) )

###### box plot generation #####
  Box.pl<-list()
  
  for(hind in 1:length(h) )
  {
    pl<-list()
    
    l.tight<-lapply(1:runs,function(x) ALL.runs[[x]]$Tight.indexes[,hind])
    l.DB<-lapply(1:runs,function(x) ALL.runs[[x]]$DB.indexes[,hind])
    
    dt.fr<-do.call(cbind, l.tight)
    if(length(dt.fr[,1])==1) row.names(dt.fr)=paste("k=",k) 
    
    dt.fr<-data.frame(clust=rep(row.names(dt.fr),length(dt.fr[1,])),y=round(c(dt.fr),digits = 3) )
    counts<-count(dt.fr,vars=c("clust","y"))
    
    pl[["Tight"]]<-ggplot()+geom_boxplot(data= dt.fr,aes(x=clust,y=y))+geom_point(data=counts, col="red",aes(x=clust,y=y,size=freq/runs) )+
      labs(title=paste("Elbow method with h=",h[hind] ),x="Cluster",y="Tightness")+
      theme(text = element_text(size=20)) +labs(size="Counts freq.") 
    
    dt.fr2<-do.call(cbind, l.DB)
    if(length(dt.fr2[,1])==1) row.names(dt.fr2)=paste("k=",k) 
    
    dt.fr2<-data.frame(clust=rep(row.names(dt.fr2),length(dt.fr2[1,])),y=round(c(dt.fr2),digits = 3) )
    counts2<-count(dt.fr2,vars=c("clust","y"))
    
    pl[["DBindex"]]<-ggplot()+geom_boxplot(data= dt.fr2,aes(x=clust,y=y))+geom_point(data=counts2, col="red",aes(x=clust,y=y,size=freq/runs) )+
      labs(title=paste("DB indexes variability with h=",h[hind] ),x="Cluster",y="DB index",col="")+
      theme(text = element_text(size=20))+labs(size="Counts freq.") 
    
    Box.pl[[paste("h= ",h[hind])]]$Boxplot<-pl
    Box.pl[[paste("h= ",h[hind])]]$Data<-list(Tight=dt.fr,DBindexes=dt.fr2)
    
    ######### Cluster Membership #########
  
    Consensus.Info<-lapply(k, function(kind) {
      
      ############# first we found the most probably clustering:
      ClustCounting<-sapply(1:runs,function(x) names(ALL.runs[[x]]$FCM_all[[paste("k=",kind)]][[paste("h=",h[hind])]]$FCM$cluster$cluster.member ) )
      
      ClustString<-sapply(1:runs,function(x) paste ( table(ClustCounting[,x]) , collapse = "") )
      IndexBestClustering<-which(ClustString==names(which.max(table(ClustString))))[1]
      
      BestClustering<-ALL.runs[[IndexBestClustering]]$FCM_all[[paste("k=",kind)]][[paste("h=",h[hind])]]
      ##########################################################
      #### Let build the consensus matrix
      
      a<-lapply(1:runs, function(x){
          consensusM<-diag(0,length(data$LenCurv))
          fcm.k<-ALL.runs[[x]]$FCM_all[[paste("k=",kind)]][[paste("h=",h[hind])]]
          cl.vector<-fcm.k$FCM$cluster$cluster.member
          for(i in 1:length(data$LenCurv)) consensusM[i,which(cl.vector%in%cl.vector[i])]<-1
          consensusM
          })
      
      consensusM<-Reduce("+",a)
      
      colnames(consensusM)<- data$LabCurv$SampleName
      row.names(consensusM)<- data$LabCurv$SampleName
      
      consensusM<-as.data.frame(consensusM)/runs
      
################
### ordering the samples to have all curves in the same cluster close with each otehrs.
      
     cl.memer<-ALL.runs[[1]]$FCM_all[[paste("k=",kind)]][[paste("h=",h[hind])]]$FCM$cluster$cluster.m
      #  mat[do.call(order, as.data.frame(mat)),]
      
      fcm<-ALL.runs[[1]]$FCM_all[[paste("k=",kind)]][[paste("h=",h[hind])]]
      curvename<-data$LabCurv$SampleName
      curve<-fcm$FCM$prediction$gpred
      rownames(curve)<-curvename
      
      grid<-fcm$FCM$fit$grid
      
      gauss.quad(10) -> gauss
      
      a <- min(grid)
      b <- max(grid)
      itempi <- (a+b)/2 + (b-a)/2*gauss$nodes
      
      match(itempi,grid) -> itimeindex 
      
      fxk <- (curve[,itimeindex] )^2
      int <- (b-a)/2 * rowSums( gauss$weights * fxk )
      dist.curve <- sqrt(int)
      
      ################## Sorting the names!!!
      
      #1) sorting depending by the l2 distance with the zero x-axis!
      curvename.ordered<-names(sort(dist.curve)) 
      
      #2) Sorting depending on the cluster membership
      names(cl.memer)<-curvename
      
      #3) Defining a dataframe in order to sort the curves depending by the cluster and distance!
      
      namefroml2<-1:length(curvename)
      names(namefroml2)<-curvename.ordered
      namefroml2<-namefroml2[names(cl.memer)]
      
      curve.name.dtfr<-data.frame(namefroml2=namefroml2,cl.mem=cl.memer,1:length(curvename.ordered))
      
      ind<-curve.name.dtfr[order(curve.name.dtfr$cl.mem,curve.name.dtfr$namefroml2,curve.name.dtfr$cl.mem),3]
      
      curvename.ordered<-names(cl.memer)[ind]
     
      consensusM <- consensusM[curvename.ordered,curvename.ordered]
      consensusM <- as.matrix(consensusM)
      consensusM[upper.tri(consensusM)] <- NA
      
      m<-melt( as.matrix(consensusM) , na.rm = TRUE)
      
      m$Var2<-factor(m$Var2,levels = curvename.ordered  )
      m$Var1<-factor(m$Var1,levels = rev(curvename.ordered )  )
      
############
      
      ConsensusPlot<- ggplot(m,aes(x = Var1, y = Var2,fill=value)) + 
        geom_tile()+labs(x="", y="",fill="Same cluster \ncounting") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1))
      

      return(list(ConsensusMatrix=consensusM,ConsensusPlot=ConsensusPlot,MostProbabilyClustering=BestClustering) )
      })
    
    names(Consensus.Info) <- c(paste("k=",k))
    ConsensusInfo[[paste("h=",h[hind])]]<-Consensus.Info
    
    #####################################
    if(save==TRUE)
    {
      if(is.null(path))
      {
        path <- getwd()
      }
      ggsave(filename=paste("BoxPlotDBindexes_H",hind,".pdf"),plot =pl$DBindex,width=29, height = 20, units = "cm",scale = 1,path=path )
      ggsave(filename=paste("BoxPlotElbowMethod_H",hind,".pdf"),plot =pl$Tight,width=29, height = 20, units = "cm",scale = 1,path=path )
      for(kind in k)
      {
        ggsave(filename=paste("ConsensusMatrix_H",hind,"K",kind,".pdf",sep=""),plot = Consensus.Info$ConsensusPlot,width=29, height = 20, units = "cm",scale = 1,path=path )
      }
      
    }

  }
  
  Box.pl$seed <- seed
  
  
  return( list(ConsensusInfo=ConsensusInfo,BoxPlots=Box.pl) )
}
