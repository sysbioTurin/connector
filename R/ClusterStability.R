#' Cluster Stability
#'
#'@description
#'
#' @examples
#'
#' @import ggplot2 reshape2 RColorBrewer statmod
#' @export

StabilityAnalysis<-function(data,k,h,p,runs=50,seed=NULL)
{
  ALL.runs<-list()
  ConsensusInfo<-list()
  
  if( !is.null(seed)) set.seed(seed)

  ALL.runs<-lapply(1:runs, function(i) ClusterChoice(CONNECTORList, k = k, h = h, p = p) )

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
      
      M <- ALL.runs[[1]]$FCM_all[[paste("k=",kind)]][[paste("h=",h[hind])]]$Cl.Info$Coefficents$emme
      if( kind !=1 )
      {
        for(x in 1:length(M[1,]))
        {
          #if( !all(!M[x,]%in%min(M[M!=0]) ))
            if( !all(!M[x,]%in%max(M[M!=0]) )) order(M[x,])-> Cl.order
        }
      }else{
        Cl.order<-1
      }
      
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
      curvename.ordered<-names(sort(dist.curve))
      
      
      names(cl.memer)<-curvename
      
     # curvename.ordered<-names(cl.memer[order(match(cl.memer,Cl.order))])
     
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
      

      return(list(ConsensusMatrix=consensusM,ConsensusPlot=ConsensusPlot) )
      })
    
    names(Consensus.Info) <- c(paste("k=",k))
    ConsensusInfo[[paste("h=",h[hind])]]<-Consensus.Info
    
    #####################################
  }
  
  Box.pl$seed <- seed
  
  return( list(ConsensusInfo=ConsensusInfo,BoxPlots=Box.pl) )
}
