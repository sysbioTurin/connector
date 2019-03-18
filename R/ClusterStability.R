#' Cluster Stability
#'
#'@description
#'
#' @examples
#'
#' @import ggplot2
#' @export

StabilityAnalysis<-function(data,k,h,p,runs=50)
{
  ALL.runs<-list()
  seed<-.Random.seed[1:runs]

  ALL.runs<-lapply(1:runs, function(i) ClusterChoice(CONNECTORList, k = k, h = h, p = p,seed=seed[i]) )

###### box plot generation #####
  Box.pl<-list()
  
  for(hind in h)
  {
    pl<-list()
    
    l.tight<-lapply(1:runs,function(x) ALL.runs[[x]]$Tight.indexes[,hind])
    l.DB<-lapply(1:runs,function(x) ALL.runs[[x]]$DB.indexes[,hind])
    
    dt.fr<-do.call(cbind, l.tight)
    dt.fr<-data.frame(clust=rep(row.names(dt.fr),length(dt.fr[1,])),y=c(dt.fr) )
    pl[["Tight"]]<-ggplot(dt.fr,aes(x=clust,y=y))+geom_boxplot()+
      labs(title=paste("Elbow method with h=",hind ),x="Cluster",y="Tightness")+
      theme(text = element_text(size=20))
    
    dt.fr2<-do.call(cbind, l.DB)
    dt.fr2<-data.frame(clust=rep(row.names(dt.fr2),length(dt.fr2[1,])),y=c(dt.fr2) )
    pl[["DBindex"]]<-ggplot(dt.fr2,aes(x=clust,y=y))+geom_boxplot()+
      labs(title=paste("DB indexes variability with h=",hind ),x="Cluster",y="DB index")+
      theme(text = element_text(size=20))
    
    Box.pl[[paste("h= ",hind)]]$Boxplot<-pl
    Box.pl[[paste("h= ",hind)]]$Data<-list(Tight=dt.fr,DBindexes=dt.fr2)
  
  }
  
  Box.pl$seed <- seed
  
  return(Box.pl)
}
