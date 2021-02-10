#' Cluster Stability
#' 
#'@description
#'
#'  Fits and clusters the data with respect to the Functional Clustering Model [Sugar and James]. The plots based on the Elbow Method and on the functional Davies Bouldin (fDB) indexes are generated to properly guide the choice of the number of clusters. As explained in [Sugar and James], to have a simple low-dimensional representation of the individual curves and to reduce the number of parameters to be estimated, h value must be equals or lower than \eqn{min(p,G-1)}.
#'  
#' @param data CONNECTORList. (see \code{\link{DataImport}} or \code{\link{DataTruncation}})
#' @param G  The vector/number of possible clusters.
#' @param h  The  vector/number representing the dimension of the cluster mean space. If NULL, ClusterChoice set the h value equals to the number of PCA components needed to explain the 95\% of variability of the natural cubic spline coefficients, but the \emph{PCAperc} is needed (see \code{\link{PCA.Analysis}}).
#' @param p The dimension of the natural cubic spline basis. (see \code{\link{BasisDimension.Choice}})
#' @param runs Number of runs.
#' @param seed Seed for the kmeans function.
#' @param save If TRUE then the growth curves plot truncated at the "TruncTime" is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' @param Cores Number of cores to parallelize computations.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return StabilityAnalysis returns a list of (i) lists, called ConsensusInfo, reporting for each G and h: the Consensus Matrix, either as a NxN matrix, where N is the number of samples, or plot, and the most probable clustering obtained from running several times the method; (ii) the box plots showing both the Elbow plot considering the total tightness and the box plots of the fDB indexes for each G; and finally, (iii) the seed. See \code{\link{BoxPlot.Extrapolation}} and \code{\link{MostProbableClustering.Extrapolation}}.
#' 
#' @details 
#'  Connector provides two different plots to properly guide the choice of the number of clusters:
#'  \itemize{
#'  \item{Elbow Plot:}{ a plot in which the total tightness against the number of clusters is plotted. A proper number of clusters can be inferred as large enough to let the total tightness drop down to relatively little values but as the smallest over which the total tightness does not decrease substantially (we look for the location of an "elbow" in the plot). }
#'  \item{Box Plot:}{ a plot in which for each number of clusters , G, the functional Davies-Buldin (fDB), the cluster separation measure index, is plotted as a boxplot. A proper number of clusters can be associated to the minimum fDB value. }
#'  }
#'  
#'    The proximities measures choosen is defined as follow
#'  \deqn{D_q(f,g) = \sqrt( \integral | f^{(q)}(s)-g^{(q)}(s) |^2 ds ), d=0,1,2}
#' where f and g are two curves and f^{(q)} and g^{(q)} are their q-th derivatives. Note that for q=0, the equation becomes the distance induced by the classical L^2-norm.
#' Hence, we can define the following indexes for obtaining a cluster separation measure:
#' \itemize{
#' \item{T:}{ the total tightness representing the dispersion measure given by
#' \deqn{ T = \sum_{k=1}^G \sum_{i=1}^n D_0(\hat{g}_i, \bar{g}^k)}; }
#'  \item{S_k:}{ 
#'  \deqn{ S_k = \sqrt{\frac{1}{G_k} \sum_{i=1}^{G_k} D_q^2(\hat{g}_i, \bar{g}^k);} }
#'  with G_k the number of curves in the k-th cluster;
#'  }
#'  \item{M_{hk}:}{ the distance between centroids (mean-curves) of h-th and k-th cluster 
#'  \deqn{M_{hk} =  D_q(\bar{g}^h, \bar{g}^k);}
#'  }
#'  \item{R_{hk}:}{  a measure of how good the clustering is,
#'  \deqn{R_{hk} =   \frac{S_h + S_k}{M_{hk}};}
#'  }
#'  \item{fDB_q:}{ functional Davies-Bouldin index, the cluster separation measure
#'  \deqn{fDB_q = \frac{1}{G} \sum_{k=1}^G \max_{h \neq k} {  \frac{S_h + S_k}{M_{hk}} } }
#'  }
#' }
#' 
#'
#' @seealso MostProbableClustering.Extrapolation, BoxPlot.Extrapolation, ConsMatrix.Extrapolation.
#' 
#' @import ggplot2 reshape2 RColorBrewer statmod parallel
#' @importFrom cowplot plot_grid
#' @export

StabilityAnalysis<-function(data,G,h,p,runs=50,seed=2404,save=FALSE,path=NULL,Cores=1)
{
  ALL.runs<-list()
  ConsensusInfo<-list()
  
  set.seed(seed)
  seed<-.Random.seed 
  assign(x = ".Random.seed", value = seed, envir = .GlobalEnv)
  
  nworkers <- detectCores()
  if(nworkers<Cores) Cores <- nworkers
  
  cl <- makeCluster(Cores)

  CONNECTORList <- data
  G = G
  h = h
  p = p
  
  ALL.runs<-parLapply(cl,1:runs, function(i) ClusterChoice(data = CONNECTORList, G = G, h = h, p = p,seed=NULL) )

  stopCluster(cl)
###### box plot generation #####
  Box.pl<-list()
  
  for(hind in 1:length(h) )
  {
    pl<-list()
    
    ##### Calculation of the tightness and DB (0,1,2) indexes
    l.tight<-lapply(1:runs,function(x) ALL.runs[[x]]$Tight.indexes[,hind])
    l.fDB<-lapply(1:runs,function(x) ALL.runs[[x]]$fDB.indexes[,hind])
    l.fDB1<-lapply(1:runs,function(x) ALL.runs[[x]]$fDB1deriv.indexes[,hind])
    l.fDB2<-lapply(1:runs,function(x) ALL.runs[[x]]$fDB2deriv.indexes[,hind])
    #####
    
    dt.fr<-do.call(cbind, l.tight)
    
    row.names(dt.fr)=paste("G=",G) 
    
    dt.fr<-data.frame(clust=rep(paste("G=",G),length(dt.fr[1,])),y=round(c(dt.fr),digits = 3) )
    counts<-count(dt.fr,vars=c("clust","y"))
    
    pl[["Tight"]]<-ggplot()+geom_boxplot(data= dt.fr,aes(x=clust,y=y))+geom_point(data=counts, col="red",aes(x=clust,y=y,size=freq/runs) )+
      labs(title=paste("Elbow plot ( h =",h[hind],")" ),x="Number of Clusters",y="Tightness (T)")+
      theme(text = element_text(size=20)) +labs(size="Counts freq.") 
    
    
    Tight<-do.call(cbind, l.tight)
    fDB<-do.call(cbind, l.fDB)
    fDB.der1<-do.call(cbind, l.fDB1)
    fDB.der2<-do.call(cbind, l.fDB2)
    
    row.names(fDB)=paste("G=",G) 
    
    dt.fr2<-data.frame( clust=rep(paste("G=",G) ,length(fDB[1,])) ,y=round(c(fDB),digits = 3) )
    counts2<-count(dt.fr2,vars=c("clust","y"))
    
    pl[["fDBindex"]]<-ggplot()+geom_boxplot(data= dt.fr2,aes(x=clust,y=y))+geom_point(data=counts2, col="red",aes(x=clust,y=y,size=freq/runs) )+
      labs(title=paste("fDB indexes ( h =",h[hind],")" ),x="Number of Clusters",y="fDB index",col="")+
      theme(text = element_text(size=20))+labs(size="Counts freq.") 
    
    boxplots<-plot_grid(plotlist=pl)
    
    Box.pl[[paste("h=",h[hind])]]$Plot<-list(Both=boxplots,fDBindex = pl[["fDBindex"]], Elbow= pl[["Tight"]])
    Box.pl[[paste("h=",h[hind])]]$Data<-list(Tight=Tight,fDB=fDB,fDB.1=fDB.der1,fDB.2=fDB.der2)
    
    ######### Cluster Membership #########
    
    ConsM.generation<-function(Gind,runs,ALL.runs,hind=hind) {
      
      ############# first we found the most probably clustering:
      ClustCounting<-sapply(1:runs,function(x) names(ALL.runs[[x]]$FCM_all[[paste("G=",Gind)]][[paste("h=",h[hind])]]$FCM$cluster$cluster.member ) )
      
      ClustString<-sapply(1:runs,function(x) paste ( table(ClustCounting[,x]) , collapse = "") )
      IndexBestClustering<-which(ClustString==names(which.max(table(ClustString))))[1]
      
      BestClustering<-ALL.runs[[IndexBestClustering]]$FCM_all[[paste("G=",Gind)]][[paste("h=",h[hind])]]
      
      ##########################################################
      #### Let build the consensus matrix
      
      consensus.runs<-lapply(1:runs, function(x){
        consensusM<-diag(0,length(data$LenCurv))
        fcm.G<-ALL.runs[[x]]$FCM_all[[paste("G=",Gind)]][[paste("h=",h[hind])]]
        cl.vector<-fcm.G$FCM$cluster$cluster.member
        for(i in 1:length(data$LenCurv)) consensusM[i,which(cl.vector%in%cl.vector[i])]<-1
        consensusM
      })
      
      consensusM<-Reduce("+",consensus.runs)
      
      colnames(consensusM)<- data$LabCurv$ID
      row.names(consensusM)<- data$LabCurv$ID
      
      consensusM<-as.data.frame(consensusM)/runs
      
      ################
      ### ordering the samples to have all curves in the same cluster close with each otehrs.
      
      cl.memer<-BestClustering$FCM$cluster$cluster.m
      #  mat[do.call(order, as.data.frame(mat)),]
      
      fcm<-BestClustering
      curvename<-data$LabCurv$ID
      curve<-fcm$FCM$prediction$gpred
      rownames(curve)<-curvename
      
      grid<-fcm$FCM$fit$grid
      
      gauss.quad(10) -> gauss
      
      a <- min(grid)
      b <- max(grid)
      itempi <- (a+b)/2 + (b-a)/2*gauss$nodes
      
      match(itempi,grid) -> itimeindex 
      
      fxG <- (curve[,itimeindex] )^2
      int <- (b-a)/2 * rowSums( gauss$weights * fxG )
      dist.curve <- sqrt(int)
      names(which.min(dist.curve))->lowest.curve
      
      m.lowercurve<-matrix(curve[lowest.curve,itimeindex],nrow = length(curve[,1]),ncol = length(itimeindex),byrow = T)
      fxG <- (curve[,itimeindex]-m.lowercurve )^2
      int <- (b-a)/2 * rowSums( gauss$weights * fxG )
      dist.curve <- sqrt(int)
      ################## Sorting the names!!!
      
      #1) sorting depending by the l2 distance with the zero x-axis!
      curvename.ordered<-names(sort(dist.curve)) 
      
      #2) Sorting depending on the cluster membership
      names(cl.memer)<-curvename
      cl.memer<-sort(cl.memer)
      #3) Defining a dataframe in order to sort the curves depending by the cluster and distance!
      
      namefroml2<-1:length(curvename)
      names(namefroml2)<-curvename.ordered
      namefroml2<-namefroml2[names(cl.memer)] # ordering the curve ordered by the L2 distance
      
      curve.name.dtfr<-data.frame(namefroml2=namefroml2,cl.mem=cl.memer,1:length(curvename.ordered))
      
      # ind<-curve.name.dtfr[order(curve.name.dtfr$namefroml2,curve.name.dtfr$cl.mem,curve.name.dtfr$namefroml2),3]
      ind<-curve.name.dtfr[order(curve.name.dtfr$cl.mem,curve.name.dtfr$namefroml2),3]
      
      curvename.ordered<-names(cl.memer)[ind]
      
      consensusM <- consensusM[curvename.ordered,curvename.ordered]
      consensusM <- as.matrix(consensusM)
      consensusM[upper.tri(consensusM)] <- NA
      
      m<-melt( as.matrix(consensusM) , na.rm = TRUE)
      
      m$Var2<-factor(m$Var2,levels = curvename.ordered  )
      m$Var1<-factor(m$Var1,levels = rev(curvename.ordered )  )
      
      ############
      length.cl<-rev(as.numeric(table(cl.memer[ind])[unique(cl.memer[ind])]))
      coor<-c(0,cumsum(length.cl))
      coor.2<-sum(length.cl)-coor
      
      cluster.lines<-data.frame(x=c(coor[-length(coor)],coor[-length(coor)]),
                                y=c(coor.2[-1],coor.2[-length(coor.2)]),
                                xend=c(coor[-1],coor[-length(coor)]),
                                yend=c(coor.2[-1],coor.2[-1] ) ) +.5
      
      x.text<-(cluster.lines$xend-cluster.lines$x)/2 + cluster.lines$x
      y.text<-(cluster.lines$yend-cluster.lines$y)/2 + cluster.lines$y
      x.text<-x.text[1:Gind]
      y.text<-y.text[-(1:Gind)]
      
      ConsensusPlot<- ggplot() + 
        geom_tile(data=m,aes(x = Var1, y = Var2,fill=value))+
        labs(x="", y="",fill="Same cluster \ncounting") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_fill_gradient2(midpoint=0.5, low="blue", mid="yellow",
                              high="red")+
        geom_segment(data=cluster.lines, aes(x,y,xend=xend, yend=yend), size=1.5, inherit.aes=F)+
        labs(title="Consensus Matrix",subtitle = "Black line for the most probable clustering" )+
        annotate(geom="text",x=x.text, y=y.text, label= rev(BestClustering$FCM$cluster$cluster.names),size=6)
      
      return(list(ConsensusMatrix=consensusM,ConsensusPlot=ConsensusPlot,MostProbableClustering=BestClustering) )
    } 
    
    
    Consensus.Info<-lapply(G,runs=runs,ALL.runs=ALL.runs,hind=hind, ConsM.generation)

    names(Consensus.Info) <- c(paste("G=",G))
    
    ConsensusInfo[[paste("h=",h[hind])]]<-Consensus.Info
    
    #####################################
    if(save==TRUE)
    {
      if(is.null(path))
      {
        path <- getwd()
      }
      
      ggsave(filename=paste0("BoxPlotClusterGoodness_p",p,"_h",h[hind],"_runs",runs,"_G_",min(G),"_",max(G),".pdf"),plot =Box.pl[[hind]]$Plot,width=50, height = 20, units = "cm",scale = 1,path=path )

      # for(Gind in G)
      # {
      #   ggsave(filename=paste0("ConsensusMatrix_H",h[hind],"G",Gind,".pdf",sep=""),plot = Consensus.Info$ConsensusPlot[[paste("G=",Gind)]]$ConsensusPlot,width=29, height = 20, units = "cm",scale = 1,path=path )
      # }
      
    }

  }
  #Box.pl$seed <- seed
  return( list(ConsensusInfo=ConsensusInfo,BoxPlots=Box.pl,seed=seed) )
}





