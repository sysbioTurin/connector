
####################

WithCluster_MeanDist <- function(ClustCurve,MeanCurves,i)
{
  #' Computes the Hausdorff distance between the curves belonging to the i-th cluster and the corresponding mean curve.
  
  TimeGrid <- sort(unique(ClustCurve[,2]))
  ### i-th cluster curves data
  ClustCurve.i <- ClustCurve[ClustCurve[,4]==i,]
  ### i-th cluster curves ID
  index <- sort(unique(ClustCurve.i[,1]))
  ### i-th cluster magnitudo
  ni <- length(index)
  ### i-th cluster max obs time
  tmax <- max(ClustCurve.i[,2])
  tmin <- min(ClustCurve.i[,2])
  ### i-th cluster meancurve truncated at tmax
  MeanCurve.i <- MeanCurves[,i][which(TimeGrid <= tmax & TimeGrid >= tmin )]
  ### i-th cluster obs time grid
  TimeGrid.i <- TimeGrid[TimeGrid <= tmax & TimeGrid >= tmin]
  ### i-th cluster meancurve
  A <- matrix(c(TimeGrid.i,MeanCurve.i),ncol=2)
  
  ff<-function(k){
    ### i-th cluster curves
    B <- ClustCurve.i[ClustCurve.i[,1]==index[k],2:3]
    ### Hausdorff distance between i-th cluster meancurve and k-th curve respectively
    hausdorff(A,B)
  }
  withmean.i<-sapply(1:ni, ff)
  #cat("cluster ",i, " ID:",index," distanze:\n ",withmean.i,"\n sum:",sum(withmean.i),"\n sum(^2):",sum(withmean.i^2),"\n sum^2:",sum(withmean.i)^2,"\n mean:",mean(withmean.i), "\n\n" )
  return(withmean.i)
  ### Returns a vector reporting the hausdorff distance of the i-th cluster curves from the corresponding cluster meancurve.
}



##################

DataFrameWithinness.2 <- function(clusterdata,Meancurves,cluster,shift,ClustSymbol,K)
{
  if(isS4(clusterdata))
  {
    dist<-dist2centers.new(clusterdata@data, Meancurves)
  }
  else {
    dist<-dist2centers.new(clusterdata$Summary[,c(1,3,2)], Meancurves)
  }
  
  min.dist<-makeClMat(dist)[,1]
  Withinness.i <-lapply( 1:K, function(x){
    with<-min.dist[which(cluster==x)]
    mean.dist <- mean(with)
    st.dev <- sd(with)
    sh<-shift[x]
    
    ### Data frame for ggplot
    circles.i <- data.frame(
      x0 <- rep(0+sh,3),
      y0 <- rep(0,3),
      r <- rep(mean.dist, 3)+c(st.dev,0,-st.dev),
      distance <- c("sd","mean","sd"),
      Cluster <- factor(ClustSymbol[x],levels=ClustSymbol)
    )
    colnames(circles.i) <- c("x0","y0","r","distance","Cluster")
    
    WithDist.i <- data.frame(
      x1 <- with + sh,
      y1 <- numeric(length(with)),
      Cluster <- factor((x-1),levels=c(0:(K-1)))
    )
    
    colnames(WithDist.i) <- c("x1","y1","Cluster")
    
    return(DataFrame.i=list(circles=circles.i,WithDist=WithDist.i))
    
  })
  
  return(Withinness.i)
}


With_coeff<-function(clusterdata,data)
{
  if(isS4(clusterdata))
  { 
    ### Hausdorff withiness
    out.fit<-clusterdata@models$fitfclust@fit
    within<-c()
    database<-data$Dataset
    FCM.ClustCurve<-data.frame(ID=database[,1],Times=database[,3],Vol=database[,2],Cluster=rep(clusterdata@allClusters,times=data$LenCurv))
    
    if(clusterdata@reg==1)
    {
      fitfclust.curvepred(out.fit)$meancurves->meancurves
    } else{
      fitfclust.curvepredIrreg(out.fit)$meancurves->meancurves
    }
    clusterdata@k->k
    ### i-th cluster curves and meancurve distance
    within <-sapply(1:k, function(i){ sum( WithCluster_MeanDist(FCM.ClustCurve,meancurves,i) ) } )
    
    HausWith<-sum(within)
    
    ### Euclidian withiness
    dist<-dist2centers.new(clusterdata@data, meancurves)
    minDist<-makeClMat(dist)[,1]
    EucWith<-sum(minDist)
    
  }else{
    EucWith<-clusterdata$Tot.within
    ### i-th cluster curves and meancurve distance
    MeanCurves <- clusterdata$meancurves
    ClustCurve <- clusterdata$Summary
    k<-length(MeanCurves[1,])
    within <-sapply(1:k, function(i){ sum( WithCluster_MeanDist(ClustCurve,MeanCurves,i) ) } )
    
    HausWith<-sum(within)
  }
  
  coeff<-(HausWith/(HausWith+EucWith)-.5)*2
  return(list(EucTight=EucWith,HausTight=HausWith,coeff=coeff))
}

dist2centers.new <- function(data, centers){
  chf <- checkFormat(data)
  reg <-chf$reg
  data <- chf$data
  if(reg){
    nc <- dim(data)[2]
    nt <- dim(data)[1]
    res <- t(sapply(1:nc, function(x) colSums(sqrt((data[,x]-centers)^2)))/nt)
  }else{
    dataNew <- formatFuncy(data, format="Format3")
    res <- with(dataNew, makeSparse(Yin, Tin, isobs,
                                    t_all))
    out_indx <- match(sort(unique(data[,3])), dataNew$t_all)
    data <- as.matrix(res$longYin)
    isobs <- as.matrix(res$longIsobs)
    if(dim(centers)[2]==1)
      sumFct <- sum
    else
      sumFct <- colSums
    
    nc <- dim(data)[1]
    
    if(!is.null(time))
      centers <- centers[out_indx,]
    res <- t(sapply(1:nc, function(x){
      dist <- sqrt((data[x,]-centers)[isobs[x,]==1,]^2)
      if(is.null(dim(dist)) | is.null(dim(centers)))
        return(dist/sum(isobs[x,]))
      else
        return(colSums(dist)/sum(isobs[x,]))
    }
    )
    )
    
  } 
  return(res)
}