#' Shuffle Analysis
#'
#' @import  statmod
#' 
#'
ShuffleAnalysis<-function(clusterdata,nruns=10^3)
{
  
  grid <- clusterdata$FCM$TimeGrid
  
  
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
  
  #################################
  
  ncurve<-max(unique(clusterdata$FCM$fit$data$curve))
  
  real.cluster<- clusterdata$FCM$cluster$cluster.member
  fcm.prediction<-clusterdata$FCM$prediction
  fcm.fit<-clusterdata$FCM$fit
 
  DBshuffle<-list()
  DBshuffle[["First"]][["DB.index"]]<-40
  DBshuffle[["Second"]][["DB.index"]]<-40
  
  for(i in 1:nruns){
      fake.indexes<-sample.int(ncurve, ncurve, replace = F)
      fake.cluster<-real.cluster[fake.indexes]
    
      Coefficents<-Rclust.coeff(clust=fake.cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 0)
      
      Deriv.Coefficents<-Rclust.coeff(clust=fake.cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 1)
      
      DB.temp <- Coefficents
      
      if( DB.temp$DB.index < DBshuffle[["Second"]][["DB.index"]]) {
        if( DB.temp$DB.index < DBshuffle[["First"]][["DB.index"]]  ){
          DBshuffle[["First"]] <- DB.temp
          DBshuffle[["First"]][["Cluster"]]<-fake.cluster
        } 
        else{
          DBshuffle[["Second"]] <- DB.temp
          DBshuffle[["Second"]][["Cluster"]]<-fake.cluster
        } 
      }
  }
  
return(DBshuffle)

}
