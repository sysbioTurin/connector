###########################################################################
###################### Davies–Bouldin index (DBI) #########################
###########################################################################
#' @import  splines statmod
#' @export
#' 
###############
# L2dist.curve2mu calculates the L2 (or with derivatives) distance between all the spline estimated for each sample with
# respect to the corresponding cluster center curve.
##############

L2dist.curve2mu <- function(clust,fcm.curve,gauss.info,fcm.fit=NULL,deriv=0){

  itimeindex <- gauss.info$itimeindex
  weights <- gauss.info$gauss$weights 
  a<-gauss.info$a
  b<-gauss.info$b
  
  n.curves<-length(fcm.curve$gpred[,1])
  
  if(deriv==0)
  {
    fxk <- (fcm.curve$gpred[1:n.curves,itimeindex]-t(fcm.curve$meancurves[itimeindex,clust[1:n.curves]]) )^2
  
  }else{
    
    if(is.null(fcm.fit)) warning("The fcm.fit is needed to calculate the derivatives!! ")
    
      dspl <-basis.derivation(fcm.fit,deriv)
      u.dspl<-dspl$u.dspl
      dmeancurves<-dspl$dmeancurves
      
      etapred <- fcm.curve$etapred
      
      matrix(0,n.curves,nrow(u.dspl)) -> dgpred
      
      for (ind in 1:n.curves){
        dgpred[ind,] <- as.vector(u.dspl %*% etapred[ind,])
      }

      fxk <- ( dgpred[1:n.curves,itimeindex]- t(dmeancurves[itimeindex,clust[1:n.curves]]) )^2

  }
  
  int <- (b-a)/2 * rowSums( weights * fxk )
  dist.curve2mu <- sqrt(int)
  
  return(dist.curve2mu)
}

###############
# L2dist.mu2mu calculates the L2 distance between the cluster centers curves.
##############

L2dist.mu2mu <- function(fcm.curve,gauss.info,fcm.fit=NULL,deriv=0){
  
  itimeindex <- gauss.info$itimeindex
  weights <- gauss.info$gauss$weights 
  a<-gauss.info$a
  b<-gauss.info$b
  
  fcm.curve$meancurves->meancurves
  
  k<-length(meancurves[1,])
  dist.mu2mu<-matrix(0,ncol = k,nrow = k)
  
 #clust.ind <- combn(1:k,2)
  clust.ind <- expand.grid(1:k,1:k)
  
  if(deriv==0)
  {
    fxk <- (meancurves[itimeindex,clust.ind$Var1]-meancurves[itimeindex,clust.ind$Var2])^2
  }else{
    dspl <-basis.derivation(fcm.fit,deriv)
    u.dspl<-dspl$u.dspl
    dmeancurves<-dspl$dmeancurves
    
    fxk <- (dmeancurves[itimeindex,clust.ind$Var1]-dmeancurves[itimeindex,clust.ind$Var2])^2

  }
  
  if(length(as.matrix(fxk)[1,] ) == 1 ){
    int <- (b-a)/2 * sum(weights * fxk)
  }else int <- (b-a)/2 * colSums(weights * fxk)
  
  dist.mu2mu[,] <- sqrt(int)

  return(dist.mu2mu)
}

###############
# L2dist.mu20 calculates the L2 distance between the cluster centers curves and zero.
##############

L2dist.mu20 <- function(fcm.curve,gauss.info,fcm.fit=NULL,deriv=0){
  
  itimeindex <- gauss.info$itimeindex
  weights <- gauss.info$gauss$weights 
  a<-gauss.info$a
  b<-gauss.info$b
  
  fcm.curve$meancurves->meancurves
  
  k<-length(meancurves[1,])
  
  if(deriv==0)
  {
    fxk <- (meancurves[itimeindex,])^2
  }else{
    dspl <-basis.derivation(fcm.fit,deriv)
    u.dspl<-dspl$u.dspl
    dmeancurves<-dspl$dmeancurves
    
    fxk <- (dmeancurves[itimeindex,])^2
    
  }
  
  if(length(as.matrix(fxk)[1,] ) == 1 ){
    int <- (b-a)/2 * sum(weights * fxk)
  }else int <- (b-a)/2 * colSums(weights * fxk)
  
  dist.mu20 <- sqrt(int)
  
  return(dist.mu20)
}
###############
# Sclust.coeff calculates the S coefficents for each cluster.
##############

Sclust.coeff <- function(clust,fcm.curve,gauss.info,fcm.fit=NULL,deriv=0){

  distances <- L2dist.curve2mu(clust,fcm.curve,gauss.info,fcm.fit,deriv)
  k<-length(fcm.curve$meancurves[1,])

  out <- sapply(1:k,
   function(x){
    curve <- which(clust == x)
    if(length(curve)==0) warning("Empty clusters imply a NaN coefficents",immediate. = TRUE)
    sqrt(1/length(curve)*sum(distances[curve]^2))
   }
 )
  
  return(out)
}

###############
# Rclust.coeff calculates the R coefficents for each cluster.
###############

Rclust.coeff <- function(clust,fcm.curve,gauss.info,fcm.fit=NULL,deriv=0){
  k<-length(fcm.curve$meancurves[1,])
  
  emme <- L2dist.mu2mu(fcm.curve,gauss.info,fcm.fit,deriv)
  cl   <- L2dist.mu20(fcm.curve,gauss.info,fcm.fit,deriv)
  ######### Let name the cluster with A->Z from the lower mean curve to the higher.
 
  M <- emme
  
  
  if( k !=1 )
  {
    Cl.order<-order(cl)
  }else{
    Cl.order<-1
  }
  
  symbols<-LETTERS[order(Cl.order)]
  
  
  row.names(emme)<-symbols
  colnames(emme)<-symbols
  
  esse <- Sclust.coeff(clust,fcm.curve,gauss.info,fcm.fit,deriv) 
  names(esse)<-symbols
  
  erre <- matrix(0,nrow=k,ncol=k,dimnames = list(symbols,symbols))

  S.matrix<-matrix(rep(esse,length(esse)),ncol=length(esse))
  sum.2S <- S.matrix + t(S.matrix) 
  erre <- sum.2S / emme
  diag(erre)<-0
  errei <- apply(erre,1,max)
  errebar <- mean(errei)
  
  return(list(erre=erre,errei=errei,DB.index=errebar,emme=emme,esse=esse))
}

#######################################################################################à


###############
# basis.derivation calculates the derivatives of the spline basis considering the matrix decomposition.
###############

basis.derivation<- function(fcm.fit,deriv)
{
  Lambda <- fcm.fit$par$Lambda
  alpha <- fcm.fit$par$alpha
  lambda.zero <- as.vector(fcm.fit$par$lambda.zero)
  Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
  q<-length(fcm.fit$FullS[1,])
  
  grid <- fcm.fit$grid
  
  cbind(1,ns(grid, df = (q - 1))) -> spl
  singular.spl <- svd(spl)
  cbind(0,ns.deriv(grid, df = (q-1), deriv = deriv)) -> dspl
  utrasf <- function(dspl){
    out <- matrix(0, nrow = dim(dspl)[1], ncol = dim(dspl)[2])
    for (colonna in 1:dim(dspl)[2]){
      out[,colonna] <- 1/singular.spl$d[colonna] * dspl %*% singular.spl$v[,colonna]
    }
    return(out)
  }
  u.dspl <- utrasf(dspl)
  dmeancurves <- u.dspl %*% Lambda.alpha
  return(list(u.dspl=u.dspl,dmeancurves=dmeancurves) )
}

#################
## ns.deriv is the ns function modified in order to make the function splineDesign calculates the derivatives!
#################

ns.deriv <- function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), deriv) {
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  }
  else {
    if (length(x) == 1L) 
      Boundary.knots <- x * c(7, 9)/8
    FALSE
  }
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d", 
                       1L + intercept), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                         2L)[-c(1L, nIknots + 2L)]
      quantile(x[!outside], knots)
    }
  }
  else nIknots <- length(knots)
  
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    warning(gettextf("I cannot derive"), domain = NA)
    #     basis <- array(0, c(length(x), nIknots + 4L))
    #     if (any(ol)) {
    #         k.pivot <- Boundary.knots[1L]
    #         xl <- cbind(1, x[ol] - k.pivot)
    #         tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
    #             1))
    #         basis[ol, ] <- xl %*% tt
    #     }
    #     if (any(or)) {
    #         k.pivot <- Boundary.knots[2L]
    #         xr <- cbind(1, x[or] - k.pivot)
    #         tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
    #             1))
    #         basis[or, ] <- xr %*% tt
    #     }
    #     if (any(inside <- !outside)) 
    #         basis[inside, ] <- splineDesign(Aknots, x[inside], 
    #             4)
  }
  else basis <- splineDesign(Aknots, x, ord = 4L, derivs = deriv)
  const <- splineDesign(Aknots, Boundary.knots, ord = 4L, derivs = c(2L, 
                                                                     2L))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
                                                     drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots, 
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}
