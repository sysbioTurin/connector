#' FCM functions
#'
#' @import  fda splines 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#' 

"fitfclust" <-
  function(x=NULL,curve=NULL,timeindex=NULL,data=NULL, q = 5, h = 1, p = 5,
           K = 2, tol = 0.001, maxit = 20, pert = 0.01, grid = seq(0, 1, length = 100), hard = F, plot= F,trace=F,seed=NULL)
  {
    # This is the main function to implement the FClust procedure.
    if (is.null(data))  data <- list(x=x,curve=curve,timeindex=timeindex)
    initfit <- fclustinit(data = data, pert = pert, grid = grid, h = h,
                          p = p, q = q, K = K,seed)
    # Initialize the parameters
    parameters <- initfit$parameters
    vars <- initfit$vars
    S <- initfit$S
    FullS <- initfit$FullS
    sigma.old <- 0
    sigma.new <- 1
    ind <- 1
    # Main loop. Iterates between M and E steps and stops when
    # sigma  has converged.
    while(abs(sigma.old - sigma.new)/sigma.new > tol & (ind <= maxit)) {
      parameters <- fclustMstep(parameters, data, vars, S, tol, p, hard)
      vars <- fclustEstep(parameters, data, vars, S, hard)
      sigma.old <- sigma.new
      sigma.new <- parameters$sigma[1]
      if (trace)
        print(paste("Iteration", ind,": Sigma = ",sigma.new))
      # Plot cluster mean curves.  
      if(plot){
        cluster.mean <- matrix(0,K,dim(FullS)[1])
        for(k in 1:K)
          cluster.mean[k,] <- FullS %*% (parameters$lambda.zero + parameters$Lambda %*% parameters$alpha[k,  ])
        plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
        for(k in 1:K)
          lines(grid,cluster.mean[k,], col = 4, lwd = 2)
      }
      ind <- ind + 1
    }
    # Enforce parameter constraint (7)
    cfit <- fclustconst(data, parameters, vars, FullS)
    list(data=data,parameters = cfit$parameters, vars = cfit$vars, FullS = FullS,grid=grid)
  }

"fclustinit" <-
  function(data, pert = 0, grid = seq(0.01, 1, length = 100), h = 1, p = 2, q = 5,K = K,seed){
    S <- FullS <- NULL
    # This function initializes all the parameters.
    # Produce spline basis matrix
    FullS <- cbind(1, ns(grid, df = (q - 1)))
    FullS <- svd(FullS)$u
    
    #####################################################
    ### modifica di simone
    #bObj <-  create.bspline.irregular(c(grid[1],grid[length(grid)]),
    #nbasis=q,
    #norder=min(q, 4))
    
    #tempBase <- eval.basis(grid, bObj)
    #FullS <- svd(tempBase)$u
    
    #################################################### 
    
    S <- FullS[data$timeindex,  ]
    N <- length(table(data$curve))
    # Compute initial estimate of basis coefficients.
    points <- matrix(0,N,sum(q))
    for (i in 1:N){
      Si <- S[data$curve==i,]
      xi <- data$x[data$curve==i]
      points[i,] <- solve(t(Si) %*% Si + pert * diag(q)) %*% t(Si) %*%xi}
    # Use k-means to get initial cluster memberships from points.
    if(K > 1)
    {
      if( !is.null(seed) ) set.seed(seed)
      
      class <- kmeans(points, K, 10)$cluster
    }
    else class <- rep(1, N)
    # Initial estimates for the posterior probs of cluster membership.
    
    piigivej <- matrix(0, N, K)
    piigivej[col(piigivej) == class] <- 1
    # Calculate coefficeints for cluster means.
    classmean <- matrix(0,K,q)
    
    ########
    ## modified by Simone 
    
    for (l in 1:K)
    {
      if(sum(class==l)>1)
        classmean[l,] <- apply(points[class==l,], 2, mean)
      else
        classmean[l,] <- points[class==l,]
    }
    #classmean[k,] <- apply(points[class==k,],2,mean)
    
    # Initialize lambda.zero, Lambda and alpha as defined in paper.
    lambda.zero <- apply(classmean, 2, mean)
    Lambda <- as.matrix(svd(scale(classmean, scale = F))$v[, 1:h])
    alpha <- scale(classmean, scale = F)%*%Lambda
    # Calculate estimates for gamma and gprod.
    gamma <- t(t(points) - lambda.zero - (Lambda %*% t(alpha[class,  ])))
    gprod <- NULL
    for(i in 1:N)
      gprod <- cbind(gprod, (gamma[i,  ]) %*% t(gamma[i,  ]))
    
    gamma <- array(gamma[, rep(1:sum(q), rep(K, sum(q)))], c(N, K, sum(q)))
    gcov <- matrix(0, sum(q), N * sum(q))
    list(S = S, FullS = FullS, parameters = list(lambda.zero = lambda.zero,
                                                 Lambda = Lambda, alpha = alpha), vars = list(gamma = gamma,
                                                                                              gprod = gprod, gcov = gcov, piigivej = piigivej))
  }

"fclustMstep" <-
  function(parameters, data, vars, S, tol, p, hard)
  {
    # This function implements the M step of the EM algorithm.
    K <- dim(parameters$alpha)[1]
    alpha <- parameters$alpha
    Lambda <- parameters$Lambda
    gamma <- vars$gamma
    gcov <- vars$gcov
    curve <- data$curve
    piigivej <- vars$piigivej
    N <- dim(gamma)[1]
    K <- dim(alpha)[1]
    h <- dim(alpha)[2]
    n <- length(curve)
    q <- dim(S)[2]
    sigma.old <- 2
    sigma <- 1
    # Compute pi.
    if(hard)
      parameters$pi <- rep(1/K, K)
    else parameters$pi <- apply(vars$piigivej, 2, mean)
    # Compute rank p estimate of Gamma
    ind <- matrix(rep(c(rep(c(1, rep(0, q - 1)), N), 0), q)[1:(N*q^2)], N * q, q)
    gsvd <- svd(vars$gprod %*% ind/N)
    gsvd$d[ - (1:p)] <- 0
    parameters$Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
    # This loop iteratively calculates alpha and then Lambda and stops
    # when they have converged.
    loopnumber <- 1
    while((abs(sigma.old[1] - sigma[1])/sigma[1] > tol) & (loopnumber <10)){
      sigma.old <- sigma
      # Calculate lambda.zero.
      gamma.pi <- diag(S %*% t(apply(gamma * as.vector(piigivej),
                                     c(1, 3), sum)[curve,  ]))
      alpha.pi <- t(matrix(apply(alpha, 2, function(x, piigivej,K)
      {rep(1, K) %*% (piigivej * x)}
      , t(piigivej), K), N, h)[curve,  ])
      lambda.zero <- solve(t(S) %*% S) %*% t(S) %*% (data$x - diag(
        S %*% Lambda %*% alpha.pi) - gamma.pi)
      x <- data$x - S %*% lambda.zero
      # Calculate alpha.
      for(i in 1.:K) {
        S.Lam <- S %*% Lambda
        S.Lam.pi <- S.Lam * piigivej[curve, i]
        if(sum(piigivej[, i]) > 10^(-4))
          alpha[i,  ] <- solve(t(S.Lam.pi) %*% S.Lam) %*%
          t(S.Lam.pi) %*% (x - diag(S %*% t(gamma[curve, i,  ])))
        else print("Warning: empty cluster")
      }
      # Calculate Lambda given alpha. This is done by iterating
      # through each column of Lambda holding the other columns fixed.
      for(m in 1:h) {
        pi.alphasq <- apply(t(piigivej) * (alpha^2)[, m], 2,sum)[curve]
        pi.alpha <- apply(t(piigivej) * alpha[, m], 2, sum)[curve]
        S.Lambda <- t(S %*% Lambda)
        if(h != 1) {
          temp <- NULL
          for(i in 1:K) {
            temp <- cbind(temp, as.vector(rep(1, h - 1) %*% matrix((S.Lambda *
                                                                      alpha[i,  ])[ - m,  ], h - 1,dim(S)[1])) * alpha[i, m])
          }
          otherLambda <- (temp * piigivej[curve,  ])%*%rep(1, K)
        }
        else otherLambda <- 0
        gamma.pi.alpha <- apply(gamma * as.vector(piigivej) *
                                  rep(alpha[, m], rep(N, K)), c(1, 3), sum)[curve,  ]
        Lambda[, m] <- solve(t(S * pi.alphasq) %*% S) %*% t(S) %*%
          (x * pi.alpha - otherLambda - (S *gamma.pi.alpha) %*% rep(1, sum(q)))
      }
      # Calculate sigma 
      ind <- matrix(rep(c(rep(F, q), rep(T, N * q)),N)
                    [1:(N * N * q)], N, N * q, byrow = T)[rep(1:N, table(curve)),]
      mat1 <- matrix(rep(S, N), n, N * q)
      mat3 <- t(mat1)
      mat3[t(ind)] <- 0
      ind2 <- matrix(rep(c(rep(F, q), rep(T, N * q)),
                         N)[1:(N * N * q)], N, N * q, byrow = T)[rep(1:N, rep(q, N)),]
      mat2 <- matrix(rep(t(gcov), N), N * q, N * q,byrow = T)
      mat2[ind2] <- 0
      econtrib <- 0
      for(i2 in 1:K) {
        vect1 <- x - S %*% Lambda %*% alpha[
          i2,  ] - (S * gamma[curve, i2,  ]) %*% rep(1, q)
        econtrib <- econtrib + t(piigivej[curve,i2] * vect1) %*% vect1
      }
      sigma <- as.vector((econtrib + sum(diag(mat1 %*% mat2 %*% mat3)))/n)
      loopnumber <- loopnumber + 1
    }
    parameters$lambda.zero <- as.vector(lambda.zero)
    parameters$alpha <- alpha
    parameters$Lambda <- Lambda
    parameters$sigma <- sigma
    parameters
  }

"fclustEstep" <-
  function(parameters, data, vars, S, hard)
  {
    # This function performs the E step of the EM algorithm by
    # calculating the expected values of gamma and gamma %*% t(gamma)
    # given the current parameter estimates.
    par <- parameters
    N <- dim(vars$gamma)[1]
    K <- dim(vars$gamma)[2]
    q <- dim(vars$gamma)[3]
    Gamma <- par$Gamma
    Lambda.alpha <- par$lambda.zero + par$Lambda %*% t(par$alpha)
    vars$gprod <- vars$gcov <- NULL
    
    cost <- matrix(0, nrow=N, ncol=K)
    
    for(j in 1:N) {
      # Calculate expected value of gamma.
      Sj <- S[data$curve == j,  ]
      nj <- sum(data$curve == j)
      invvar <- diag(1/rep(par$sigma, nj))
      Cgamma <- Gamma - Gamma %*% t(Sj) %*% solve(diag(nj) + Sj %*%
                                                    Gamma %*% t(Sj) %*% invvar) %*% invvar %*% Sj %*% Gamma
      centx <- data$x[data$curve == j] - Sj %*% Lambda.alpha
      vars$gamma[j,  ,  ] <- t(Cgamma %*% t(Sj) %*% invvar %*% centx)
      # Calculate pi i given j.
      covx <- Sj %*% par$Gamma %*% t(Sj) + solve(invvar)
      d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * par$pi
      
      #vars$piigivej[j,  ] <- d/sum(d)
      
      #### modified by Simone: inspirated by the R package!
      if(sum(d)!=0)
        vars$piigivej[j,  ] <- d/sum(d)
      else
        vars$piigivej[j,  ] <- 0
      ######################################
      
      cost[j,] <-exp( - diag(t(centx) %*% solve(covx) %*% centx)/2)
      if(hard) {
        m <- order( - d)[1]
        vars$piigivej[j,  ] <- 0
        vars$piigivej[j, m] <- 1
      }
      # Calculate expected value of gamma %*% t(gamma).
      vars$gprod <- cbind(vars$gprod, t(matrix(vars$gamma[j,  ,  ],
                                               K, q)) %*% (matrix(vars$gamma[j,  ,  ], K, q) * vars$
                                                             piigivej[j,  ]) + Cgamma)
      vars$gcov <- cbind(vars$gcov, Cgamma)
    }
    vars$loglik <- sum(log(cost))
    vars
  }


"fclustconst" <-
  function(data, parameters, vars, S)
  {
    # This function enforces the constraint (7) from the paper on the
    # parameters. This means that the alphas can be interpreted as the
    # number of standard deviations that the groups are apart etc.
    par <- parameters
    A <- t(S) %*% solve(par$sigma * diag(dim(S)[1]) + S %*% par$Gamma %*%
                          t(S)) %*% S
    svdA <- svd(A)
    sqrtA <- diag(sqrt(svdA$d)) %*% t(svdA$u)
    negsqrtA <- svdA$u %*% diag(1/sqrt(svdA$d))
    finalsvd <- svd(sqrtA %*% par$Lambda)
    par$Lambda <- negsqrtA %*% finalsvd$u
    if(dim(par$Lambda)[2] > 1)
      par$alpha <- t(diag(finalsvd$d) %*% t(finalsvd$v) %*% t(par$alpha))
    else par$alpha <- t(finalsvd$d * t(finalsvd$v) %*% t(par$alpha))
    meanalpha <- apply(par$alpha, 2, mean)
    par$alpha <- t(t(par$alpha) - meanalpha)
    par$lambda.zero <- par$lambda.zero + par$Lambda %*% meanalpha
    list(parameters = par, vars = vars)
  }

"nummax" <-
  function(X)
  {
    ind <- rep(1, dim(X)[1])
    m <- X[, 1]
    if(dim(X)[2] > 1)
      for(i in 2:dim(X)[2]) {
        test <- X[, i] > m
        ind[test] <- i
        m[test] <- X[test, i]
      }
    list(ind = ind, max = m)
  }

"fclust.pred" <-
  function(fit,data=NULL,reweight=F)
  {
    # This function produces the alpha hats used to provide a low
    # dimensional pictorial respresentation of each curve. It also
    # produces a class prediction for each curve. It takes as
    # input the fit from fldafit (for predictions on the original data)
    # or the fit and a new data set (for predictions on new data).
    if (is.null(data))
      data <- fit$data
    FullS <- fit$FullS
    par <- fit$parameters
    curve <- data$curve
    time <- data$time
    N <- length(table(curve))
    h <- dim(par$alpha)[2]
    alpha.hat <- matrix(0, N, h)
    K <- dim(fit$par$alpha)[1]
    distance <- matrix(0, N, K)
    Calpha <- array(0, c(N, h, h))
    for(i in 1:N) {
      Sij <- FullS[time[curve == i],  ]
      xij <- data$x[curve == i]
      n <- length(xij)
      Sigma <- par$sigma * diag(n) + Sij %*% par$Gamma %*% t(Sij)
      # Calculate covariance for each alpha hat.
      InvCalpha <- t(par$Lambda) %*% t(Sij) %*% solve(Sigma) %*% Sij %*%
        par$Lambda
      Calpha[i,  ,  ] <- solve(InvCalpha,tol = 1e-30)
      # Calculate each alpha hat.
      alpha.hat[i,  ] <- Calpha[i,  ,  ] %*% t(par$Lambda) %*% t(
        Sij) %*% solve(Sigma) %*% (xij - Sij %*% par$lambda.zero)
      # Calculate the matrix of distances, relative to the
      # appropriate metric of each curve from each class centroid. 
      for (k in 1:K){
        y <- as.vector(alpha.hat[i,])-fit$par$alpha[k,]
        distance[i,k] <- t(y)%*%InvCalpha %*%y}}
    # Calculate final class predictions for each curve.
    class.pred <- rep(1, N)
    log.pi <- log(fit$par$pi)
    if (!reweight)
      log.pi <- rep(0,K)
    probs <- t(exp(log.pi-t(distance)/2))
    probs <- probs/apply(probs,1,sum)
    m <- probs[,1]
    if(K != 1)
      for(k in 2:K) {
        test <- (probs[, k] > m)
        class.pred[test] <- k
        m[test] <- probs[test, k]
      }
    list(Calpha = Calpha, alpha.hat = alpha.hat, class.pred = class.pred,
         distance = distance, m = m,probs=probs)
  }

"fclust.curvepred" <-
  function(fit, data=NULL, index=NULL, tau = 0.95, tau1 = 0.975)
  {
    if (is.null(data))
      data <- fit$data
    if (is.null(index))
      index <- 1:length(table(data$curve))
    tau2 <- tau/tau1
    sigma <- fit$par$sigma
    Gamma <- fit$par$Gamma
    Lambda <- fit$par$Lambda
    alpha <- fit$par$alpha
    lambda.zero <- as.vector(fit$par$lambda.zero)
    S <- fit$FullS
    N <- length(index)
    upci <-lowci <- uppi <- lowpi <- gpred <- matrix(0,N,nrow(S))
    etapred <- matrix(0,N,ncol(S))
    ind <- 1
    Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
    for (i in index){
      y <- data$x[data$curve == i]
      Si <- S[data$time[data$curve == i],  ]
      ni <- dim(Si)[1]
      invvar <- diag(1/rep(sigma, ni))
      covx <- Si %*% Gamma %*% t(Si) + solve(invvar)
      centx <- data$x[data$curve == i] - Si %*% Lambda.alpha
      d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * fit$par$pi
      pi <- d/sum(d)
      K <- length(pi)
      mu <- lambda.zero + Lambda %*% t(alpha * pi) %*% rep(1, K)
      cov <- (Gamma - Gamma %*% t(Si) %*% solve(diag(ni) + Si %*% Gamma %*%
                                                  t(Si)/sigma) %*% Si %*% Gamma/sigma)/sigma
      etapred[ind,] <- mu + cov %*% t(Si) %*% (y - Si %*% mu)
      ord <- order( - pi)
      numb <- sum(cumsum(pi[ord]) <= tau1) + 1
      v <- diag(S %*% (cov * sigma) %*% t(S))
      pse <- sqrt(v + sigma)
      se <- sqrt(v)
      lower.p <- upper.p <- lower.c <- upper.c <- matrix(0, nrow(S), numb)
      for(j in 1:numb) {
        mu <- lambda.zero + Lambda %*% alpha[ord[j],  ]
        mean <- S %*% (mu + cov %*% t(Si) %*% (y - Si %*% mu))
        upper.p[, j] <- mean + qnorm(tau2) * pse
        lower.p[, j] <- mean - qnorm(tau2) * pse
        upper.c[, j] <- mean + qnorm(tau2) * se
        lower.c[, j] <- mean - qnorm(tau2) * se
      }
      upci[ind,] <- nummax(upper.c)$max
      lowci[ind,] <-  - nummax( - lower.c)$max
      uppi[ind,] <- nummax(upper.p)$max
      lowpi[ind,] <-  - nummax( - lower.p)$max
      gpred[ind,] <- as.vector(S %*%etapred[ind,])
      ind <- ind+1
    }
    meancurves <- S%*%Lambda.alpha
    list(etapred = etapred, gpred = gpred,  upci = upci,lowci = lowci,  uppi = uppi, lowpi = lowpi,index=index,grid=fit$grid,data=data,meancurves=meancurves)
  }

"fclust.discrim" <-
  function(fit,absvalue=F){
    S <- fit$FullS
    sigma <- fit$par$sigma
    n <- nrow(S)
    Gamma <- fit$par$Gamma
    Sigma <- S%*%Gamma%*%t(S)+sigma*diag(n)
    Lambda <- fit$par$Lambda
    discrim <- solve(Sigma)%*%S%*%Lambda
    if (absvalue)
      discrim <- abs(discrim)
    n <- ncol(discrim)
    nrows <- ceiling(sqrt(n))
    par(mfrow=c(nrows,nrows))
    for (i in 1:n){
      plot(fit$grid,discrim[,i],ylab=paste("Discriminant Function ",i),xlab="Time",type='n')
      lines(fit$grid,discrim[,i],lwd=3)
      abline(0,0)}}

"fclust.plotcurves" <-
  function(object=NULL,fit=NULL,index=NULL,ci=T,pi=T,clustermean=F,timei=F){
    if (is.null(object))
      object <- fclust.curvepred(fit)
    if (is.null(index))
      index <- 1:length(table(object$data$curve))
    r <- ceiling(sqrt(length(index)))
    par(mfrow=c(r,r))
    for (i in index){
      grid <- object$grid
      upci <- object$upci[i,]
      uppi <- object$uppi[i,]
      lowci <- object$lowci[i,]
      lowpi <- object$lowpi[i,]
      gpred <- object$gpred[i,]
      meancurves <- (object$mean)
      if (clustermean)
        yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
      else
        yrange <- c(min(lowpi),max(uppi))
      plot(grid,grid,ylim=yrange,ylab="Predictions",xlab="Time",type='n',
           main=paste("Curve ",i))
      if (clustermean)
        for (k  in 1:ncol(meancurves))
          lines(grid,meancurves[,k],col=6,lty=2,lwd=2)
      if (ci){
        lines(grid,upci,col=3)
        lines(grid,lowci,col=3)}
      if (pi){
        lines(grid,uppi,col=4)
        lines(grid,lowpi,col=4)}
      lines(grid,gpred,col=2,lwd=2)
      if(timei){
        lines(object$data$time[object$data$curve==i],object$data$x[object$data$curve==i],lwd=2)
        points(object$data$time[object$data$curve==i],object$data$x[object$data$curve==i],pch=19,cex=1.5)
      }
      else{
        lines(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],lwd=2)
        points(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],pch=19,cex=1.5)
      }
    }
  }

