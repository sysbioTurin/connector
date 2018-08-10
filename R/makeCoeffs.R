makeCoeffs <- function(data, base=NULL, reg, dimBase, grid=NULL, pert, baseType=baseType){
  if(is.null(base)){
    time<-grid
    nbasis <- dimBase
    m <- length(time)
    bObj <-  create.bspline.irregular(c(time[1],time[m]),
                                      nbasis=nbasis,
                                      norder=min(nbasis, 4))
    tempBase<- eval.basis(time, bObj)
    base <- svd(tempBase)$u
  }else{
    base <- base[,1:dimBase]
  }
  if(reg){
    coeffs <- t((solve(t(base) %*% base + pert *diag(dimBase))%*%t(base))%*% data)
    fullBase <- base
  }else{
    curveIndx <- data[,1]
    timeIndx <- match(data[,3],grid)
    n <- max(curveIndx)
    fullBase <- base[timeIndx,  ]
    coeffs <- matrix(0,nrow=n,ncol=sum(dimBase))
    for (i in 1:n){
      if(is.null(dim(base)[1]))
        base <- t(t(base))
      basei <- fullBase[curveIndx==i,]
      yi <- data[curveIndx==i,2]
      if(length(yi)>1){
        coeffs[i,] <- solve(t(basei) %*% basei + pert * diag(dimBase)) %*% t(basei) %*%yi
      }else{
        coeffs[i,] <- ((basei) * basei + pert )^(-1) * (basei)*yi
      }
    }
  }
  return(list(coeffs=coeffs, base=base, fullBase=fullBase, dimBase=dimBase))
}
