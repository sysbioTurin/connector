#' PCA Analysis
#'
#'@description
#'
#'Generates a bar plot reporting the percentage values associated with each component identified by the Principal Component Analysis (PCA). Considering the Functional Clustering Model [Sugar and James] the PCA is applied on the natural cubic spline coefficients to estimate the dimension of the subspace, h, in which the cluster means lie. This is due to a parameterization of the cluster means for producing low-dimensional representations of the curves, for this reason \eqn{h\le p}.
#'
#' @param data CONNECTORList.  (see \code{\link{DataImport}})
#' @param p The dimension of the natural cubic spline basis.
#' @param save  If TRUE the bar plot of the PCA components is saved into a pdf file.
#' @param path The folder path   where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @return  PCA Analysis returns the bar plot of the PCA components and the vector of percentage values associated with the PCA components.
#' @examples
#'
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#' CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))
#'
#' PCA<-PCA.Analysis(CONNECTORList)
#' PCA$plot
#' 
#' @references
#' Gareth M. James and Catherine A. Sugar, (2003). Clustering for Sparsely Sampled Functional Data. Journal of the American Statistical Association.
#' 
#' @seealso \code{\link[fda]{create.bspline.irregular}}.
#' 
#' @import ggplot2 fda
#' @export
PCA.Analysis <- function(data,p=5,save=FALSE,path=NULL)
{
  database<-data$Dataset
  
  data <-matrix(c(database$ID,database$Vol,database$Time),ncol=3,byrow=F)
  
  
  chf <- checkFormat(data)
  reg <- chf$reg
  data<-chf$data
  if(reg){
    TimeGrid<-1:dim(data)[1]
  }else{
    TimeGrid <- unique(sort(data[,3]))
  }
  
  # curves splines basis coefficients
  m <- length(TimeGrid)
  bObj <-  create.bspline.irregular(c(TimeGrid[1],TimeGrid[m]),
                                    nbasis=p,
                                    norder=min(p, 4))
  
  tempBase <- eval.basis(TimeGrid, bObj)
  
  base <- svd(tempBase)$u
  pert<-0.01
  if(reg){
    coeffs <- t((solve(t(base) %*% base + pert *
                         diag(p))%*%t(base)) %*% data)
  }else{
    curveIndx <- data[,1]
    timeIndx <- match(data[,3],TimeGrid)
    n <- max(curveIndx)
    fullBase <- base[timeIndx,  ]
    coeffs <- matrix(0,nrow=n,ncol=sum(p))
    for (i in 1:n){
      if(is.null(dim(base)[1]))
        base <- t(t(base))
      basei <- fullBase[curveIndx==i,]
      yi <- data[curveIndx==i,2] 
      if(length(yi)>1){
        coeffs[i,] <- solve(t(basei) %*% basei + pert * diag(p)) %*% t(basei) %*%yi
      }else{
        coeffs[i,] <- ((basei) * basei + pert )^(-1) * (basei)*yi
      }
    }
  }
  
  # Principal Components Analysis

  princomp(as.matrix(coeffs)) -> pca
  # Number of principal components
  ncomp <- length(names(pca$sdev))
  # Principal components variances
  eigs <- pca$sdev^2
  # Percentage of variances explained by each component
  percentage <- eigs/sum(eigs)*100

  # PCA bar plot
  PCA_barplot<-ggplot(data=data.frame(comp=paste("Comp.",1:ncomp),Variances=eigs,perc=paste(signif(percentage,4),"%",sep="")), aes(x=comp, y=Variances)) +
    geom_bar(stat="identity", fill="steelblue")+
    geom_text(aes(label=perc), vjust=-.3,  size=3.5)+
    labs(title="PCA barplot", x="Components", y = "Variances")+
    theme(plot.title = element_text(hjust = 0.5))

  if(save==TRUE)
  {
    ggsave(filename="PCA_Analysis.pdf",plot =PCA_barplot,width=29, height = 20, units = "cm",scale = 1,path = path)
  }
  return(list(plot=PCA_barplot,perc=percentage))
  }

