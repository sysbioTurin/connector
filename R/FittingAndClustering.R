#' Fitting and Clustering
#'
#'
#' @description
#' Fits and clusters the data with respect to the Malthus, Gompertz and Logistic model. The clustering method used is the \emph{k-means} one.
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param k   The number of clusters.
#' @param model Model name, i.e. Malthus, Gompertz and Logistic. See \sQuote{Details}.
#' @param fitting.method The method to be used to fit the data with respect to the Malthus, Gompertz and Logistic model. Three methods are proposed: \code{\link[optimr]{opm}}, \code{\link[GenSA]{GenSA}}, \code{\link[DEoptim]{DEoptim}}. See \sQuote{Details}.
#' @param init A vector containing the initial values for the parameters that must be optimized. In the Malthus case the length of this vector has to be two, instead in the Gompertz and Logistic cases it has to be three. It is not necessary when  \emph{DEoptim} is used.
#' @param lower Lower bound associated with each parameter.
#' @param upper Upper bound associated with each parameter.
#' @param seed  Seed for the k-means.
#' 
#' @return FittingAndClustering returns a list containing all the information related to the clustering, i.e. the parameters estimated for each sample, the parameters of the cluster centers, the mean curve values, the samples cluster membership and the seed used in the k-means. Finally, it is also returned a data frame with four arguments: ID, time, volume and cluster membership for each sample.
#' 
#' 
#' @details 
#' The model proposed are the following:
#' \itemize{
#'  \item \emph{Malthus}:  \eqn{f(t) = V_0 e^{a\cdot t}}{f(t) = V_0 exp{a t}}, 
#'  \item \emph{Gompertz}:  \eqn{ f(t) = V_0 e^{\frac{a}{b} (1-e^{-b\cdot t})} }{ f(t) = V_0 exp{ a / b ( 1-exp{-b t} ) } },
#'  \item \emph{Logistic}:  \eqn{ f(t)= \dfrac{V_0  K}{ V_0 + ( K-V_0 )e^{-b\cdot t}} }{ f(t) = (V_0  K) / ( V_0 + ( K-V_0 )exp{-b t} ) }.
#' }
#' These are parametric models characterized by three (or two in the Malthus case) unknow parameters. In order to identify the parameters that minimize the sums of squared error (SSE), three methods are proposed:
#'\itemize{
#' \item "\emph{optimr}": the function \code{\link[optimr]{opm}}  from the package \emph{optimr} (see the \emph{optimr} manual) is used. Given that the result is a dataframe having one row for each method included in this package for which results are obtained. 
#'  \item "\emph{GenSA}": the function \code{\link[GenSA]{GenSA}} searches for global minimum of a very complex non-linear objective function with a very large number of optima.
#'  \item "\emph{DEoptim}": the function \code{\link[DEoptim]{DEoptim}} performs evolutionary global optimization via  Differential Evolution algorithm.
#' }
#' 
#' 
#' @references 
#' 
#' Benzekry, Sébastien AND Lamont, Clare AND Beheshti, Afshin AND Tracz, Amanda AND Ebos, John M. L. AND Hlatky, Lynn AND Hahnfeldt, Philip. (2014) Classical Mathematical Models for Description and Prediction of Experimental Tumor Growth. PLOS Computational Biology.
#' 
#' @examples
#' ### Data files
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#' 
#' ### Merge curves and target file
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#' 
#' ############################## MALTHUS ##############################
#' lower<-c(10^(-5),0)
#' upper<-c(10^2,10^3)
#' init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=1)
#' 
#' 
#' Malthus1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",feature="Progeny",fitting.method="optimr",lower=lower,upper=upper,init=init)
#' Malthus2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Malthus",feature="Progeny",fitting.method="GenSA",lower=lower,upper=upper,init=init)
#' 
#' ############################## LOGISTIC ##############################
#'
#'lower<-c(10^(-5),0,0)
#'upper<-c(10^2,10^5,1)
#'init<- list(V0=max(0.1,min(CONNECTORList$Dataset$Vol)),a=.5, b=.5)
#'
#'Logistic1<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",feature="Progeny",fitting.method="optimr",lower=lower,upper=upper,init=init)
#'Logistic2<- FittingAndClustering(data = CONNECTORList, k = 4, model="Logistic",feature="Progeny",fitting.method="DEoptim",lower=lower,upper=upper)
#'
#' @import GenSA DEoptim optimr
#' @export
FittingAndClustering <- function(data,k,model, fitting.method="optimr",init=NULL,lower,upper,seed=NULL)
{

  data->dati
  time <- dati$Dataset$Time
  t <- sort(unique(time))
  x <- dati$Dataset$Vol
  curve <- dati$Dataset$ID
  number_curves<-max(curve)
  
  data.fit <- matrix(c(curve,x,time),ncol=3,byrow=F)
  meancurves <- matrix(numeric(length(t)*k),ncol=k)
  
  
if((length(lower)!=3 || length(upper)!=3 ))
{
    if(model=="Malthus" & (length(lower)!=2 || length(upper)!=2 ) ) 
    {
      stop("Considering the Malthus model, the length of buonds should be 2.")
    }
    else{
      stop("Considering the Logistic or Gompertz model, the length of buonds should be 3.")
    }
}
 
if(fitting.method!="DEoptim")
{
  if((length(init)!=3 ))
  {
    if(model=="Malthus" & (length(init)!=2 ) ) 
      {
        stop("Considering the Malthus model, the length of the initial vector should be 2.")
      }else{
  stop("Considering the Logistic or Gompertz model, the length of the initial vector should be 3.")
}
}
}


    
  
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  
  
  # fitgompertz
  if(model=="Gompertz")
  {
    lsfunction <- function(x,t,y)
    {
      sum((x[1]*exp(x[2]/x[3]*(1-exp(-x[3]*t))) - y )^2)
      # sum(((x[1]*exp(-exp( (x[2]*exp(1)/x[1])*(x[3]-t) + 1) ) )-y)^2)
    }
    
  }
  
  # fitlogistic
  if(model=="Logistic")
  {
    lsfunction <- function(x,t,y)
    { 
      sum( ( ( x[1]*x[2]/( x[1] + (x[2]-x[1])*exp(- x[3]*t) ) ) - y )^2)
      # sum(((x[1]/(1+exp((4*x[2]/x[1])*(x[3]-t)+2)))-y)^2)
    }
    
  }
  
  # fitMalthus
  if(model=="Malthus")
  {
    
    lsfunction <- function(x,t,y)
    {
      sum(((x[1]*exp(x[2]*t))-y)^2)
    } 
    
  }
  
  
  ######################  FITTING THE DATA
  
  if(fitting.method == "GenSA"){
    
    fn<-function(i){
      
      m <- which(data.fit[,1] == i)
      tdata <- time[m]
      Voldata <- x[m]
      GenSA(unlist(init), lsfunction,lower=lower,upper=upper,t=tdata,y=Voldata)$par
    }
    
  }
  
  if(fitting.method == "optimr"){
    
    fn<-function(i){
      
      m <- which(data.fit[,1] == i)
      tdata <- time[m]
      Voldata <- x[m]
      # optim(init,lsfunction,method="L-BFGS-B",lower=lower,t=tdata,y=Voldata)$par
      params<-opm(unlist(init),lsfunction,lower=lower,upper=upper,method = "ALL",t=tdata,y=Voldata)
      ##### let take the params combination with minor value of SSE
      params<-summary(params,order=value)
      ##### just the convent method
      covergent<-params[params$convergence==0,]
      
      if(empty(covergent)) {
        warning("Using optmir, no methods are successful convergent!")
      }else{   
        covergent[1,1:length(init)] 
      }
      
    }
    
  }
  

  if(fitting.method == "DEoptim"){
    
    fn<-function(i){
      
      m <- which(data.fit[,1] == i)
      tdata <- time[m]
      Voldata <- x[m]
      DEoptim( lsfunction, lower=lower,upper=upper,t=tdata,y=Voldata,control = DEoptim.control(trace = F, itermax = 2000))$optim$bestmem
    }
    
  }
  
  parameters<- t(sapply(1:number_curves,fn))
  
  x<-parameters[1,]
  # set.seed(2404)
  group <- kmeans(parameters,k,iter.max = 50)
  cluster <- group$cluster
  center <- group$centers
  curves <- data.frame(ID=data$Dataset$ID, Times=data$Dataset$Time, Vol=data$Dataset$Vol, Cluster= rep(cluster,data$LenCurv))
  
  
  seed <- .Random.seed
  ###### meancurves calculating ######
  
  for(i in c(1:k))
  {
    center[i,]->x
    if(model=="Gompertz")
    {
      # meancurves[,i] <- (x[1]*exp(-exp((x[2]*exp(1)/x[1])*(x[3]-t)+1)))
      meancurves[,i] <-x[1]*exp(x[2]/x[3]*(1-exp(-x[3]*t)))
    }
    if(model=="Logistic")
    {
      # meancurves[,i] <-(x[1]/(1+exp((4*x[2]/x[1])*(x[3]-t)+2)))
      meancurves[,i] <-x[1]*x[2]/( x[1] + (x[2]-x[1])*exp(-x[3]*t) )
    }
    if(model=="Malthus")
    {
      meancurves[,i] <- x[1]*exp(x[2]*t)
    }
  }
  
 OUT<-list(parameters=parameters,cluster=cluster,center=center,meancurves=meancurves,Summary=curves,seed=seed)

  return(OUT)
}
