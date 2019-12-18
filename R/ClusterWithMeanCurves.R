#' Cluster Curves
#'
#'@description
#' Visualizes the plots regarding the fitted and clustered data by exploiting the FCM.
#'
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{MostProbableClustering.Extrapolation}}).
#' @param data CONNECTORList. (see \code{\link{DataImport}} or \code{\link{DataTruncation}}).
#' @param feature The column name reported in the AnnotationFile containing the feature  to be investigated.
#' @param title The string containing  the plot title. 
#' @param labels  The vector containing the axis names. 
#' @param save If TRUE then the following objects are saved: (i) the mean curves plot, (ii) the plots of each cluster showing the correspondive mean curve and the samples belonging to the cluster, (iii) one plot storing all the clustering plots, and (iv) the tables reporting the M, S, R and fDB indexes considering the 0, 1 and 2 derivatives.
#' @param path  The folder path where the plots will be saved. If it is missing, the plots are saved in the current working  directory.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#' 
#' @return ClusterWithMeanCurve returns a list with three objects:
#' \itemize{
#' \item{PlotsCluster:}{a list storing the growth curves plots partitioned according to the cluster membership; }
#' \item{PlotMeanCurve:}{the cluster mean curves plot; }
#' \item{spline.plots:}{a list of N plots, where N is the number of samples, showing: (i) in blue the sample curve, (ii) in red the cubic spline estimated from the FCM, (iii) in black the correspondive cluster mean curve, and finally (iv) the grey area represents the confidence interval. } 
#' }
#'
#' @details We define the following indexes for obtaining a cluster separation measure:
#' \itemize{
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
#'  \deqn{fDB_q = \frac{1}{G} \sum_{k=1}^G \max_{h \neq k} { R_{hk} }. }
#'  }
#' }
#' Where the proximities measures choosen is defined as follow
#'  \deqn{D_q(f,g) = \sqrt( \integral | f^{(q)}(s)-g^{(q)}(s) |^2 ds ), d=0,1,2}
#' with f and g are two curves and f^{(q)} and g^{(q)} are their q-th derivatives. Note that for q=0, the equation becomes the distance induced by the classical L^2-norm.
#' 
#' 
#' @seealso MostProbableClustering.Extrapolation, StabilityAnalysis, Spline.plots.
#'
#' @examples
#' 
#' @import ggplot2  gridExtra
#' @importFrom cowplot plot_grid add_sub ggdraw
#' @importFrom knitr kable
#' @export
ClusterWithMeanCurve<-function(clusterdata, data, feature ,title="", labels=c("","") ,save=FALSE,path=NULL )
{
axis.x<-labels[1]
axis.y<-labels[2]

clusterdata.info<-clusterdata$Cl.Info
clusterdata<-clusterdata$FCM

    if(!is.null(clusterdata$fit))
    {
      G<-length(clusterdata$prediction$meancurves[1,])
      
      clusterdata$cluster$cluster.member -> classes
      clusterdata$prediction$meancurves -> meancurves
      clusterdata$fit$grid -> grid
      symbols<-clusterdata$cluster$cluster.names
      
    }else{
      warning("An object returned by Cluster_choice is required.",immediate. = TRUE)
}
    
    
    

    classificate <- rep(classes,data$LenCurv)
    
    curves <- data.frame(Times=data$Dataset$Time,Vol=data$Dataset$Vol,ID=data$Dataset$ID,Cluster=classificate, Info=rep(t(data$LabCurv[feature]),data$LenCurv))
    

    
    # cut the meancurves at the growth curves' maximum time
    
    meancurves_truncated<-c()
    time3<-c()
    cluster<-c()
    
    for(clust in 1:G)
    {
      time2<-sort(unique(curves[curves$Cluster==clust,]$Times))
      m<-meancurves[,clust][grid%in%time2]
      time3<-c(time3,grid[grid%in%time2])
      meancurves_truncated<-c(meancurves_truncated,m)
      cluster<-c(cluster,rep(symbols[clust],length(grid[grid%in%time2])))
    }
    cl.names<-clusterdata$cluster$cluster.names
    

    plot_data<-data.frame(time=time3,means=meancurves_truncated,clusters=cluster)
    
    PlotMeanCurve<-ggplot()+
      geom_line(data=plot_data, aes(x=time,y=means,group=clusters,linetype= as.factor(clusters)),size=1 )+  
      scale_linetype_manual(values =1:G ,limits=sort(symbols),breaks=sort(symbols),name="Cluster") +
      labs(title=title, x=axis.x, y = axis.y,linetype="Cluster")+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())

    
      col<- as.character(unique(curves$Info))
      col1<- data$ColFeature
      plots<-list()
      plots_spline<-list()
      ymax<-max(plot_data$means,curves$Vol)
      xmax<-max(grid)
      
      esse<-clusterdata.info$Coefficents$esse
      essed1<-clusterdata.info$Deriv.Coefficents$esse
      essed2<-clusterdata.info$Deriv2.Coefficents$esse
      fDB<-clusterdata.info$Coefficents$fDB.index
      fDB1<-clusterdata.info$Deriv.Coefficents$fDB.index
      fDB2<-clusterdata.info$Deriv2.Coefficents$fDB.index
      errei<-clusterdata.info$Coefficents$errei
      erreid1<-clusterdata.info$Deriv.Coefficents$errei
      erreid2<-clusterdata.info$Deriv2.Coefficents$errei
      
      emme<-clusterdata.info$Coefficents$emme
      row.names(emme)=paste("Cluster",row.names(emme))
      colnames(emme)=paste("Cluster",colnames(emme))
      
      
      ##### Build the S data frame  
      S.indexes<-data.frame(S=esse,S_1=essed1,S_2=essed2)
      row.names(S.indexes)=paste("Cluster",row.names(S.indexes))
      cat("\n######## S indexes ############\n")
      print(kable(S.indexes))
      cat("\n##############################################################")
      
      ##### Build the M data frame
      cat("############################################################## \n######## M indexes ############\n")
      print(kable(emme))
      cat("\n##############################################################")
      
      ##### Build the R data frame  
      R.indexes<-data.frame(R=errei,R_1=erreid1,R_2=erreid2)
      row.names(R.indexes)=paste("Cluster",row.names(R.indexes))
      cat("\n######## R indexes ############\n")
      print(kable(R.indexes))
      cat("\n##############################################################")
      
      ##### Build the fDB data frame 
      fDB.indexes<-data.frame(fDB=fDB,fDB_1=fDB1,fDB_2=fDB2)
      cat("\n######## fDB indexes ############\n")
      print(kable(fDB.indexes))
      cat("\n##############################################################")
      ###########################
      
      
      for(i in 1:G)
      {
        order(symbols)[i]->index.symb # sorting from the lower to the higher mean curve w.r.t. the zero axis
        
        plots[[i]]<- ggplot()+
          geom_line(data=plot_data[plot_data$clusters==symbols[index.symb],], aes(x=time,y=means,linetype= as.factor(clusters)),size = 1.2 )+
          scale_linetype_manual(values =1:G ,limits=sort(symbols),breaks=sort(symbols),name="Cluster") +
          labs(title=paste("Cluster",symbols[index.symb]), x=axis.x, y = axis.y)+
          geom_line(data = curves[curves$Cluster==index.symb,],aes(x=Times,y=Vol,group=ID,color=factor(Info)))+
          scale_colour_manual(values = col1,limits=col,breaks=col,name=feature)+
          theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())+ ylim(0,ymax)+xlim(0,xmax)

      }
      
      ### grouping all the plot per cluster and print it
       allCl.plot<-plot_grid(plotlist = plots) 
       
       h <- length(clusterdata$fit$parameters$Lambda[1,])
       p <- length(clusterdata$fit$parameters$Lambda[,1])
       
       plots[["ALL"]]<-ggdraw(add_sub(allCl.plot, paste("Other parameters: p = ", p, ", h = ", h, ", G = ", G  ,sep = ""))  )
       
       print(plots[["ALL"]])
       
       ### 
       
       if(save)
       {
         if(is.null(path)) path <- getwd()
         
         for(i in 1:G)
         {
           #### Save the clusters and mean curves as pdf
           
           ggsave(filename = paste(symbols[i],"Cluster.pdf",sep=""),plot=plots[[paste(symbols[i],"Cluster")]],width=29, height = 20, units = "cm",scale = 1,path = path)
         }
           ggsave(filename = "ALLCluster.pdf",plot=plots$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
           ggsave(filename = "MeanCurve.pdf",plot=PlotMeanCurve,width=29, height = 20, units = "cm",scale = 1,path = path)
           
           #### Save the fDB indexes as pdf
           # Transform tables into grobs
           tt3 <- ttheme_default(
             rowhead=list(fg_params=list(fontface="bold"))
           )
           
           Table1 <- arrangeGrob(top="fDB indexes",tableGrob(fDB.indexes,theme=tt3))
           Table2 <- arrangeGrob(top="R indexes",tableGrob(R.indexes,theme=tt3))
           Table3 <- arrangeGrob(top="S indexes",tableGrob(S.indexes,theme=tt3))
           Table4 <- arrangeGrob(top="M indexes",tableGrob(emme,theme=tt3))

           ml <- marrangeGrob(list(Table1, Table2, Table3, Table4),nrow=4,ncol=1)
           
           ggsave(paste0(path,"Indexes.pdf"),width=29, height = 20, units = "cm", ml)
           
       }
         
######## Let's plot the spline fitting for each sample 
         
    if(!is.null(clusterdata$fit))
     {
      
      allsplineEstimation<-curve_prediction(cluster=classes ,object=clusterdata$prediction)

      out<-list(plotMeanCurve=PlotMeanCurve,plotsCluster=plots,spline.plots=allsplineEstimation)
    }
    else{
      out<-list(plotMeanCurve=PlotMeanCurve,plotsCluster=plots)
    }

  return(out)
}

curve_prediction<-function(cluster,object)
{
  index <- 1:length(table(object$data$curve)) #pdx
  
  timeIndx <- object$data$timeindex
  curveIndx <- object$data$curve
  
 plotCreation<-function(i,cluster,object)
 {
    cl<-cluster[i]
    grid <- object$grid
    upci <- object$upci[i,]
    uppi <- object$uppi[i,]
    lowci <- object$lowci[i,]
    lowpi <- object$lowpi[i,]
    gpred <- object$gpred[i,]
    meancurves <- (object$mean)[,cl]

        yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
    
    data.ggplot<-data.frame(grid=grid,upci=upci,uppi=uppi,lowci=lowci,lowpi=lowpi,gpred=gpred,meancurves=meancurves)
    data.real<-data.frame(time=grid[timeIndx[curveIndx==i]],vol=object$data$x[curveIndx==i] )
    
    gpl<-
      ggplot()+
      geom_ribbon(data=data.ggplot,aes(x=grid,ymin=lowci, ymax=upci), alpha=0.1)+
      geom_line(data=data.ggplot,aes(x=grid,y=gpred,linetype="Spline estimated",col="Spline estimated"))+
      geom_line(data=data.ggplot,aes(x=grid,y=meancurves,linetype="Cluster mean",col="Cluster mean"))+
      #geom_line(data=data.real,aes(x=time,y=gpredF,linetype="SplineFiltering",col="SplineFiltering"))+
      geom_line(data=data.real,aes(x=time,y=vol,col="Real points",linetype="Real points"))+
      geom_point(data=data.real,aes(x=time,y=vol),col="blue")+
      labs(title=paste("Sample ",i), x="Time", y="Growth value")+
      scale_colour_manual(values = c("black","red","blue") ,limits =c("Cluster mean","Spline estimated","Real points"),breaks= c("Cluster mean","Spline estimated","Real points") , name=" ")+
      guides(
             linetype = FALSE,
             colour = guide_legend(override.aes = list(linetype = c("solid","dashed","dashed"))))+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.key.width = unit(1, "cm"))
    
    
    return(gpl)
 }
 
 plot_list<-lapply(index,plotCreation,cluster,object)
 names(plot_list)<-paste(paste("Sample ",index))

 return(plot_list)
}

