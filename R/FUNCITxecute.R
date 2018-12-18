#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

funcit.simo <- function(data,  k, 
                   methods=c("fitfclust", "distclust", "iterSubspace", "funclust",
                       "funHDDC", "fscm", "waveclust"),
                   seed=NULL,
                   regTime=NULL, clusters=NULL,
                   funcyCtrl=NULL, fpcCtrl=NULL,
                   parallel=FALSE, save.data=TRUE, ...){
    
    ##Check for missing arguments
    if(missing(methods))
        stop("Please select one method or methods='ALL'.")
    else if(length(methods)==1)
        if(methods=="ALL")
            methods <- 1:7
    ##Method names
    allMethods <- c("fitfclust", "distclust", "iterSubspace", "funclust",
                    "funHDDC", "fscm", "waveclust")
    
    if(is.numeric(methods))
        usedMethods <- allMethods[methods]
    else
        usedMethods <- match.arg(methods, allMethods, several.ok=TRUE)
    nrMethods <- length(usedMethods)
    
    if(missing(k))
        stop(paste(k, "is missing"))
    if(!(class(data) %in% c("matrix", "data.frame")))
        stop(paste(data,
                   "must be given in matrix or data.frame format."))
    
    ##Check if data is in the right format
    chf <- checkFormat(data)
    data <- chf$data
    reg <- chf$reg

    ##Check if funcyCtrl class is given
    if(is.null(funcyCtrl))
        funcyCtrl <- new("funcyCtrl")
    funcyCtrl@seed <- seed

    ##Convert funcyCtrl automatically to funcyCtrl if model based cluster
    ##algorithm but fpcCtrl was chosen 
    if(sum(usedMethods%in%allMethods[c(1,3,4,5,6,7)])>0 & class(funcyCtrl)=="funcyCtrl")
        funcyCtrl <- as(funcyCtrl, "funcyCtrlMbc")
    
    ##Check if fpcCtrl object if defined for eigenbasis
    if(funcyCtrl@baseType!="eigenbasis" & !is.null(fpcCtrl))
        warning("fpcCtrl is ignored since it controls only eigenbasis.")
    else if(funcyCtrl@baseType=="eigenbasis")
        fpcCtrl <- fpcCtrlCheck(fpcCtrl=fpcCtrl, data=data, reg=reg)
    
    ##Check if correct method for the given dataset was chosen.
    if(reg==0 & sum(usedMethods%in%c("fscm", "funclust", "funHDDC"))>0){
        notWork <- usedMethods[which(usedMethods%in%allMethods[-c(1:3)])]
        stop(paste("Algorithm", notWork,
                   "works only on regular data!\n Please choose one of fitfclust, distclust or iterSubspace."))}
    
    
    # check if parallel computing is needed
    if(nrMethods == 1)  parallel <- FALSE
    
    # check if parallel computing is available
    if(.Platform$OS.type!="unix" & parallel) {
    	warning("Parallel computing is only supported on Unix platforms.")
    	parallel <- FALSE
    }
    
    if(parallel) {
    	parallelFct <- parallel::mcparallel
    	coresNr <- detectCores()-1
    	options("cores"=coresNr)
    	
    } else {
    	parallelFct <- identity
    }

    
    RES <-  list()
    ##Method1--------------------
    if("fitfclust" %in% usedMethods){
        indx <- match("fitfclust",usedMethods)
        RES[[indx]] <-
            parallelFct(fitfclustWrapper(data=data, k=k, 
                                         reg=reg, regTime=regTime, fpcCtrl=fpcCtrl,
                                         funcyCtrlMbc=funcyCtrl,
                                         ...))
    }
    ##Method2----------------------
    if("distclust" %in% usedMethods){
        indx <- match("distclust", usedMethods)
        RES[[indx]] <- 
            parallelFct(distclustWrapper(data=data, k=k,
                                         reg=reg, regTime=regTime, fpcCtrl=fpcCtrl,
                                         funcyCtrl=funcyCtrl, ...))           
    }
    ##Method 3----------------------
    if("iterSubspace" %in% usedMethods){
        indx <- match("iterSubspace",usedMethods)
        RES[[indx]] <- 
            parallelFct(iterSubspaceWrapper(data=data, k=k, reg=reg, regTime=regTime,
                                            fpcCtrl=fpcCtrl,
                                            funcyCtrlMbc=funcyCtrl,  ...))
    }
    ##Method 4-----------
    if("funclust" %in% usedMethods){
    	if(!requireNamespace("Funclustering"))
    		stop("Please install package 'Funclustering' to use method 'funclust'.")
    	
        indx <- match("funclust",usedMethods)
        RES[[indx]] <-
            parallelFct(funclustWrapper(data=data, k=k, 
                                        reg=reg, regTime=regTime,
                                        funcyCtrlMbc=funcyCtrl,
                                        ...))
    }
    ##Method 5-----------
    if("funHDDC" %in% usedMethods){
    	if(!requireNamespace("funHDDC"))
    		stop("Please install package 'funHDDC' to use method 'funHDDC'.")
    	
        indx <- match("funHDDC", usedMethods)
        RES[[indx]] <-
            parallelFct(funHDDCWrapper(data=data, k=k,
                                       reg=reg, regTime=regTime,
                                       funcyCtrlMbc=funcyCtrl, ...))
    }
    ##Method 6-----------
    if("fscm" %in% usedMethods){
        indx <- match("fscm", usedMethods)
        RES[[indx]] <-
            parallelFct(fscmWrapper(data=data, k=k, reg=reg,
                                    regTime=regTime,
                                    funcyCtrlMbc=funcyCtrl, ...))
    }
    ##Method 7-----------
    if("waveclust" %in% usedMethods){
        indx <- match("waveclust", usedMethods)
        RES[[indx]] <-
            parallelFct(waveclustWrapper(data=data, k=k, reg=reg,
                                         regTime=regTime,
                                         funcyCtrlMbc=funcyCtrl, ...))
    }
    
    FRES <- new("funcyOutList")
    FRES@call <- match.call()
    
    if(parallel)
        FRES@models <- parallel::mccollect(RES)
    else
        FRES@models <- RES
    names(FRES@models) <- usedMethods

    ##Check if error appeard (only for parallel computing)----
    error <- which(sapply(FRES@models, class) == "try-error")
    if(sum(error)!=0)
        stop(paste("Method", usedMethods[error[1]], ":",
                   attributes(FRES@models[[error[1]]])$condition$message))
    
    
    allClusters <- sapply(FRES@models, function(x) x@cluster)
    allCenters <- lapply(FRES@models, function(x) x@centers)
    names(allCenters) <- colnames(allClusters) <- usedMethods
    rI <- rIMethods(methodNames=usedMethods, cls=allClusters, trueCluster=clusters)

    ##Relabel cluster output for better comparability in plots
    if(nrMethods>1){
        rel <- relabelMethods(methodNames=usedMethods, cls=allClusters,
                              ctrs=allCenters)
        allClusters <- rel$allClusters
        allCenters <- rel$allCenters
        for(i in 1:nrMethods){
            FRES@models[[i]]@cluster <- allClusters[,i]
            FRES@models[[i]]@centers <- allCenters[[i]]
            FRES@models[[i]]@correctCl <- rI[i,i]
        }
    }

    ##Warning if cluster size is smaller than 3
    smallCl <-  which(apply(allClusters, 2, function(x)
        min(table(x)))<2)

    if(length(smallCl)!=0){
        warning(paste("Method", usedMethods[smallCl],
                      "has clusters with less than 3 obervations!\n"), immediate.=TRUE)
    }
    
    accord <- accordance(cls=allClusters, relabel=FALSE)

    if(save.data)
        FRES@data <- data
    else
        FRES@data <- as.matrix(NULL)
    
    FRES@timeNr <- calcTimeNr(data, reg)
    FRES@reg <- reg
    FRES@k <- k
    FRES@methodName <- usedMethods
    FRES@allClusters <- allClusters
    FRES@randIndex <-  rI
    FRES@votedCluster <- accord$votedCluster
    FRES@accordance <- accord$accordance

    return(FRES)
}
