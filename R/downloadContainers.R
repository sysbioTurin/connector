#' @title Download for the first time all containers embedded in the workflows
#' @description This is a function that prepares the docker environment to be used for the first time the conector package is installed.
#' @param containers.file, a character string with the name of the file which indicate which are the initial set of containers to be downloaded. If NULL then the set is given by a file called "containersNames.txt" located in the folder inst/Containers of connector package.
#' @author Beccuti Marco, Pernice Simone
#'
#' @examples
#'\dontrun{
##'     #running runDocker
#'      downloadContainers()
#'
#' }
#' @export
#' 
downloadContainers <- function(containers.file=NULL){
    if(is.null(containers.file)){
        containers.file=paste(path.package(package="connector"),"Containers/containersNames.txt",sep="/")
        containers <-read.table(containers.file,header=TRUE,row.names = 1)
    }
    else{
        containers <-read.table(containers.file,header=TRUE,row.names = 1)
        #writeLines(containers, paste(path.package(package="connector"),"Containers/containersNames.txt",sep="/"))
        
    }


    for(i in 1:dim(containers)[1]){
        system(paste("docker pull ",containers[i,1], sep=""))
    }
    
 
}
