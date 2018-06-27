#'Calculate turnover in predicted species interaction networks
#'
#'This function takes a network based on predicted species' interactions and
#'returns the beta diversity metric \emph{B}'os,
#'defined by Poisot et al. 2012. Please note, this function was adapted from code
#'provided by Timothee Poisot in the \code{\link[betalink]{betalink}} package. For details,
#'see \code{\link[betalink]{beta_os_prime}}
#'
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param adjacency_list The input data should be a \code{list} of
#'\code{\link[igraph]{graph.adjacency}} matrices, ideally as returned by
#'\code{\link[MRFcov]{predict_MRFnetworks}} using \code{metric == 'adjacency'}
#'
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@details This function calculates \emph{B}'os,
#'i.e. the distance between a local network's realized interactions and the
#'potential interactions (calculated by integrating interactions from all local
#'networks into a regional \emph{metaweb}). Values close to \code{0} suggest there is a low local
#'selection of interactions (nearly all regional interactions are found in the local
#'network), whereas values closer to \code{1} reflect local selection of species' interactions
#'
#'@return A \code{vector} containing predicted \emph{B}'os values for each network in \code{adjacency_list}
#'
#'@references Poisot, T., Canard, E., Mouillot, D., Mouquet, N. & Gravel, D. (2012)
#'The dissimilarity of species interaction networks. Ecology Letters, 15, 1353-1361.
#'
#'@seealso \code{\link[betalink]{beta_os_prime}}\cr
#'\code{\link[betalink]{metaweb}}
#'@export
#'
#'
networkBetaOS <- function(adjacency_list, n_cores){

#Make sure that betalink is installed
  if(!require(betalink)){
    devtools::install_github('PoisotLab/betalink')
  }

#Convert the list of adjacency matrices to a metaweb
adjacency_list_no.na <- adjacency_list[lapply(adjacency_list, length) > 0]
metaweb <- betalink::metaweb(adjacency_list_no.na)

if(missing(n_cores)){
  n_cores <- parallel::detectCores() - 1
}

#### If n_cores > 1, check parallel library loading ####
if(n_cores > 1){
  #Initiate the n_cores parallel clusters
  cl <- makePSOCKcluster(n_cores)
  setDefaultCluster(cl)

  #### Check for errors when directly loading a library on each cluster ####
  test_load1 <- try(clusterEvalQ(cl, library(betalink)), silent = TRUE)

  #If errors produced, iterate through other options for library loading
  if(class(test_load1) == "try-error") {

    #Try finding unique library paths using system.file()
    pkgLibs <- unique(c(sub("/betalink$", "", system.file(package = "betalink"))))
    clusterExport(NULL, c('pkgLibs'), envir = environment())
    clusterEvalQ(cl, .libPaths(pkgLibs))

    #Check again for errors loading libraries
    test_load2 <- try(clusterEvalQ(cl, library(betalink)), silent = TRUE)

    if(class(test_load2) == "try-error"){

      #Try loading the user's .libPath() directly
      clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
      test_load3 <- try(clusterEvalQ(cl, library(betalink)), silent = TRUE)

      if(class(test_load3) == "try-error"){

        #Give up and use lapply instead!
        parallel_compliant <- FALSE
        stopCluster(cl)

      } else {
        parallel_compliant <- TRUE
      }

    } else {
      parallel_compliant <- TRUE
    }

  } else {
    parallel_compliant <- TRUE
  }
} else {
  #If n_cores = 1, set parallel_compliant to FALSE
  parallel_compliant <- FALSE
  warning('Parallel loading failed')
}

if(parallel_compliant){

  #Export necessary data and variables to each cluster
  clusterExport(NULL, c('adjacency_list', 'metaweb'),
                envir = environment())

  #Export necessary libraries
  clusterEvalQ(cl, library(betalink))

  #Calculate B'os for each local network
  os_prime <- unlist(pbapply::pblapply(seq_along(adjacency_list), function(x){
    if(is.null(adjacency_list[[x]])){
      'NA'
    } else{
      betalink::betalink(adjacency_list[[x]], metaweb)$OS
    }
  }, cl = cl))
  stopCluster(cl)

} else{
  os_prime <- unlist(lapply(seq_along(adjacency_list), function(x){
    if(is.null(adjacency_list[[x]])){
      'NA'
    } else{
      betalink::betalink(adjacency_list[[x]], metaweb)$OS
    }
}))
}
return(os_prime)
}
