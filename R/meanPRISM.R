#'Extract mean monthly PRISM climate values for buffers around lat/long points
#'
#'
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param points A \code{dataframe} of points with columns \code{Latitude} and \code{Longitude}
#'@param var_name A \code{character} string specifying the name of the climate variable. If not provided,
#'the variable's column name in the returned \code{dataframe} will be \code{mean.value}. Note that some
#'particular point * month measurements may be missing, and so \code{NA}s will be returned for these. A
#'warnign message will be printed if this occurs.
#'@param buffer Positive numeric value stating the length (in meters) of the desired buffer radius.
#'Mean values of the climate variable will be calculated for each point from within this buffer
#'@param filepath A \code{character} string stating the path to the folder where unzipped
#'PRISM files are stored (downloaded using \code{\link[prism]{get_prism_monthlys}})
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'@param delete_files Logical describing whether or not to delete raw PRISM data folders once the job
#'is complete. Default is \code{FALSE}
#'
#'@export
meanPRISM = function(points, var_name, buffer, filepath, n_cores, delete_files){

  if(missing(n_cores)){
    n_cores <- detectCores() - 1
  }

  if(missing(delete_files)){
    delete_files <- FALSE
  }

  if(missing(var_name)){
    rename_output <- FALSE
  } else{
    rename_output <- TRUE
  }

  #### Create a SpatialPointsDataFrame from the specified coordinates ####
  cat('Converting coordinates to a SpatialPointsDataFrame ...\n')
  points_spdf <- sp::SpatialPointsDataFrame(coords = points[, c('Longitude','Latitude')],
                                            data = points,
                                            proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +no_defs"))

  #### If n_cores > 1, check parallel library loading ####
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/purrr$", "", system.file(package = "purrr"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(purrr)), silent = TRUE)

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
    warning('Parallel loading failed, calculations may crash!')
  }

  if(parallel_compliant){
    cat('Parallel loading successful', "\n", sep = "")

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('points', 'points_spdf', 'buffer', 'filepath'),
                  envir = environment())

    #Export necessary libraries
    clusterEvalQ(cl, library(prism))
    clusterEvalQ(cl, library(sp))
    clusterEvalQ(cl, library(raster))
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(purrr))
    clusterEvalQ(cl, library(stringr))


  #### Extract mean values around each point for each raster in the downloaded PRISM data ####
  cat('Processing raster files in parallel...', "\n", sep = "")

  # Set the directory for prism access
  options(prism.path = filepath)

  extracted.means <- pbapply::pblapply(seq_len(length(prism::ls_prism_data()[,1])), function(x){

    # Make sure the directory is set on each cluster
    options(prism.path = filepath)

    # Extract the prism raster stack file (check if it exists first!)
    test_RS <- prism::ls_prism_data(absPath = TRUE)

    if(file.exists(test_RS$abs_path[x])){
    RS <- prism::prism_stack(prism::ls_prism_data()[x, 1])

    # Make sure to set the coordinate reference system
    proj4string(RS) <- sp::CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

    # Extract mean climate variable for each point using the supplied buffer
    clim.buffer <- as.data.frame(raster::extract(RS, points_spdf,
                                                 buffer = buffer, df = TRUE)) %>%
      purrr::set_names(c('site', 'value')) %>%
      dplyr::group_by(site) %>%
      dplyr::summarise(mean.value = mean(value, na.rm = T)) %>%
      dplyr::mutate(Year = as.numeric(substr(stringr::str_extract(names(RS@layers[[1]]),
                                                                  "\\d+(?=_[a-zA-Z]+.+$)"),
                                             start = 1, stop = 4)),
                    Month = as.numeric(substr(stringr::str_extract(names(RS@layers[[1]]),
                                                                   "\\d+(?=_[a-zA-Z]+.+$)"),
                                              start = 5, stop = 6))) %>%
      dplyr::bind_cols(points) %>%
      dplyr::select(-site)
    } else {
      # Return NAs for climate variables if file does not exist
      clim.buffer <- data.frame(mean.value = 'NA',
                                Year = as.numeric(substr(stringr::str_extract(test_RS$files[x],
                                                                              "\\d+(?=_[a-zA-Z]+.+$)"),
                                                         start = 1, stop = 4)),
                                Month = as.numeric(substr(stringr::str_extract(test_RS$files[x],
                                                                               "\\d+(?=_[a-zA-Z]+.+$)"),
                                                          start = 5, stop = 6)),
                                Latitude = points$Latitude,
                                Longitude = points$Longitude)
    }
    clim.buffer
  },
  cl = cl)
  stopCluster(cl)
  }

  if(!parallel_compliant){
    cat('Parallel loading failed or not specified, using a single core instead (may take a while)...', "\n", sep = "")

    # Set the directory for prism access
    options(prism.path = filepath)

    extracted.means <- pbapply::pblapply(seq_len(length(prism::ls_prism_data()[,1])), function(x){

      options(prism.path = filepath)

      # Extract the prism raster stack file (check if it exists first!)
      test_RS <- prism::ls_prism_data(absPath = TRUE)

      if(file.exists(test_RS$abs_path[x])){
        RS <- prism::prism_stack(prism::ls_prism_data()[x, 1])

      # Make sure to set the coordinate reference system
      proj4string(RS) <- sp::CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

      # Extract mean climate variable for each point using the supplied buffer
      clim.buffer <- as.data.frame(raster::extract(RS, points_spdf,
                                                   buffer = buffer, df = TRUE)) %>%
        purrr::set_names(c('site', 'value')) %>%
        dplyr::group_by(site) %>%
        dplyr::summarise(mean.value = mean(value, na.rm = T)) %>%
        dplyr::mutate(Year = as.numeric(substr(stringr::str_extract(names(RS@layers[[1]]),
                                                                    "\\d+(?=_[a-zA-Z]+.+$)"),
                                               start = 1, stop = 4)),
                      Month = as.numeric(substr(stringr::str_extract(names(RS@layers[[1]]),
                                                                     "\\d+(?=_[a-zA-Z]+.+$)"),
                                                start = 5, stop = 6))) %>%
        dplyr::bind_cols(points) %>%
        dplyr::select(-site)
      } else {
        clim.buffer <- data.frame(mean.value = 'NA',
                                  Year = as.numeric(substr(stringr::str_extract(test_RS$files[x],
                                                                                "\\d+(?=_[a-zA-Z]+.+$)"),
                                                           start = 1, stop = 4)),
                                  Month = as.numeric(substr(stringr::str_extract(test_RS$files[x],
                                                                                 "\\d+(?=_[a-zA-Z]+.+$)"),
                                                            start = 5, stop = 6)),
                                  Latitude = points$Latitude,
                                  Longitude = points$Longitude)
      }
      clim.buffer
    })
  }

  # Delete folders and files in filepath if specified
  if(delete_files){
    unlink(filepath, recursive = T)
    cat('Done processing, raw PRISM files deleted')
  } else{
    cat('Done processing, raw PRISM files remain in ', filepath, sep = "")
  }

  # Return a dataframe of all mean observations
  if(rename_output){
    output <- do.call(rbind, extracted.means)
    colnames(output) <- c(var_name, 'Year', 'Month', 'Latitude', 'Longitude')
  } else{
    output <- do.call(rbind, extracted.means)
  }

  # Make sure that all cols are numeric, will return a useful warning if are NAs present
  output %>%
    dplyr::mutate_all(as.numeric) -> output

  return(output)
}
