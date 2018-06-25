#'Extract mean altitude within a 40km buffer for lat and long coordinates
#'
#'This function calculates mean elevation with 40km buffers of specified points. Note,
#'the function will temporarily save the global digital elevation map in a new folder
#'in the working directory, so please ensure there is sufficient memory (around 37MB) available
#'
#'@importFrom magrittr %>%
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param points A \code{dataframe} containing the coordinates, with colnames 'Latitude' and
#''Longitude'
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'@param delete_files Logical describing whether or not to delete raw elevation map once the job
#'is complete. Default is \code{FALSE}
#'
#'@export
processAltitudes = function(points, n_cores, delete_files){

  if(missing(n_cores)){
    n_cores <- parallel::detectCores() - 1
  }

  if(missing(delete_files)){
    delete_files <- FALSE
  }

  #### Download the geotif raster of global digital elevation ####
  cat('Downloading the global digital elevation map ...\n')
  altitude_raster <- raster::raster("https://asterweb.jpl.nasa.gov/images/GDEM-10km-BW.tif")

  # Save the raster in the working directory as a .grd object for faster computations
  cat('Saving the global digital elevation map as a .grd object in the working directory ...\n')

  dir.create('./Elevation_map')
  raster::writeRaster(altitude_raster, "./Elevation_map/altitude_raster.grd")
  altitude_raster <- raster::raster("./Elevation_map/altitude_raster.grd")

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
    test_load1 <- try(clusterEvalQ(cl, library(raster)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/raster$", "", system.file(package = "raster"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(raster)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(raster)), silent = TRUE)

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
    warning('Parallel loading failed, calculations will take a while')
  }

  cat('Extracting mean elevation for each coordinate (40km buffer) ...\n')

  if(parallel_compliant){
  #### Export necessary data / libraries ####
  clusterExport(NULL, c('points_spdf', 'altitude_raster', 'points'),
                envir = environment())

  clusterEvalQ(cl, library(dplyr))

  #### Extract mean altitude (40km buffer) for each site ####
  altitudes <- pbapply::pblapply(seq_len(nrow(points)), function(i){
    raster::extract(altitude_raster, points_spdf[i,],
                    buffer = 40000, fun = mean,
                    df = TRUE) %>%
      dplyr::bind_cols(points[i,]) %>%
      purrr::set_names('ID', 'Altitude', 'Latitude', 'Longitude') %>%
      dplyr::select(-ID)

  }, cl = cl)
  stopCluster(cl)

  } else {
    #### Use lapply if parallel loading fails ####
    altitudes <- pbapply::pblapply(seq_len(nrow(points)), function(i){
      raster::extract(altitude_raster, points_spdf[i,],
                      buffer = 40000, fun = mean,
                      df = TRUE) %>%
        dplyr::bind_cols(points[i,]) %>%
        purrr::set_names('ID', 'Altitude', 'Latitude', 'Longitude') %>%
        dplyr::select(-ID)

    })

  }

  #### Remove the saved altitude map files ####
  # Delete folders and files in filepath if specified
  if(delete_files){
    unlink('./Elevation_map/', recursive = TRUE)
    cat('Done processing, raw elevation map deleted')

  } else{
    cat('Done processing, raw elevation map remains in the ./Elevation_map/ directory')
  }

  #### Bind coordinates and altitudes into a single dataframe and return ####
  altitude_df <- do.call(rbind, altitudes)
  return(altitude_df)
}

