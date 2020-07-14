#'Extract mean summer NDVI values for buffers around lat/long points
#'
#'
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param filepath A \code{character} string stating the path to the folder where raw
#'Gimmsfiles are stored (downloaded using \code{\link[gimms]{downloadGimms}})
#'@param points A \code{dataframe} of points with columns \code{Latitude} and \code{Longitude}
#'@param buffer Positive numeric value stating the length (in meters) of the desired buffer radius.
#'Mean values of the climate variable will be calculated for each point from within this buffer
#'@param northern_summer Logical. If TRUE, calculate NDVI for the northern hemisphere summer using
#'July observations. If FALSE, calculate NDVI for southern hemisphere summer using January observations.
#'Default is TRUE
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@export
processNDVI = function(filepath, points, buffer, northern_summer = T, n_cores){

if(missing(n_cores)){
  n_cores <- parallel::detectCores() - 1
}

#### Convert coordinates to a SpatialPointsDataFrame object ####
cat('Converting coordinates to a SpatialPointsDataFrame ...\n')
points_spdf <- sp::SpatialPointsDataFrame(coords = points[, c('Longitude','Latitude')],
                                          data = points,
                                          proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +no_defs"))


#### Find all Gimms files within the specified folder ####
gimms_files <- list.files(filepath, full.names = TRUE)

if(northern_summer){
# Keep only files containing summer NDVI readings (remove winter readings)
gimms_files <- gimms_files[!grepl('0106', gimms_files)]

} else {
  gimms_files <- gimms_files[grepl('0106', gimms_files)]
}

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

if(parallel_compliant){

# Export data and libraries to the clusters
clusterExport(NULL, c('gimms_files', 'points_spdf', 'buffer'),
              envir = environment())

clusterEvalQ(cl, library(gimms))

#### Calculate mean NDVI for all points across each file in gimms_files ####
extract_gimms_data <- pbapply::pblapply(seq_along(gimms_files), function(x){

  # Extract mean values within the specified buffer range
  gimmsRaster <- gimms::rasterizeGimms(gimms_files[x], keep = c(1,2,3))
  ndvi <- raster::extract(gimmsRaster, points_spdf, buffer = buffer)
  rm(gimmsRaster)

  ndvi <- as.numeric(lapply(ndvi, mean, na.rm = TRUE))

  # Extract the Year value and return a dataframe
  year <- as.numeric(substr(basename(gimms_files[x]), 15, 18))

  data.frame(Year = year, NDVI = ndvi,
             Latitude = points$Latitude,
             Longitude = points$Longitude,
             stringsAsFactors = FALSE)
},
cl = cl)

stopCluster(cl)

} else {

  #### Use lapply if parallel loading fails ####
  extract_gimms_data <- pbapply::pblapply(seq_along(gimms_files), function(x){

    # Extract mean and sd values within the specified buffer range
    gimmsRaster <- gimms::rasterizeGimms(gimms_files[x], keep = c(1,2,3))
    ndvi <- raster::extract(gimmsRaster, points_spdf, buffer = buffer)
    rm(gimmsRaster)

    ndvi_mean <- as.numeric(lapply(ndvi, mean, na.rm = TRUE))
    ndvi_sd <- as.numeric(lapply(ndvi, sd, na.rm = TRUE))

    # Extract the Year value and return a dataframe
    year <- as.numeric(substr(basename(gimms_files[x]), 15, 18))

    data.frame(Year = year, NDVI_mean = ndvi_mean,
               NDVI_sd = ndvi_sd,
               Latitude = points$Latitude,
               Longitude = points$Longitude,
               stringsAsFactors = FALSE)
  })
}

#### Bind all outputs and return as a dataframe ####
output <- do.call(rbind, extract_gimms_data)
return(output)
}
