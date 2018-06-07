#'Use hierarchical clustering to group sites into biogeographical regions
#'based on similarity in avian species composition
#'
#'@importFrom magrittr %>%
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores
#'
#'@param dist_weights Optional \code{list} giving relative weights (as integers)
#'for the six different distance matrices when creating the final weighted matrix. The
#'six matrices are: species composition, latitude, longitude, mean annual precipitation,
#'mean minimum temperature and mean maximum temperature. Default is to only use
#'species composition, latitude, longitude with equal weights
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@details Distance matrices are created for each of the six site-level parameters. Species
#'composition is calculated by considering the presence-absence of species at each site, where
#'a presence is recorded only if that species has been observed in that site on at least five
#'of the seven total BBS surveys. A mean distance is then calculated from the six different matrices
#'using the specified weights in \code{dist_weights}. Hierachical clustering is then performed and the
#'resulting dendrogram is cut to group sites into broader biogeographical regions. Separate grouping
#'is performed for each flyway in \code{\link{Site.descriptors}}
#'
#'@return A \code{dataframe} that contains the \code{\link{Site.descriptors}} along with a new column
#'describing the bioregions
#'
#'@export
createRegions = function(dist_weights, n_cores){

  #### Function to identify optimal heights for cutting dendrograms ####
  best.cutree = function(hc, min, max, loss = FALSE){
      max <- min(max, length(hc$height))
      inert.gain <- rev(hc$height)
      intra <- rev(cumsum(rev(inert.gain)))
      relative.loss <- intra[min:(max)] / intra[(min - 1):(max - 1)]
      best <- which.min(relative.loss)
      names(relative.loss) <- min:max

      if (loss)
        relative.loss
      else
        best + min - 1
  }

#### Function to merge outlier communities with neighboring communities ####
  replace.outliers = function(outlier_region, regions_dat, clim_dat){
  find_nearest = data.frame(Region = regions_dat) %>%
    dplyr::bind_cols(clim_dat[, c('Latitude', 'Longitude')]) %>%
    dplyr::group_by(Region) %>%
    dplyr::summarise_all(funs(median(.))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Latitude, Longitude) %>%
    dplyr::mutate(Diff.lat = c(0, diff(Latitude)),
                  Diff.long = c(0, diff(Longitude))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Diff.tot = mean(c(Diff.lat, Diff.long), na.rm = T))

  if(find_nearest[find_nearest$Region == outlier_region, 'Diff.tot'] == 0){
    replacement <- as.numeric(find_nearest[2, 'Region'])

  } else {
    row_index <- which(find_nearest$Region == outlier_region)
    nearest <- which.min(c(as.numeric(find_nearest[row_index - 1, 'Region']),
                           as.numeric(find_nearest[row_index + 1, 'Region'])))
    replacement <- as.numeric(nearest)
  }
  return(replacement)
  }

  #### Set default values####
  if(missing(dist_weights)){
    # Default is to give equal weights to species composition, latitude and longitude
    # and not include climate variables
    dist_weights <- c(100, 100, 100, 0, 0, 0)
  }

  if(missing(n_cores)){
    n_cores <- parallel::detectCores() - 1
  }

  #### Load datasets from the package ####
  data("BBS.occurrences")
  data("Site.descriptors")

  # Gather names of flyways
  flyways <- as.character(unique(Site.descriptors$Site.group))

  #### If n_cores > 1, check parallel library loading ####
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(dplyr)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/dplyr$", "", system.file(package = "dplyr"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(dplyr)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(dplyr)), silent = TRUE)

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
  }

  if(parallel_compliant){

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('BBS.occurrences', 'Site.descriptors'),
                  envir = environment())

    #Export necessary functions to each cluster
    clusterExport(NULL, c('best.cutree'), envir = environment())
    clusterExport(NULL, c('replace.outliers'), envir = environment())

    #Export necessary libraries
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(ade4))

  #### Cluster sites into biogeographic regions for each flyway ####
  cluster_sites <- pbapply::pblapply(seq_len(length(flyways)), function(x){
  flyway_rows <- which(Site.descriptors$Site.group == flyways[x])

  # Calculate mean climate variables
  clim_dat = Site.descriptors %>%
    dplyr::filter(Site.group == flyways[x]) %>%
    dplyr::select(Latitude, Longitude, Tot.precip.yr,
                  Mean.min.temp.yr:Mean.max.temp.yr) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Latitude, Longitude) %>%
    dplyr::summarise_all(funs(mean(., na.rm = TRUE))) %>%
    dplyr::ungroup()

  clim_dat[is.na(clim_dat)] <- 0

  # Extract species binary data for each unique site
  sp_binary_dat = BBS.occurrences[flyway_rows, ] %>%
    dplyr::bind_cols(Site.descriptors[flyway_rows, c('Latitude', 'Longitude')]) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Latitude, Longitude) %>%
    dplyr::summarise_all(funs(sum(.))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Latitude, -Longitude)

  sp_binary_dat[sp_binary_dat >= 5] <- 1
  sp_binary_dat <- sp_binary_dat[, colSums(sp_binary_dat) > 0]

  #Prep species binary dataset for distance calculation
  sp.dist.prep <- ade4::prep.binary(data.frame(sp_binary_dat),
                                    col.blocks = ncol(sp_binary_dat))

  #Gather species composition and climate datasets into a 'ktab.list.df' object
  ktab.regions <- ade4::ktab.list.df(list(sp.dist.prep, clim_dat),
                                     w.col = list(rep(1 / ncol(sp.dist.prep),
                                                      ncol(sp.dist.prep)),
                                                  rep(1 / (ncol(clim_dat)),
                                                      (ncol(clim_dat)))))

  #Calculate site pairwise distances for each variable in ktab.regions
  dist.regions <- ade4::ldist.ktab(ktab.regions, type = c('B', 'Q'),
                                   option = c('noscale', 'scaledBYsd'))

  #Compute weighted mean pairwise distances between sites
  dist.regions.weighted <- Reduce(`+`, Map(`*`, dist.regions, dist_weights))

  #Cluster the distance matrix using hierarchical clustering
  hc_clim <- stats::hclust(dist.regions.weighted, method = "complete")

  #Calculate optimal number of clusters based on the relative loss of inertia
  cut_height <- best.cutree(hc_clim,
                            min = ceiling(max(hc_clim$height) / 3),
                            max = max(hc_clim$height))

  #Group sites into hierarchical regions
  regions_dat <- stats::cutree(hc_clim, h = cut_height)

  if(length(unique(regions_dat)) > 4){
    regions_dat <- stats::cutree(hc_clim, k = 4)

    # Merge outlier regions (less than 10 sites) with largest clusters
    if(any(table(regions_dat) < 10)){
      outlier_regions <- as.numeric(names(which(table(regions_dat) < 10)))

      index_replacement <- vector()
      for(i in 1:length(outlier_regions)){
        index_replacement[i] <- replace.outliers(outlier_region = outlier_regions[i],
                                                 regions_dat = regions_dat,
                                                 clim_dat = clim_dat)
        regions_dat[regions_dat == outlier_regions[i]] <- as.numeric(index_replacement[i])
      }
    }

  } else if(length(unique(regions_dat)) < 3){
    regions_dat <- stats::cutree(hc_clim, k = 4)

    # Merge outlier regions (less than 10 sites) with largest clusters
    if(any(table(regions_dat) < 10)){
      outlier_regions <- as.numeric(names(which(table(regions_dat) < 10)))

      index_replacement <- vector()
      for(i in 1:length(outlier_regions)){
        index_replacement[i] <- replace.outliers(outlier_region = outlier_regions[i],
                                                 regions_dat = regions_dat,
                                                 clim_dat = clim_dat)
        regions_dat[regions_dat == outlier_regions[i]] <- as.numeric(index_replacement[i])
      }
    }

  }

  regions_dat = data.frame(Flyway.region = regions_dat) %>%
    dplyr::bind_cols(clim_dat[, c('Latitude', 'Longitude')]) %>%
    dplyr::ungroup()

  regions_dat$Site.group <- flyways[x]
  regions_dat
  }, cl = cl)
  stopCluster(cl)

  } else {
    #### If parallel loading fails, use lapply instead ####
    cluster_sites <- lapply(seq_len(length(flyways)), function(x){
      flyway_rows <- which(Site.descriptors$Site.group == flyways[x])

      # Calculate mean climate variables
      clim_dat = Site.descriptors %>%
        dplyr::filter(Site.group == flyways[x]) %>%
        dplyr::select(Latitude, Longitude, Tot.precip.yr,
                      Mean.min.temp.yr:Mean.max.temp.yr) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Latitude, Longitude) %>%
        dplyr::summarise_all(funs(mean(., na.rm = TRUE))) %>%
        dplyr::ungroup()

      clim_dat[is.na(clim_dat)] <- 0

      # Extract species binary data for each unique site
      sp_binary_dat = BBS.occurrences[flyway_rows, ] %>%
        dplyr::bind_cols(Site.descriptors[flyway_rows, c('Latitude', 'Longitude')]) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Latitude, Longitude) %>%
        dplyr::summarise_all(funs(sum(.))) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Latitude, -Longitude)

      sp_binary_dat[sp_binary_dat >= 5] <- 1
      sp_binary_dat <- sp_binary_dat[, colSums(sp_binary_dat) > 0]

      #Prep species binary dataset for distance calculation
      sp.dist.prep <- ade4::prep.binary(data.frame(sp_binary_dat),
                                        col.blocks = ncol(sp_binary_dat))

      #Gather species composition and climate datasets into a 'ktab.list.df' object
      ktab.regions <- ade4::ktab.list.df(list(sp.dist.prep, clim_dat),
                                         w.col = list(rep(1 / ncol(sp.dist.prep),
                                                          ncol(sp.dist.prep)),
                                                      rep(1 / (ncol(clim_dat)),
                                                          (ncol(clim_dat)))))

      #Calculate site pairwise distances for each variable in ktab.regions
      dist.regions <- ade4::ldist.ktab(ktab.regions, type = c('B', 'Q'),
                                       option = c('noscale', 'scaledBYsd'))

      #Compute weighted mean pairwise distances between sites
      dist.regions.weighted <- Reduce(`+`, Map(`*`, dist.regions, dist_weights))

      #Cluster the distance matrix using hierarchical clustering
      hc_clim <- stats::hclust(dist.regions.weighted, method = "complete")

      #Calculate optimal number of clusters based on the relative loss of inertia
      cut_height <- best.cutree(hc_clim,
                                min = ceiling(max(hc_clim$height) / 3),
                                max = max(hc_clim$height))

      #Group sites into hierarchical regions
      regions_dat <- stats::cutree(hc_clim, h = cut_height)

      if(length(unique(regions_dat)) > 4){
        regions_dat <- stats::cutree(hc_clim, k = 4)

        # Merge outlier regions (less than 10 sites) with largest clusters
        if(any(table(regions_dat) < 10)){
          outlier_regions <- as.numeric(names(which(table(regions_dat) < 10)))

          index_replacement <- vector()
          for(i in 1:length(outlier_regions)){
            index_replacement[i] <- replace.outliers(outlier_region = outlier_regions[i],
                                                     regions_dat = regions_dat,
                                                     clim_dat = clim_dat)
            regions_dat[regions_dat == outlier_regions[i]] <- as.numeric(index_replacement[i])
          }
        }

      } else if(length(unique(regions_dat)) < 3){
        regions_dat <- stats::cutree(hc_clim, k = 4)

        # Merge outlier regions (less than 10 sites) with largest clusters
        if(any(table(regions_dat) < 10)){
          outlier_regions <- as.numeric(names(which(table(regions_dat) < 10)))

          index_replacement <- vector()
          for(i in 1:length(outlier_regions)){
            index_replacement[i] <- replace.outliers(outlier_region = outlier_regions[i],
                                                     regions_dat = regions_dat,
                                                     clim_dat = clim_dat)
            regions_dat[regions_dat == outlier_regions[i]] <- as.numeric(index_replacement[i])
          }
        }

      }

      regions_dat = data.frame(Flyway.region = regions_dat) %>%
        dplyr::bind_cols(clim_dat[, c('Latitude', 'Longitude')]) %>%
        dplyr::ungroup()

      regions_dat$Site.group <- flyways[x]
      regions_dat
    })
  }

  #### Collect dataframes from each flyway and join together ####
  all_regions <- do.call(rbind, cluster_sites) %>%
    dplyr::mutate(Flyway.region.unique.id = paste0(Site.group, Flyway.region))

  joined_data <- all_regions %>%
    dplyr::right_join(Site.descriptors)

  return(joined_data)
}

