#'Process mean monthly PRISM climate values into binned variables
#'
#'
#'@importFrom magrittr %>%
#'
#'@param monthly_climate_data A \code{dataframe} of PRISM data for
#'points. Individual datasets should be created using \code{\link{meanPRISM}} and joined
#'together (e.g. with \code{\link[dplyr]{left_join}}), resulting in a single \code{dataframe}
#' with columns \code{Latitude}, \code{Longitude}, \code{Year}, \code{Month},
#'\code{Max.temp}, \code{Min.temp}, \code{Mean.temp} and \code{Precipitation}
#'
#'@export
processPRISM = function(monthly_climate_data){

  #### Calculate temp and precip variables for the three months prior to sampling ####
  threemonths_prior <- data.frame(Month = 1:12, Quarter = c(0,0,1,1,1,0,0,0,0,0,0,0))

  bioclim_data_3months = monthly_climate_data %>%
    dplyr::left_join(threemonths_prior, by = 'Month') %>%
    dplyr::group_by(Latitude, Longitude, Year, Quarter) %>%
    dplyr::mutate(Tot.Precip.Spr = sum(Precipitation),
                  Mean.Temp.Spr = mean(Mean.temp, na.rm = T),
                  Max.Temp.Spr = max(Max.temp, na.rm = T),
                  Min.Temp.Spr = min(Min.temp, na.rm = T),
                  Mean.Min.Temp.Spr = mean(Min.temp, na.rm = T),
                  Mean.Max.Temp.Spr = mean(Max.temp, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Quarter == 1) %>%
    dplyr::select(-Quarter, -Month, -Min.temp, -Max.temp, -Mean.temp, -Precipitation)

    #Infinite values ned to be set to NA
    bioclim_data_3months <- do.call(data.frame,
                                  lapply(bioclim_data_3months,
                                         function(x) replace(x, is.infinite(x), NA)))


    #### Now calculate for the year prior to sampling ####
    #First, we offset the year by 6 months so that the window for calculating bioclim variables
    #will be July 1 - June 30. See https://github.com/weecology/bbs-forecasting/issues/114
    monthly_climate_data$Year <- with(monthly_climate_data, ifelse(Month %in% 7:12, Year+1, Year))

    bioclim_data_1year = monthly_climate_data %>%
      dplyr::group_by(Latitude, Longitude, Year) %>%
      dplyr::mutate(Tot.Precip.Yr = sum(Precipitation, na.rm = T),
                    Precip.Seas.Yr = sd(Precipitation, na.rm = T),
                    Mean.Temp.Yr = mean(Mean.temp, na.rm = T),
                    Mean.Min.Temp.Yr = mean(Min.temp, na.rm = T),
                    Mean.Max.Temp.Yr = mean(Max.temp, na.rm = T),
                    Temp.Seas.Yr = sd(Mean.temp, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-Month, -Min.temp, -Max.temp, -Precipitation, -Mean.temp)

 bioclim_data = bioclim_data_1year %>%
   dplyr::left_join(bioclim_data_3months) %>%
   dplyr::distinct()

return(bioclim_data)
}
