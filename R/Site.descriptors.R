#' Site - level covariates for the NA BBS observations
#'
#' A dataset on site-level predictors for observations in \code{BBS.occurrences} and
#' \code{BBS.abundances} was compiled from multiple sources. The workflows for generating this dataset
#' are shown in \code{Appendices 1, 2, 3, 5}
#'
#' @format A \code{data.frame} containing 49 variables describing
#' covariates and other descriptors for the individual BBS observations. Variables
#' include standardised information gathered by BBS (i.e. route IDs, temperatures, times,
#' etc...) as well as GPS coordinates and percent cover for a range of
#' landcover covariates. Landcover variables were accessed at 40km resolution, downloaded
#' from MODIS Landcover Type 1, product MCD12Q1. See scripts at
#' https://github.com/nicholasjclark/LandcoverMODIS for more details on MODIS downloads. PRISM
#' climate variables were downloaded at 40km resolution using
#' \code{\link[prism]{get_prism_monthlys}} and coerced into lagged variables representing
#' climate conditions over the previous Spring (March - May) and over the previous year
#' (July - June) before sampling occurred. Altitude data was downloaded at 40km resolution using
#' the Aster Global Digital Elevation map. NDVI data was downloaded at 40km resolution using
#' \code{\link[gimms]{downloadGimms}}
#'
#' @docType data
#' @usage data(Site.descriptors)
#' @keywords datasets
#'
"Site.descriptors"
