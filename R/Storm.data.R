#' Descriptions and locations of damaging storms in the USA
#'
#' A dataset on storm events was accessed from NOAA and processed using regex
#' to categorise storms into categories. The entire workflow is shown in the \code{Appendix 6} script.
#' @format A \code{data.frame} containing 15 variables describing
#' location, severity and timing of the top third most damaging weather events recorded in the
#' NOAA storms database. Variables
#' include property damage, crop damage and total damage (all in US dollars), Year, Month,
#' Storm.type and geographrical coordinates
#'
#' @docType data
#' @usage data(Storm.data)
#' @keywords datasets
#'
"Storm.data"
