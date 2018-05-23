#' Compiled data on avian species occurrences at NA BBS sites from
#' 2003 to 2009.
#'
#'
#' BBS data for each year was extracted using functions
#' in the 'pullBBS' package
#' (i.e. species.df <- pullBBS(year = 2003, country = 840, useCache = T)). The
#' entire workflow is shown in the \code{Appendix S1} script. Names were then
#' cleaned and filtered to only include species with names matching to
#' BirdTree.org using \code{validateSpecies.R}
#'
#' @format A \code{data.frame} object with 15863 rows, where each row lists the
#' binary occurrences of 303 avian species (stored in columns) across sites. Row order
#' matches the row orders in the \code{Site.descriptors} and \code{BBS.abundaneces} datafiles
#'
#' @docType data
#' @usage data("BBS.abundances")
#' @keywords datasets
"BBS.occurrences"
