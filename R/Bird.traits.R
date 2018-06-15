#' Information on avian species' migratory status (extracted from BirdLife International),
#' habitat and diet diversity (extracted from EltonTraits) and rarity
#'
#'
#' A \code{dataframe} containing variables describing the broad migratory status of each species in
#' \code{BBS.occurrences} (extracted from BirdLife International using web scraping) as well as
#' derived indices of species' diet and habitat diversities (calculated as a Shannon index
#' using data in the EltonTraits database) and species rarity (scaled variable representing
#' the proportional number of observations in \code{BBS.occurrences} in which a species
#' was observed). The
#' entire workflow is shown in the \code{Appendix S7} script.
#'
#' @format A \code{data.frame} object with 303 rows and five columns
#'
#' @docType data
#' @usage data("Bird.traits")
#' @keywords datasets
"Bird.traits"
