#' Functional pairwise distance matrix for all avian species
#' included in \code{BBS.abundances} and \code{BBS.occurrences}
#'
#'
#' A pairwise functional distance matrices was created
#' using functions in the \code{ade4} package. Species proportional use of
#' seven different foraging habitat categories and proportional use of ten diet categories
#' were accessed from the Eltontraits database and used to calculate a mixed variable coefficient
#' of distance using a generalised Gower's distance. The pairwise distances were then scaled to
#' range \code{[0, 1]}. The entire workflow is shown in the \code{Appendix S4} script.
#'
#' @format A \code{matrix} object with 303 rows and 303 columns, where each cell lists the
#' scaled pairwise functional distance between species in the corresponding row and column.
#' Row and column orders match the column orders in the \code{BBS.abundances}
#' and \code{BBS.occurrences} datafiles
#'
#' @docType data
#' @usage data("Bird.ecol.distances")
#' @keywords datasets
"Bird.ecol.distances"
