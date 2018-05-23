#' Phylogenetic pairwise distance matrix for all avian species
#' included in \code{BBS.abundances} and \code{BBS.occurrences}
#'
#'
#' A pairwise phylogenetic distance matrices was created
#' using functions in the \code{ape} package. 100 phylogenetic trees were
#' downloaded from Birdtree.org and used to calculate mean pairwise distance. The pairwise
#' distances were then scaled to range \code{[0, 1]}. The
#' entire workflow is shown in the \code{Appendix S4} script.
#'
#' @format A \code{matrix} object with 303 rows and 303 columns, where each cell lists the
#' scaled pairwise phylogenetic distance between species in the corresponding row and column.
#' Row and column orders match the column orders in the \code{BBS.abundances}
#' and \code{BBS.occurrences} datafiles
#'
#' @docType data
#' @usage data("Bird.phy.distances")
#' @keywords datasets
"Bird.phy.distances"
