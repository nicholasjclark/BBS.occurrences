% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processAltitudes.R
\name{processAltitudes}
\alias{processAltitudes}
\title{Extract mean altitude within a 40km buffer for lat and long coordinates}
\usage{
processAltitudes(points, n_cores, delete_files)
}
\arguments{
\item{points}{A \code{dataframe} containing the coordinates, with colnames 'Latitude' and
'Longitude'}

\item{n_cores}{Positive integer stating the number of processing cores to split the job across.
Default is \code{parallel::detect_cores() - 1}}

\item{delete_files}{Logical describing whether or not to delete raw elevation map once the job
is complete. Default is \code{FALSE}}
}
\description{
This function calculates mean elevation with 40km buffers of specified points. Note,
the function will temporarily save the global digital elevation map in a new folder
in the working directory, so please ensure there is sufficient memory (around 37MB) available
}
