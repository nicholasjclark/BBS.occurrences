% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processNDVI.R
\name{processNDVI}
\alias{processNDVI}
\title{Extract mean summer NDVI values for buffers around lat/long points}
\usage{
processNDVI(filepath, points, buffer, northern_summer = T, n_cores)
}
\arguments{
\item{filepath}{A \code{character} string stating the path to the folder where raw
Gimmsfiles are stored (downloaded using \code{\link[gimms]{downloadGimms}})}

\item{points}{A \code{dataframe} of points with columns \code{Latitude} and \code{Longitude}}

\item{buffer}{Positive numeric value stating the length (in meters) of the desired buffer radius.
Mean values of the climate variable will be calculated for each point from within this buffer}

\item{northern_summer}{Logical. If TRUE, calculate NDVI for the northern hemisphere summer using
July observations. If FALSE, calculate NDVI for southern hemisphere summer using January observations.
Default is TRUE}

\item{n_cores}{Positive integer stating the number of processing cores to split the job across.
Default is \code{parallel::detect_cores() - 1}}
}
\description{
Extract mean summer NDVI values for buffers around lat/long points
}
