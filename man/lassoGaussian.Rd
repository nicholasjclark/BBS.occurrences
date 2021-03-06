% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lassoGaussian.R
\name{lassoGaussian}
\alias{lassoGaussian}
\alias{lassoAbund_comm}
\title{Perform LASSO regularized gaussian regressions to model the process
that governs species' abundances}
\usage{
lassoGaussian(response, covariates, cutoff, n_reps)

lassoAbund_comm(outcome_data, binary_data, outcome_indices, covariates,
  cutoff, n_reps, n_cores)
}
\arguments{
\item{response}{A \code{vector} of count observations for the
focal species}

\item{covariates}{A \code{matrix} of covariates, with \code{nrow(covariates) == nrow(response)}}

\item{cutoff}{Positive numeric value representing the proportion of models in which
a predictor must be retained in order to be treated as meaningful. If the predictor is
retained in fewer than \code{cutoff} proportion of \code{\link[glmnet]{cv.glmnet}}
regularized models, its coefficient is forced to be zero. If the predictor is retained
in at least \code{cutoff} proportion of models, its mean coefficient is returned.
Default is \code{0.80}}

\item{n_reps}{Positive \code{integer} representing the number of times to repeat
10-fold \code{\link[glmnet]{cv.glmnet}} regularized regressions (default is \code{10})}

\item{outcome_data}{A \code{dataframe} containing count observations
for species (each column representing a different species)}

\item{binary_data}{A \code{dataframe} containing binary presence-absence observations
for species (each column representing a different species)}

\item{outcome_indices}{A sequence of positive integers representing the column indices in
\code{outcome_data} that are to be modelled as poisson outcome variables (i.e. species
counts). Each one of these columns will be treated as a separate species whose
abundance is to be modelled using \code{lassoGaussian}. Default is to run models
for all columns in \code{outcome_data}}

\item{n_cores}{Positive integer stating the number of processing cores to split the job across.
Default is \code{parallel::detect_cores() - 1}}
}
\value{
\code{lassoGaussian} returns a single \code{vector} of coefficients for
predictors in \code{covariates}. \cr\cr
\code{lassoAbund_comm} binds these coefficient vectors into a
\code{dataframe} with rownames matching species names in \code{outcome_data}. It then
returns a \code{list} containing coefficients and scaling factors, which are used
in predictive functions
}
\description{
\code{lassoGaussian} runs cross-validated regularized regressions for a single species
and returns important coefficients
}
\details{
Regularized regressions are performed to identify meaningful predictors of
the species' scaled abundance using \code{\link[glmnet]{cv.glmnet}}. These models
use coordinated gradient descent, applied to training sets of the data, to identify
regression parameters. These parameters are predicted on the remaining subset of the data
(the test set) to assess model fit. The process is repeated until a best-fitting model
is identified (minimising the loss function, which is cross-validated deviance in this case).
By replicating the process \code{n_reps} times, we account for uncertainty in the fold
generating process and can more confidently identify meaningful predictors (i.e. those that
are retained in at least \code{cutoff} proportion of \code{n_reps} models)
}
\seealso{
\code{\link[glmnet]{cv.glmnet}}
}
