#'Perform LASSO regularized logistic regressions to model the process
#'that generates zeros and non-zeros in species abundance data
#'
#'\code{lassoBinomial} runs cross-validated regularized regressions for a single species
#'and returns important coefficients
#'
#'@importFrom magrittr %>%
#'
#'@param response A \code{vector} of binary presence-absence observations for the
#'focal species
#'@param covariates A \code{matrix} of covariates, with \code{nrow(covariates) == nrow(response)}
#'@param cutoff Positive numeric value representing the proportion of models in which
#'a predictor must be retained in order to be treated as meaningful. If the predictor is
#'retained in fewer than \code{cutoff} proportion of \code{\link[glmnet]{cv.glmnet}}
#'regularized models, its coefficient is forced to be zero. If the predictor is retained
#'in at least \code{cutoff} proportion of models, its mean coefficient is returned.
#'Default is \code{0.80}
#'@param n_reps Positive \code{integer} representing the number of times to repeat
#'10-fold \code{\link[glmnet]{cv.glmnet}} regularized regressions (default is \code{10})
#'
#'@seealso \code{\link[glmnet]{cv.glmnet}}
#'
#'@details Regularized regressions are performed to identify meaningful predictors of
#'the species' occurrence probability using \code{\link[glmnet]{cv.glmnet}}. These models
#'use coordinated gradient descent, applied to training sets of the data, to identify
#'regression parameters. These parameters are predicted on the remaining subset of the data
#'(the test set) to assess model fit. The process is repeated until a best-fitting model
#'is identified (minimising the loss function, which is cross-validated deviance in this case).
#'By replicating the process \code{n_reps} times, we account for uncertainty in the fold
#'generating process and can more confidently identify meaningful predictors (i.e. those that
#'are retained in at least \code{cutoff} proportion of \code{n_reps} models)
#'
#'@return \code{lassoBinomial} returns a single \code{vector} of coefficients for
#'predictors in \code{covariates}. \cr\cr
#'\code{lassoBinomial_comm} binds these vectors into a
#'\code{dataframe} with rownames matching species names in \code{outcome_data}
#'
#'@export
lassoBinomial = function(response, covariates, cutoff, n_reps){

  #### Specify default values ####
  if(missing(n_reps)){
    n_reps <- 10
  }

  if(missing(cutoff)){
    cutoff <- 0.8
  }

  #### Perform cross-validated LASSO regression using the binary response ####
  #Note, this function is replicated n_reps times to account for uncertainty in
  #fold generation
  cv.lasso <- lapply(seq_len(n_reps), function(x){
    cv.mod <- glmnet::cv.glmnet(y = response,
                                x = covariates,
                                nfolds = 10, family = "binomial",
                                weights = rep(1, nrow(covariates)),
                                intercept = TRUE, standardize = TRUE)
    coef(cv.mod)

    # Gather coefficients from the best-fitting model into a list
    coefs <- as.vector(t(as.matrix(coef(cv.mod, s = 'lambda.min'))))
    list(coef = coefs,
         parameter = dimnames(t(as.matrix(coef(cv.mod,
                                               s = 'lambda.min'))))[[2]])
  })

  #### Keep important coefficients (retained by more than 80% of models) ####
  #Regularize remaining coefficients to zero
  coefs.keep <- t(as.matrix(do.call(rbind, lapply(cv.lasso,
                                                  as.data.frame)) %>%
                              dplyr::group_by(parameter) %>%
                              dplyr::mutate(prop.retained = length(which(coef != 0)) / n_reps,
                                     mean.coef = mean(coef)) %>%
                              dplyr::rowwise() %>%
                              dplyr::mutate(mean.coef = ifelse(prop.retained >= cutoff,
                                                        mean.coef, 0)) %>%
                              dplyr::ungroup() %>%
                              dplyr::select(parameter, mean.coef) %>%
                              dplyr::mutate_at(dplyr::vars(parameter), as.character) %>%
                              dplyr::distinct()))

  #Convert coefficients to a named numeric vector to save memory
  coef.names <- coefs.keep[1, ]
  coefs.keep <- as.numeric(coefs.keep[-1, ])
  names(coefs.keep) <- coef.names

  return(coefs.keep)
}



#'\code{lassoBinomial_comm} runs the lassoBinomial function across multiple species
#'in a community dataset
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param outcome_data A \code{dataframe} containing binary presence-absence observations
#'for species (each column representing a different species)
#'@param outcome_indices A sequence of positive integers representing the column indices in
#'\code{outcome_data} that are to be modelled as binary outcome variables (i.e. species
#'occurrence observations). Each one of these columns will be treated as a separate species whose
#'occurrence probability is to be modelled using \code{lassoBinomial}. Default is to run models
#'for all columns in \code{outcome_data}
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@inheritParams lassoBinomial
#'@rdname lassoBinomial
#'
#'
#'@export
lassoBinomial_comm = function(outcome_data, outcome_indices,
                              covariates, cutoff, n_reps, n_cores){

  #### Specify default values ####
  if(missing(outcome_indices)){
    outcome_indices <- seq_len(ncol(outcome_data))
  }

  if(missing(cutoff)){
    cutoff <- 0.8
  }

  if(missing(n_cores)){
    n_cores <- detectCores() - 1
  }

  #### If n_cores > 1, check parallel library loading ####
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(glmnet)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/glmnet$", "", system.file(package = "glmnet"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(glmnet)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(glmnet)), silent = TRUE)

        if(class(test_load3) == "try-error"){

          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)

        } else {
          parallel_compliant <- TRUE
        }

      } else {
        parallel_compliant <- TRUE
      }

    } else {
      parallel_compliant <- TRUE
    }
  } else {
    #If n_cores = 1, set parallel_compliant to FALSE
    parallel_compliant <- FALSE
    warning('Parallel loading failed, calculations may crash!')
  }

  if(parallel_compliant){

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('outcome_indices', 'outcome_data',
                          'covariates', 'n_reps'),
                  envir = environment())

    #Export necessary functions to each cluster
    clusterExport(NULL, c('lassoBinomial'))

    #Export necessary libraries
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(glmnet))

    #### Perform lasso_binomial function for each outcome variable ####
    cv_mods <- pbapply::pblapply(outcome_indices, function(x){

      #Make sure the focal species is not present in the covariates data
      predictors <- covariates[, -which(grepl(colnames(outcome_data)[x], colnames(covariates)) == T)]
      result <- lassoBinomial(response = outcome_data[, x],
                              covariates = predictors,
                              cutoff = cutoff,
                              n_reps = n_reps)
    })
    stopCluster(cl)
    names(cv_mods) <- colnames(outcome_data)[outcome_indices]

  } else {
    #### If parallel loading fails, use lapply instead (may crash!)
    cv_mods <- pbapply::pblapply(outcome_indices, function(x){
      predictors <- covariates[, -which(grepl(colnames(outcome_data)[x], colnames(covariates)) == T)]
      result <- lassoBinomial(response = outcome_data[, x],
                              covariates = predictors,
                              cutoff = cutoff,
                              n_reps = n_reps)
    })
    names(cv_mods) <- colnames(outcome_data)[outcome_indices]
  }

  #### Bind resulting model coefficient vectors into a dataframe and return ####
  results <- dplyr::bind_rows(lapply(cv_mods, function(x){
    as.data.frame(t(x))
    }))
  results[is.na(results)] <- 0
  rownames(results) <- names(cv_mods)

  #Make sure column order of results matches that of covariates (+ the intercept column)
  reorder_ids <- match(c('(Intercept)', colnames(covariates)),
                         colnames(results))
  results <- results[ , reorder_ids]

  return(results)
}
