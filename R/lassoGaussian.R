#'Perform LASSO regularized gaussian regressions to model the process
#'that governs species' abundances
#'
#'\code{lassoGaussian} runs cross-validated regularized regressions for a single species
#'and returns important coefficients
#'
#'@importFrom magrittr %>%
#'
#'@param response A \code{vector} of count observations for the
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
#'the species' scaled abundance using \code{\link[glmnet]{cv.glmnet}}. These models
#'use coordinated gradient descent, applied to training sets of the data, to identify
#'regression parameters. These parameters are predicted on the remaining subset of the data
#'(the test set) to assess model fit. The process is repeated until a best-fitting model
#'is identified (minimising the loss function, which is cross-validated deviance in this case).
#'By replicating the process \code{n_reps} times, we account for uncertainty in the fold
#'generating process and can more confidently identify meaningful predictors (i.e. those that
#'are retained in at least \code{cutoff} proportion of \code{n_reps} models)
#'
#'@return \code{lassoGaussian} returns a single \code{vector} of coefficients for
#'predictors in \code{covariates}. \cr\cr
#'\code{lassoAbund_comm} binds these coefficient vectors into a
#'\code{dataframe} with rownames matching species names in \code{outcome_data}. It then
#'returns a \code{list} containing coefficients and scaling factors, which are used
#'in predictive functions
#'
#'@export
lassoGaussian = function(response, covariates, cutoff, n_reps){

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
                                x = as.matrix(covariates),
                                nfolds = 10, family = "gaussian",
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


#'\code{lassoAbund_comm} runs the lassoGaussian function across multiple species
#'in a community dataset
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param outcome_data A \code{dataframe} containing count observations
#'for species (each column representing a different species)
#'@param binary_data A \code{dataframe} containing binary presence-absence observations
#'for species (each column representing a different species)
#'@param outcome_indices A sequence of positive integers representing the column indices in
#'\code{outcome_data} that are to be modelled as poisson outcome variables (i.e. species
#'counts). Each one of these columns will be treated as a separate species whose
#'abundance is to be modelled using \code{lassoGaussian}. Default is to run models
#'for all columns in \code{outcome_data}
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@inheritParams lassoGaussian
#'@rdname lassoGaussian
#'
#'
#'@export
lassoAbund_comm = function(outcome_data, binary_data, outcome_indices,
                              covariates, cutoff, n_reps, n_cores){

  #### Specify default values ####
  if(missing(outcome_indices)){
    outcome_indices <- seq_len(ncol(outcome_data))
  }

  if(missing(n_reps)){
    n_reps <- 10
  }

  if(missing(cutoff)){
    cutoff <- 0.8
  }

  if(missing(n_cores)){
    n_cores <- detectCores() - 1
  }

  n_covariates <- ncol(covariates)

  #### Use sqrt mean transformation for Poisson variables ####
    square_root_mean = function(x) {sqrt(mean(x ^ 2))}
    poiss_sc_factors <- apply(outcome_data, 2, square_root_mean)
    outcome_data <- apply(outcome_data, 2,
                               function(x) x / square_root_mean(x))
    family <- 'gaussian'

  covariates <- cbind(outcome_data, covariates)

  #### For each outcome species, remove those species that infrequently (or never)
  # co-occur as possible predictors ####
  remain_vars <- lapply(seq_len(ncol(binary_data)), function(x){
    co.occurs <- colSums(binary_data[which(binary_data[,x] == 1),])
    co.occurs <- ifelse(co.occurs > (nrow(binary_data[which(binary_data[,x] == 1), ]) / 5),
                        TRUE, FALSE)
  })
  prepped_covariates <- MRFcov::prep_MRF_covariates(covariates,
                                                    n_nodes = ncol(outcome_data))

 cov_names <- colnames(prepped_covariates)[(max(outcome_indices) + 1):ncol(prepped_covariates)]

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
                          'covariates', 'n_reps', 'remain_vars',
                          'prepped_covariates', 'cutoff'),
                  envir = environment())

    #Export necessary functions to each cluster
    clusterExport(NULL, c('lassoGaussian'))

    #Export necessary libraries
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(glmnet))

    #### Perform lasso_gaussian function for each outcome variable ####
    cv_mods <- pbapply::pblapply(outcome_indices, function(x){

       #Make sure the focal species and infrequent co-occurring species are not present in covariates
      names_drop <- names(which(remain_vars[[x]] == FALSE))
       vars_drop <- lapply(seq_along(names_drop), function(j){
        vars_drop <- which(grepl(names_drop[j],
                                 colnames(prepped_covariates)) == TRUE)
      })
      all_vars_drop <- unique(unlist(vars_drop))

      if(is.null(all_vars_drop)){
        predictors <- prepped_covariates
      } else {
        predictors <- prepped_covariates[, -all_vars_drop]
      }

      predictors <- predictors[, -which(grepl(colnames(outcome_data)[x],
                                              colnames(predictors)) == T)]

      # Filter outcome and predictor datasets to only asses observations in which
      # the target species' abundance was greater than zero
      rows_keep <- which(outcome_data[, x] > 0)
      outcome_vector <- outcome_data[rows_keep, x]
      predictors <- predictors[rows_keep, ]

      result <- lassoGaussian(response = outcome_vector,
                              covariates = predictors,
                              cutoff = cutoff,
                              n_reps = n_reps)
    }, cl = cl)
    stopCluster(cl)
    names(cv_mods) <- colnames(outcome_data)[outcome_indices]

  } else {
    #### If parallel loading fails, use lapply instead (may crash!)
    cv_mods <- pbapply::pblapply(outcome_indices, function(x){

      # Remove species that do not co-occur frequently as possible predictors
      names_drop <- names(which(remain_vars[[x]] == FALSE))
      vars_drop <- lapply(seq_along(names_drop), function(j){
        vars_drop <- which(grepl(names_drop[j],
                                 colnames(prepped_covariates)) == TRUE)
      })
      all_vars_drop <- unique(unlist(vars_drop))

      if(is.null(all_vars_drop)){
        predictors <- prepped_covariates
      } else {
        predictors <- prepped_covariates[, -all_vars_drop]
      }

      predictors <- predictors[, -which(grepl(colnames(outcome_data)[x],
                                                  colnames(predictors)) == T)]

      # Filter outcome and predictor datasets to only assess observations in which
      # the target species' abundance was greater than zero
      rows_keep <- which(outcome_data[, x] > 0)
      outcome_vector <- outcome_data[rows_keep, x]
      predictors <- predictors[rows_keep, ]

      # Run the models
      result <- lassoGaussian(response = outcome_vector,
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

  # Add extra row with all columns in case some predictors are missing from all
  # models
  all_preds <- data.frame(t(rep(0, length(c('(Intercept)', colnames(prepped_covariates))))))
  colnames(all_preds) <- c('(Intercept)', colnames(prepped_covariates))

  results <- dplyr::bind_rows(all_preds, results)[-1, ]
  rownames(results) <- names(cv_mods)

  column_order <- c('(Intercept)', colnames(prepped_covariates))
  results <- results[, column_order]
  results[is.na(results)] <- 0

  #### Function to symmetrize corresponding coefficients ####
  symmetr <- function(coef_matrix, check_directs = FALSE, direct_upper = NULL,
                      symmetrise){
    coef_matrix_upper <- coef_matrix[upper.tri(coef_matrix)]
    coef_matrix.lower <- t(coef_matrix)[upper.tri(coef_matrix)]

    if(missing(symmetrise)){
      symmetrise <- 'mean'
    }

    if(symmetrise == 'min'){
      # If min, keep the coefficient with the smaller absolute value
      coef_matrix_upper_new <- ifelse(abs(coef_matrix_upper) < abs(coef_matrix.lower),
                                      coef_matrix_upper, coef_matrix.lower)
    }

    if(symmetrise == 'mean'){
      # If mean, take the mean of the two coefficients
      coef_matrix_upper_new <- (coef_matrix_upper + coef_matrix.lower) / 2
    }

    if(symmetrise == 'max'){
      # If max, keep the coefficient with the larger absolute value
      coef_matrix_upper_new <- ifelse(abs(coef_matrix_upper) > abs(coef_matrix.lower),
                                      coef_matrix_upper, coef_matrix.lower)
    }

    if(check_directs){
      # For indirect interactions, conditional relationships can only occur if
      # a direct interaction is found
      direct_upper <- direct_upper[upper.tri(direct_upper)]
      direct_upper[direct_upper > 0] <- 1
      coef_matrix_upper_new <- coef_matrix_upper_new * direct_upper
    }

    coef_matrix_sym <- matrix(0, max(outcome_indices), max(outcome_indices))
    coef_matrix_sym[upper.tri(coef_matrix_sym)] <- coef_matrix_upper_new
    coef_matrix_sym <- t(coef_matrix_sym)
    coef_matrix_sym[upper.tri(coef_matrix_sym)] <- coef_matrix_upper_new
    coef_matrix_sym
  }

  #### Create matrices of symmetric interaction coefficient estimates ####
  interaction_matrix <- results[, c(1 + outcome_indices)]

  #Symmetrize corresponding species interaction coefficients
  interaction_matrix_sym <- symmetr(interaction_matrix)
  rownames(interaction_matrix_sym) <- colnames(outcome_data)[outcome_indices]
  colnames(interaction_matrix_sym) <- colnames(outcome_data)[outcome_indices]

  #Replace unsymmetric direct interaction coefficients with the symmetric version
  results[, c(1 + outcome_indices)] <- interaction_matrix_sym

  #Create covariate coefficient matrices
    covariate_matrices <- lapply(seq_len(n_covariates), function(x){
      cov_matrix <- matrix(0, max(outcome_indices), max(outcome_indices))
      for(i in outcome_indices){
        cov_names_match <- grepl(paste('^', cov_names[x], '_', sep = ''),
                                 names(results))
        cov_matrix <- results[, cov_names_match]
        cov_matrix[is.na(cov_matrix)] <- 0
      }
      cov_matrix
    })

    #Symmetrize corresponding interaction coefficients
    indirect_coefs <- lapply(seq_along(covariate_matrices), function(x){
      matrix_to_sym <- covariate_matrices[[x]]
      sym_matrix <- symmetr(matrix_to_sym, check_directs = TRUE,
                            direct_upper = interaction_matrix_sym)
      rownames(sym_matrix) <- rownames(outcome_data)[outcome_indices]
      colnames(sym_matrix) <- colnames(outcome_data)[outcome_indices]
      list(sym_matrix)
    })
    names(indirect_coefs) <- cov_names[1:n_covariates]

    #Replace unsymmetric indirect interaction coefficients with symmetric versions
    covs_to_sym <- ncol(results) - (1 + max(outcome_indices) + n_covariates)
    covs_to_sym_end <- seq(max(outcome_indices), covs_to_sym,
                           by = max(outcome_indices)) + (1 + max(outcome_indices) + n_covariates)
    covs_to_sym_beg <- covs_to_sym_end - (max(outcome_indices) - 1)

    for(i in seq_len(n_covariates)){
      results[, covs_to_sym_beg[i] :
                     covs_to_sym_end[i]] <- data.frame(indirect_coefs[[i]])
    }

  return(list(direct_coefs = results,
              graph = interaction_matrix_sym,
              indirect_coefs = indirect_coefs,
              poiss_sc_factors = poiss_sc_factors,
              mod_type = "MRFcov"))
}

