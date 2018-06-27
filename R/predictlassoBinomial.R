#'Predict binomial outcomes from a fitted lassoBinomial_comm model
#'
#'This function calculates linear predictors for species' binary observations
#'using equations from a \code{\link{lassoBinomial_comm}} model.
#'
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param outcome_data \code{dataframe}. The input data to be predicted, where the
#'columns represent binary occurrences of the species.
#'Colnames from this sample dataset must exactly match the colnames in the dataset that
#'was used to fit the \code{lassoBinomial_comm}
#'@param count_data A \code{dataframe} containing count observations
#'for species (each column representing a different species)
#'@param covariates A \code{matrix} of covariates, with \code{nrow(covariates) == nrow(response)}
#'@param lassoBinomial A fitted \code{\link{lassoBinomial_comm}} model object
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'@return A \code{matrix} containing predictions for each observation in \code{outcome_data}.
#'
#'@export
predictlassoBinomial = function(outcome_data, lassoBinomial, count_data,
                                covariates, n_cores){

  n_nodes <- nrow(lassoBinomial$coefficients)

  # Function to back-transform logistic coefficients using inverse logit
  inverse_logit = function(x){
    exp(x) / (1 + exp(x))
  }

  # Extract relevant facets of data
  outcome_data <- data.frame(outcome_data)
  n_obs <- nrow(outcome_data)
  node_names <- colnames(outcome_data[, 1:n_nodes])

  # Scale prediction count variables using the same
  # scale factors as in the original lassoBinomial model
    for(i in seq_len(n_nodes)){
      count_data[, i] <- count_data[, i] / lassoBinomial$poiss_sc_factors[[i]]
    }

  #### If n_cores > 1, check parallel library loading ####
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a library on each cluster ####
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

  # Calculate linear predictions using the coefficients element from the model
  covariates <- cbind(count_data, covariates)

  if(parallel_compliant){

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('covariates', 'n_nodes', 'lassoBinomial'),
                  envir = environment())

    predictions <- do.call(cbind, pbapply::pblapply(seq_len(n_nodes), function(i){
      apply(covariates, 1, function(j) sum(j %*% t(lassoBinomial$coefficients[i, -1])) +
              lassoBinomial$intercepts[i])
    }, cl = cl))
    stopCluster(cl)

    } else {

    predictions <- do.call(cbind, lapply(seq_len(n_nodes), function(i){
      apply(covariates, 1, function(j) sum(j %*% t(lassoBinomial$coefficients[i, -1])) +
              lassoBinomial$intercepts[i])
    }))

  }
    colnames(predictions) <- node_names

    # Convert linear predictions to probability scale
    predictions <- inverse_logit(predictions)

    binary_predictions <- ifelse(predictions >= 0.5, 1, 0)

    return(list(Probability_predictions = round(predictions, 4),
                Binary_predictions = binary_predictions))
}

#' Replicate predictlassoBinomial across folds of the data
#'
#' \code{lassoBinomial_metrics} splits the data into ten folds and tests how well
#' the \code{lassoBinomial} model predicts binary observations in each fold.
#' \cr
#' \cr
#' Matrices of predictive metrics are returned.
#'
#' @inheritParams predictlassoBinomial
#' @rdname predictlassoBinomial
#'
#' @export
lassoBinomial_metrics = function(outcome_data, lassoBinomial, count_data,
                                 covariates, n_cores){

  if(missing(n_cores)){
    n_cores <- parallel::detectCores() - 1
  }

  folds <- caret::createFolds(rownames(outcome_data), 10)
  n_folds <- 10

   all_predictions <- predictlassoBinomial(outcome_data, lassoBinomial, count_data,
                                           covariates, n_cores = n_cores)

   n_nodes <- nrow(lassoBinomial$coefficients)

  cv_predictions <- lapply(seq_len(n_folds), function(k){
    test_data <- outcome_data[folds[[k]], ]
    predictions <- all_predictions[[2]][folds[[k]], ]

    #Calculate positive and negative predictive values
    true_pos <- (predictions == test_data[, c(1:n_nodes)])[test_data[, c(1:n_nodes)] == 1]
    false_pos <- (predictions == 1)[predictions != test_data[, c(1:n_nodes)]]
    true_neg <- (predictions == test_data[, c(1:n_nodes)])[test_data[, c(1:n_nodes)] == 0]
    false_neg <- (predictions == 0)[predictions != test_data[, c(1:n_nodes)]]

    #Calculate diagnostic predictive values
    pos_pred <- sum(true_pos, na.rm = TRUE) /
      (sum(true_pos, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
    neg_pred <- sum(true_neg, na.rm = TRUE) /
      (sum(true_neg, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
    sensitivity <- sum(true_pos, na.rm = TRUE) /
      (sum(true_pos, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
    specificity <- sum(true_neg, na.rm = TRUE) /
      (sum(true_neg, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
    tot_pred <- (sum(true_pos, na.rm = TRUE) + sum(true_neg, na.rm = TRUE)) /
      (length(as.matrix(test_data[, c(1:n_nodes)])))

    list(mean_pos_pred = mean(pos_pred, na.rm = TRUE),
         mean_neg_pred = mean(neg_pred, na.rm = TRUE),
         mean_tot_pred = mean(tot_pred, na.rm = TRUE),
         mean_sensitivity = mean(sensitivity, na.rm = TRUE),
         mean_specificity = mean(specificity, na.rm = TRUE))
  })

  pred_dat <- purrr::map_df(cv_predictions, magrittr::extract,
                                 c('mean_pos_pred','mean_tot_pred',
                                   'mean_sensitivity',
                                   'mean_specificity'))
  return(Metrics = pred_dat)
}

