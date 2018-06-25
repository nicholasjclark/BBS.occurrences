#' Run linear mixed models to identify predictors of species' eigencentralities
#' across sites
#'
#'@importFrom magrittr %>%
#'
#'@param flyway A \code{character} string representing the name of the specified flyway to model
#'@param n_reps Positive \code{integer} representing the number of times to repeat
#'10-fold \code{\link[glmnet]{cv.glmnet}} regularized regressions (default is \code{10})
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@details This function runs the hurdle model by first estimating predictors of species' occurrence
#'probability using \code{lassoBinomial} and then fitting a conditional random field to species'
#'scaled abundances using \code{lassoGaussian}. Each of these models is run \code{n_reps} times
#'to allow for uncertainty in coefficient estimates. Predictions are then used to estimate model fits and
#'predict community network \emph{B}'os and species' eigencentrality in each site
#'
#'@export
#'
hurdleModel = function(flyway, n_reps, n_cores){

  if(missing(n_reps)){
    n_reps <- 10
  }

# Split into regions
unique_regions = as.vector(all_regions %>%
                             dplyr::filter(Site.group == flyway) %>%
                             dplyr::select(Flyway.region.unique.id) %>%
                             dplyr::distinct())[,1]

#### Run separate models for each region within a flyway ####
region_mods <- lapply(seq_len(length(unique_regions)), function(j){

  region_rows <- which(all_regions$Site.group == flyway &
                         all_regions$Flyway.region.unique.id == unique_regions[j])

  # Extract abundance and binary occurrence data for the specified region
  region_abund <- BBS.abundances[region_rows, ]
  region_bin <- BBS.occurrences[region_rows, ]

  # Extract coordinates and altitude to ensure their inclusion in spatial network models
  coords <- all_regions %>%
    dplyr::select(ID, Latitude, Longitude, Year, Altitude)
  coords <- coords[region_rows, ]

  # Remove species occurring in fewer than 15% of observations
  low_occur_cols <- which((colSums(region_bin) / nrow(region_bin)) < 0.15)
  region_bin <- region_bin[, -low_occur_cols]
  region_abund <- region_abund[, -low_occur_cols]

  # Species that are too common cannot be assessed in occurrence models
  high_occur_cols <- which((colSums(region_bin) / nrow(region_bin)) > 0.9)
  region_bin_outcome <- region_bin[, -high_occur_cols]
  region_abund_predict <- region_abund[, -high_occur_cols]

  # Remove un-needed rows from covariates data and extract climate / landcover covariates
  # note, year is not included as temporal autocorrelation will be
  # difficult to detect without at least 20 - 25 years worth of data
  # see (Teller et al., 2016) at
  # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12486)
  covs_test <- c("Latitude", "Prop.forest",
                 "Prop.shrubland",
                 "Prop.grassland",
                 "Prop.wetland",
                 "Prop.anthro.urban",
                 "Prop.anthro.cropland",
                 "Tot.precip.yr",
                 "Precip.seas.yr",
                 "Mean.temp.yr",
                 "Mean.min.temp.yr",
                 "Mean.max.temp.yr",
                 "Temp.seas.yr",
                 "Tot.precip.spr",
                 "Mean.temp.spr",
                 "Max.temp.spr",
                 "Min.temp.spr",
                 "Mean.min.temp.spr",
                 "Mean.max.temp.spr",
                 "Altitude", "NDVI")

  covariates <- all_regions[region_rows, covs_test]

  # Remove predictors with strong collinearity (based on pearson correlations)
  covs_remove <- caret::findCorrelation(cor(covariates,
                                            use = "na.or.complete"),
                                        cutoff = .95, exact = T)
  removed_covs <- names(covariates)[covs_remove]
  covariates <- covariates[, -covs_remove]

  # A small number of year * site observations do not have landcover values
  # here we remove these observations from all datasets
  row_has_na <- apply(covariates, 1, function(x){any(is.na(x))})
  covariates <- as.data.frame(covariates[!row_has_na, ])
  region_bin <- as.data.frame(region_bin[!row_has_na, ])
  region_bin_outcome <- as.data.frame(region_bin_outcome[!row_has_na, ])
  region_abund <- as.data.frame(region_abund[!row_has_na, ])
  region_abund_predict <- as.data.frame(region_abund_predict[!row_has_na, ])
  coords <- coords[!row_has_na, ]

  # Bin very large counts at the 99th percentile (these likely add no extra insight, only noise)
  top_count <- ceiling(quantile(as.vector(as.matrix(region_abund)), probs = 0.99))

  region_abund[region_abund > top_count] <- top_count
  region_abund_predict[region_abund_predict > top_count] <- top_count

  # Scale all continuous covariates so magnitudes of effects can be compared
  covariates <- data.frame(covariates %>%
                            dplyr::mutate_if(is.numeric, funs(as.vector(scale(.)))))

  #### Run LASSO occurrence and abundance models ####
  occurrence_mod <- lassoBinomial_comm(outcome_data = region_bin_outcome,
                                       count_data = region_abund_predict,
                                       covariates = covariates,
                                       n_reps = n_reps, n_cores = n_cores)

  abundance_mod <- lassoAbund_comm(outcome_data = region_abund,
                                   binary_data = region_bin,
                                   covariates = covariates,
                                   n_reps = n_reps, n_cores = n_cores)

  #### Calculate predictive metrics for the binomial model ####
  occurrence_mets <- lassoBinomial_metrics(outcome_data = region_bin_outcome,
                                           count_data = region_abund_predict,
                                           covariates = covariates,
                                           lassoBinomial = occurrence_mod,
                                           n_cores = n_cores)

  #### Predict network centrality and local network beta diversity (B'os) ####
  abundance_cent <- MRFcov::predict_MRFnetworks(data = cbind(region_abund, covariates),
                                                MRF_mod = abundance_mod, metric = "eigencentrality",
                                                cutoff = 2,
                                                omit_zeros = TRUE,
                                                n_cores = n_cores)

  network_metrics <- cbind(abundance_cent, covariates)

  abundance_adj <- MRFcov::predict_MRFnetworks(data = cbind(region_abund, covariates),
                                               MRF_mod = abundance_mod, metric = 'adjacency',
                                               cutoff = 2,
                                               omit_zeros = TRUE,
                                               n_cores = n_cores)

  beta_os_primes <- networkBetaOS(abundance_adj, n_cores = n_cores)

  network_metrics$beta_os <- beta_os_primes

  rm(abundance_cent, abundance_adj)

  #### Assess model fit from abundance model ####
  abund_predictions <- MRFcov::predict_MRF(data = cbind(region_abund, covariates),
                                           MRF_mod = abundance_mod, n_cores = n_cores)

  abund_metrics <- MRFcov::cv_MRF_diag(data = cbind(region_abund, covariates),
                                       n_nodes = nrow(abundance_mod$direct_coefs),
                                       n_folds = 10, n_cores = 1, family = 'poisson',
                                       compare_null = FALSE, plot = FALSE,
                                       cached_model = list(mrf = abundance_mod),
                                       cached_predictions = list(predictions = abund_predictions),
                                       sample_seed = 1)

  list(occurrence_mod = occurrence_mod,
       occurrence_spec = quantile(occurrence_mets$mean_specificity,
                                  probs = c(0.025, 0.5, 0.975)),
       occurrence_sens = quantile(occurrence_mets$mean_sensitivity,
                                  probs = c(0.025, 0.5, 0.975)),
       abundance_mod = abundance_mod,
       network_metrics = network_metrics,
       abund_Rsquared = quantile(abund_metrics$Rsquared,
                                 probs = c(0.025, 0.5, 0.975)),
       coordinates = coords,
       removed_covs = removed_covs,
       top_count = top_count)
})

names(region_mods) <- unique_regions
return(region_mods)
}
