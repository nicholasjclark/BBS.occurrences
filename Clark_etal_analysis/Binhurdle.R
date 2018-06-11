library(dplyr)
library(BBS.occurrences)
data("BBS.occurrences")
data("BBS.abundances")
data("Site.descriptors")

# Convert sites within flyways into biogeographical regions
all_regions <- createRegions(n_cores = 4)

# Test models on a single flyway (Mississippi)
# Gather names of flyways
flyways <- as.character(unique(all_regions$Site.group))

# For later lapply calls
x <- 1

# Split into regions
unique_regions = as.vector(all_regions %>%
  dplyr::filter(Site.group == flyways[x]) %>%
  dplyr::select(Flyway.region.unique.id) %>%
  dplyr::distinct())[,1]

#### Run separate models for each region within a flyway ####
region_mods <- lapply(seq_len(length(unique_regions)), function(j){

  region_rows <- which(all_regions$Site.group == flyways[x] &
                         all_regions$Flyway.region.unique.id == unique_regions[j])

  # Extract abundance and binary occurrence data for the specified region
  region_abund <- BBS.abundances[region_rows, ]
  region_bin <- BBS.occurrences[region_rows, ]

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
  covariates <- all_regions[region_rows, c(2, 31:48)]

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

  # Scale all continuous covariates so magnitudes of effects can be compared
  covariates = data.frame(covariates %>%
                            dplyr::mutate_if(is.numeric, funs(as.vector(scale(.)))))

  #### Run LASSO occurrence and abundance models ####
  occurrence_mod <- lassoBinomial_comm(outcome_data = region_bin_outcome,
                                     count_data = region_abund_predict,
                                     covariates = covariates,
                                     n_reps = 25, n_cores = 24)

  abundance_mod <- lassoAbund_comm(outcome_data = region_abund,
                                   binary_data = region_bin,
                                   covariates = covariates,
                                   n_reps = 25, n_cores = 24)

  #### Calculate predictive metrics for the binomial model ####
  occurrence_mets <- lassoBinomial_metrics(outcome_data = region_bin_outcome,
                                           count_data = region_abund_predict,
                                           covariates = covariates,
                                           lassoBinomial = occurrence_mod,
                                           n_cores = 24)

  #### Predict network centrality and assess model fit from abundance model ####
  abundance_cent <- MRFcov::predict_MRFnetworks(data = cbind(region_abund, covariates),
                                        MRF_mod = abundance_mod, metric = "eigencentrality",
                                        n_cores = 24)

  abundance_cent <- cbind(abundance_cent, covariates)

  abund_predictions <- MRFcov::predict_MRF(data = cbind(region_abund, covariates),
                                           MRF_mod = abundance_mod, n_cores = 24)

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
       abundance_cent = abundance_cent,
       abund_Rsquared = quantile(abund_metrics$Rsquared,
                                 probs = c(0.025, 0.5, 0.975)),
       removed_covs = removed_covs)
})

names(region_mods) <- unique_regions

save(region_mods,
     all_regions,
     file = './BBS_results/Binhurdle.rda')
