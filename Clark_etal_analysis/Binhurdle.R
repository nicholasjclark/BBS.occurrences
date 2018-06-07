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
  region_abund <- as.data.frame(region_abund[!row_has_na, ])

  #### Run LASSO occurrence and abundance models ####
  occurrence_mod <- lassoBinomial_comm(outcome_data = region_bin,
                                     count_data = region_abund,
                                     covariates = covariates,
                                     n_reps = 10, n_cores = 24)

  abundance_mod <- lassoAbund_comm(outcome_data = region_abund,
                                   binary_data = region_bin,
                                   covariates = covariates,
                                   n_reps = 10, n_cores = 24)
  list(occurrence_mod, abundance_mod)
})

names(region_mods) <- unique_regions

save(region_mods,
     all_regions,
     file = './BBS_results/Binhurdle.rda')
