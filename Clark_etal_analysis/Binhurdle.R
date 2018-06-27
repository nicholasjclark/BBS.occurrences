library(dplyr)
library(BBS.occurrences)
data("BBS.occurrences")
data("BBS.abundances")
data("Site.descriptors")

# Make sure that MRFcov development version is installed
if(!require(MRFcov)){
  devtools::install_github('nicholasjclark/MRFcov@development')
}

# Convert sites within flyways into biogeographical regions
all_regions <- createRegions(n_cores = 4)

# Create unique site IDs for spatio-temporal modelling
all_regions %>%
  dplyr::select(Latitude, Longitude) %>%
  dplyr::distinct() %>%
  dplyr::mutate(ID = as.character(dplyr::id(.))) %>% #ID must be a character
  dplyr::ungroup() %>%
  dplyr::left_join(all_regions) -> all_regions

# Impute missing landcover and ndvi values for each site
all_regions %>%

  # First, impute using the mean for each site
  dplyr::group_by(ID) %>%
  dplyr::mutate_if(is.numeric, dplyr::funs(ifelse(is.na(.), mean(., na.rm=T), .))) %>%
  dplyr::ungroup() %>%

  # If still NA, impute using mean from the respective county
  dplyr::group_by(County) %>%
  dplyr::mutate_if(is.numeric, dplyr::funs(ifelse(is.na(.), mean(., na.rm=T), .))) %>%
  dplyr::ungroup() %>%

  # If still NA, impute using mean from respective region
  dplyr::group_by(Flyway.region.unique.id) %>%
  dplyr::mutate_if(is.numeric, dplyr::funs(ifelse(is.na(.), mean(., na.rm=T), .))) %>%
  dplyr::ungroup() -> all_regions

# Test models on a single flyway (Mississippi)
# Gather names of flyways
flyways <- as.character(unique(all_regions$Site.group))

# Fit the hurdle model (5 reps) for two flyways to ensure it works
Miss_mods <- hurdleModel(flyway = flyways[1], n_cores = 24, n_reps = 5)
length(Miss_mods$removed_covs)

Atl_mods <- hurdleModel(flyway = flyways[4], n_cores = 24, n_reps = 5)
length(Atl_mods$removed_covs)

save(Miss_mods,
     Atl_mods,
     all_regions,
     file = './BBS_results/Binhurdle.rda')
