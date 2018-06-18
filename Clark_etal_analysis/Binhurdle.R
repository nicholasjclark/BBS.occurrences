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

# Test models on a single flyway (Mississippi)
# Gather names of flyways
flyways <- as.character(unique(all_regions$Site.group))

# Fit the hurdle model for the first flyway
region_mods <- hurdleModel(flyway = flyway[1], n_cores = 24)

# Estimate predictors of network metrics
network_mods <- networkModel(mods_list = region_mods,
                             n_bootstraps = 100,
                             n_cores = 24)

save(region_mods,
     all_regions,
     network_mods,
     file = './BBS_results/Binhurdle.rda')
