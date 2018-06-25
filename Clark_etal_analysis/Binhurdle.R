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

# Test models on a single flyway (Mississippi)
# Gather names of flyways
flyways <- as.character(unique(all_regions$Site.group))

# Fit the hurdle model (10 reps) for each flyway
region_mods <- lapply(seq_len(length(flyways)), function(x){
  hurdleModel(flyway = flyways[x], n_cores = 24)
})
names(region_mods) <- flyways

save(region_mods,
     all_regions,
     file = './BBS_results/Binhurdle.rda')
