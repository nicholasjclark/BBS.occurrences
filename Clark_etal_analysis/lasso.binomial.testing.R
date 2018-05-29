
data("BBS.occurrences")
data("BBS.abundances")
data("Site.descriptors")
rows.MissFly <- which(Site.descriptors$Site.group == 'Mississippi Flyway')
Miss.abund <- BBS.abundances[rows.MissFly, ]

# remove species occurring in fewer than 25% of observations
abund.zero.cols <- which((colSums(Miss.abund) / nrow(Miss.abund)) < 0.25)
Miss.abund <- Miss.abund[, -abund.zero.cols]

# repeat for binary occurrence data
Miss.bin <- BBS.occurrences[rows.MissFly, ]
Miss.bin <- Miss.bin[, -abund.zero.cols]

# also remove un-needed rows and cols from covariates data
# note, year is not included as temporal autocorrelation will be
# difficult to detect without at least 20 - 25 years worth of data
# see (Teller et al., 2016) at
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12486)
Miss.covs <- Site.descriptors[rows.MissFly, c(22,24,28:33,35:46)]

# a small number of year * site observations do not have landcover values
# here we remove these observations from all datasets
row.has.na <- apply(Miss.covs, 1, function(x){any(is.na(x))})
Miss.covs <- as.data.frame(Miss.covs[!row.has.na, ])
Miss.bin <- as.data.frame(Miss.bin[!row.has.na, ])
Miss.abund <- as.data.frame(Miss.abund[!row.has.na, ])

Miss.covs.biotic <- cbind(Miss.covs, Miss.bin[,c(1:10)])

begin <- Sys.time()
test <- lassoBinomial_comm(outcome_data = Miss.bin,
                           outcome_indices = c(1:10),
                           covariates = as.matrix(Miss.covs.biotic),
                           n_reps = 5)
Sys.time() - begin
