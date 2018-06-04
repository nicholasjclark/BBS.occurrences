library(BBS.occurrences)
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
Miss.covs <- Site.descriptors[rows.MissFly, c(22,28:33,35:46)]
library(dplyr)
Miss.covs %>%
  dplyr::mutate(Prop.Veg = Prop.forest + Prop.shrubland) %>%
  dplyr::select(-Prop.forest, -Prop.shrubland) -> Miss.covs

cols_remove <- caret::findCorrelation(cor(na.omit(Miss.covs)),
                                      cutoff = .9, exact = T)
names(Miss.covs)[cols_remove]
Miss.covs <- Miss.covs[, -cols_remove]

# a small number of year * site observations do not have landcover values
# here we remove these observations from all datasets
row.has.na <- apply(Miss.covs, 1, function(x){any(is.na(x))})
Miss.covs <- as.data.frame(Miss.covs[!row.has.na, ])
Miss.bin <- as.data.frame(Miss.bin[!row.has.na, ])
Miss.abund <- as.data.frame(Miss.abund[!row.has.na, ])

## Run occurrence and abundance models
Miss.bin.mod <- lassoBinomial_comm(outcome_data = Miss.bin[,1:50],
                                   count_data = Miss.abund[,1:50],
                                   covariates = Miss.covs,
                                   n_reps = 5, n_cores = 3)

Miss.crf.mod <- lassoAbund_comm(outcome_data = Miss.abund[,1:50],
                                binary_data = Miss.bin[,1:50],
                                covariates = Miss.covs,
                                n_reps = 5,
                                n_cores = 3)

## Import phylogenetic and ecological distance datasets
data("Bird.ecol.distances")
data("Bird.phy.distances")

## Extract distances for the assessed species
sp.names <- colnames(Miss.crf.mod$graph)
Bird.ecol.distances <- Bird.ecol.distances[sp.names, sp.names]
ecol.dist.vec <- Bird.ecol.distances[upper.tri(Bird.ecol.distances)]

Bird.phy.distances <- Bird.phy.distances[sp.names, sp.names]
phy.dist.vec <- Bird.phy.distances[upper.tri(Bird.phy.distances)]

## Extract vector of occurrence influences
occur.coefs <- Miss.bin.mod$Coefficients[, sp.names]
occur.influence.vec <- occur.coefs[upper.tri(occur.coefs)]

## Run a simple linear model for occurrences
occur.results <- data.frame(coef(summary(lm(occur.influence.vec ~ phy.dist.vec + ecol.dist.vec))))

## Repeat for abundances
abund.coefs <- Miss.crf.mod$direct_coefs[, sp.names]
abund.influence.vec <- abund.coefs[upper.tri(abund.coefs)]
abund.results <- data.frame(coef(summary(lm(abund.influence.vec ~ phy.dist.vec + ecol.dist.vec))))

save(Miss.bin.mod,Miss.crf.mod,
     occur.results, abund.results,
     file = './BBS_results/Binhurdle.rda')
