Miss.bin.mod1 <- region_mods$`Mississippi Flyway1`$occurrence_mod
binplot1 <- plotBinomial(lassoBinomial = Miss.bin.mod1, cutoff = 0.1)

Miss.abund.mod1 <- region_mods$`Mississippi Flyway1`$abundance_mod
abundplot1 <- plotGaussian(lassoAbund = Miss.abund.mod1, cutoff = 0.1)
gridExtra::grid.arrange(binplot1, abundplot1, ncol = 2)

###2
Miss.bin.mod2 <- region_mods$`Mississippi Flyway2`$occurrence_mod
binplot2 <- plotBinomial(lassoBinomial = Miss.bin.mod2, cutoff = 0.1)

Miss.abund.mod2 <- region_mods$`Mississippi Flyway2`$abundance_mod
abundplot2 <- plotGaussian(lassoAbund = Miss.abund.mod2, cutoff = 0.1)
gridExtra::grid.arrange(binplot2, abundplot2, ncol = 2)

###3
Miss.bin.mod3 <- region_mods$`Mississippi Flyway3`$occurrence_mod
binplot3 <- plotBinomial(lassoBinomial = Miss.bin.mod3, cutoff = 0.1)

Miss.abund.mod3 <- region_mods$`Mississippi Flyway3`$abundance_mod
abundplot3 <- plotGaussian(lassoAbund = Miss.abund.mod3, cutoff = 0.1)
gridExtra::grid.arrange(binplot3, abundplot3, ncol = 2)

###4
Miss.bin.mod4 <- region_mods$`Mississippi Flyway4`$occurrence_mod
binplot4 <- plotBinomial(lassoBinomial = Miss.bin.mod4, cutoff = 0.1)

Miss.abund.mod4 <- region_mods$`Mississippi Flyway4`$abundance_mod
abundplot4 <- plotGaussian(lassoAbund = Miss.abund.mod4, cutoff = 0.1)
gridExtra::grid.arrange(binplot4, abundplot4, ncol = 2)

# Find covariates that were removed
region_mods$`Mississippi Flyway1`$removed_covs
region_mods$`Mississippi Flyway2`$removed_covs
region_mods$`Mississippi Flyway3`$removed_covs
region_mods$`Mississippi Flyway4`$removed_covs

# Inspect occurence model fits
region_mods$`Mississippi Flyway1`$occurrence_spec
region_mods$`Mississippi Flyway1`$occurrence_sens
region_mods$`Mississippi Flyway2`$occurrence_spec
region_mods$`Mississippi Flyway2`$occurrence_sens
region_mods$`Mississippi Flyway3`$occurrence_spec
region_mods$`Mississippi Flyway3`$occurrence_sens
region_mods$`Mississippi Flyway4`$occurrence_spec
region_mods$`Mississippi Flyway4`$occurrence_sens

# Inspect abundance model fits
region_mods$`Mississippi Flyway1`$abund_Rsquared
region_mods$`Mississippi Flyway2`$abund_Rsquared
region_mods$`Mississippi Flyway3`$abund_Rsquared
region_mods$`Mississippi Flyway4`$abund_Rsquared

# Summarise covariate influences on abundances
summariseAbund_res(Miss.abund.mod1)
summariseAbund_res(Miss.abund.mod2)
summariseAbund_res(Miss.abund.mod3)
summariseAbund_res(Miss.abund.mod4)

#### Predict centrality function ####
#function(lassoAbund, abundance_cent){
#run phylo and functional glmms to find predictors of centrality
abundance_cent <- region_mods$`Mississippi Flyway1`$abundance_cent
#abundance_cent <- lassAbund$abundance_cent

lassoAbund <- region_mods$`Mississippi Flyway1`$abundance_mod

sp_names <- rownames(lassoAbund$graph)

## Find species with high overall centrality
high_cent <- as.numeric(quantile(colMeans(abundance_cent[ , sp_names]),
                                    probs = 0.8))
keystone_sp <- names(abundance_cent[, sp_names][which(colMeans(abundance_cent[ , sp_names])
                                                      >= high_cent)])

## Gather centrality observations into tidy long format
library(dplyr)
cent_mod_data = abundance_cent %>%
  tidyr::gather('species' = sp_names, species, centrality) %>%
  dplyr::filter(centrality > 0) %>%
  dplyr::filter(species %in% keystone_sp) %>%
  #logit transform for skewed proportional data
  dplyr::mutate(centrality = log(centrality / (1 - centrality)))

hist(cent_mod_data$centrality)

fixed_effects <- paste(colnames(cent_mod_data[, !names(cent_mod_data) %in% c('species','centrality')]),
                       collapse='+')

full_formula <- paste(paste('centrality','~', fixed_effects), '(1 | species)', sep = " + ")

all_cvs <- parallel::mclapply(seq_len(10), function(j){
cat('Processing fold', j, 'of', 10, '...\n')
n_folds <- 5
folds <- caret::createFolds(rownames(cent_mod_data), n_folds)

fit_cvs <- lapply(seq_len(n_folds), function(k){
cat('Processing training validation', k, 'of', n_folds, '...\n')
train.mod <- lme4::lmer(as.formula(full_formula),
                        data = cent_mod_data[-folds[[k]], ], REML = FALSE)
test.data <- cent_mod_data[folds[[k]], ]

## calculate rsquared
test.data$predict = predict(train.mod, newdata = test.data)
rsquared <- cor.test(test.data$centrality, test.data$predict)[[4]]

## calculate rsquared
test.data$predict = predict(train.mod, newdata = test.data)
rsquared <- cor.test(test.data$centrality, test.data$predict)[[4]]

## extract coefficients
mod.results <- data.frame(coef(train.mod)$species)
intercepts <- data.frame(Species = rownames(mod.results),
                         Intercept = mod.results[,1])

mod.coefs <- data.frame(Parameter = rownames(data.frame(summary(train.mod)$coefficients[-1,])),
                        Coefficient = summary(train.mod)$coefficients[-1,1])

list(coefficients = mod.coefs,
     intercepts = intercepts,
     rsquared = rsquared)
})

list(coef.fold.summary = do.call(rbind, purrr::map(fit_cvs, 'coefficients')))
}, mc.cores = 3)

coef.summary <- do.call(rbind, purrr::map(all_cvs, 'coef.fold.summary')) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(Lower.95 = quantile(Coefficient, probs = 0.01),
                   Median = quantile(Coefficient, probs = 0.5),
                   Upper.95 = quantile(Coefficient, probs = 0.99)) %>%
  dplyr::arrange(-abs(Median))

#### Summarise covariates function ####
summariseAbund_res = function(lassoAbund){

## Summarise interaction coefficients
  graph <- lassoAbund$graph
  coefs <- graph[upper.tri(graph)]
  coefs.nonzero <- coefs[which(coefs != 0 )]
  coefs.pos <- length(coefs[which(coefs > 0 )])
  coefs.neg <- length(coefs[which(coefs < 0 )])

## Summarise direct influences of covariates on species' abundances
cov.dir.inf <- t(rev(sort(lassoAbund$direct_coefs %>%
                            dplyr::select_if(names(.) %in% names(lassoAbund$indirect_coefs)) %>%
                            dplyr::summarise_all(funs(length(.[which(. != 0 )]))) %>%
                            dplyr::select_if(. > 0))))

cov.dir.pos <- t(rev(sort(lassoAbund$direct_coefs %>%
                            dplyr::select_if(names(.) %in% rownames(cov.dir.inf)) %>%
                            dplyr::summarise_all(funs(length(.[which(. > 0 )]))))))

cov.dir.mag <- t(rev(sort(lassoAbund$direct_coefs %>%
                            dplyr::select_if(names(.) %in% rownames(cov.dir.inf)) %>%
                            dplyr::summarise_all(funs(round(mean(abs(.[which(. != 0 )])),
                                                            4))))))


cov.overview <- dplyr::bind_cols(data.frame(cov.dir.mag),
                                 data.frame(cov.dir.inf),
                                 data.frame(cov.dir.pos))
rownames(cov.overview) <- rownames(cov.dir.inf)
colnames(cov.overview) <- c('Mean.abs.magnitude', 'Number.nonzero.coefs',
                            'Number.pos.coefs')
cov.overview$Number.neg.coefs <- cov.overview$Number.nonzero.coefs -
  cov.overview$Number.pos.coefs
cov.overview$Number.species <- nrow(lassoAbund$direct_coefs)

## Summarise indirect influences of covariates
cov.inf <- lapply(seq_along(lassoAbund$indirect_coefs), function(x){
  graph <- lassoAbund$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  mean.direct.coefs <- lassoAbund$graph[upper.tri(lassoAbund$graph)]
  length(coefs[which(coefs != 0 )])
})
names(cov.inf) <- names(lassoAbund$indirect_coefs)

cov.mag <- lapply(seq_along(lassoAbund$indirect_coefs), function(x){
  graph <- lassoAbund$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  mean.direct.coefs <- lassoAbund$graph[upper.tri(lassoAbund$graph)]
  round(mean(abs(coefs[which(coefs != 0 )])), 4)
})
names(cov.mag) <- names(lassoAbund$indirect_coefs)

cov.pos <- lapply(seq_along(lassoAbund$indirect_coefs), function(x){
  graph <- lassoAbund$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  coefs.nonzero <- coefs[which(coefs != 0 )]
  length(coefs.nonzero[which(coefs.nonzero > 0)])
})
names(cov.pos) <- names(lassoAbund$indirect_coefs)

cov.indirect.overview <- dplyr::bind_cols(data.frame(unlist(cov.mag)),
                                          data.frame(unlist(cov.inf)),
                                          data.frame(unlist(cov.pos)))
rownames(cov.indirect.overview) <- names(unlist(cov.inf))
colnames(cov.indirect.overview) <- c('Mean.abs.magnitude',
                                     'Number.nonzero.coefs',
                                     'Number.pos.coefs')
cov.indirect.overview$Number.neg.coefs <- cov.indirect.overview$Number.nonzero.coefs -
  cov.indirect.overview$Number.pos.coefs

cov.indirect.overview$Number.interactions <- length(coefs.nonzero)

cov.indirect.overview <- cov.indirect.overview[order(-cov.indirect.overview$Mean.abs.magnitude),]
return(list(tot_pos_interactions = coefs.pos,
            tot_neg_interactions = coefs.neg,
            direct_effects = cov.overview,
            indirect_effects = cov.indirect.overview))
}
