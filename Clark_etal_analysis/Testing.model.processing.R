Miss.bin.mod1 <- region_mods$`Mississippi Flyway1`$occurrence_mod
binplot1 <- plotBinomial(lassoBinomial = Miss.bin.mod1, cutoff = 0.2)

Miss.abund.mod1 <- region_mods$`Mississippi Flyway1`$abundance_mod
abundplot1 <- plotGaussian(lassoAbund = Miss.abund.mod1, cutoff = 0.2)
gridExtra::grid.arrange(binplot1, abundplot1, ncol = 2)

###2
Miss.bin.mod2 <- region_mods$`Mississippi Flyway2`$occurrence_mod
binplot2 <- plotBinomial(lassoBinomial = Miss.bin.mod2, cutoff = 0.2)

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
library(dplyr)
summariseAbund_res(Miss.abund.mod1)
summariseAbund_res(Miss.abund.mod2)
summariseAbund_res(Miss.abund.mod3)
summariseAbund_res(Miss.abund.mod4)

#### Predict centrality function ####
#function(lassoAbund, abundance_cent){
#run phylo and functional glmms to find predictors of centrality
abundance_cent <- region_mods$`Mississippi Flyway1`$network_metrics
#abundance_cent <- lassAbund$abundance_cent

lassoAbund <- region_mods$`Mississippi Flyway1`$abundance_mod

sp_names <- rownames(lassoAbund$graph)

## Find species with high overall centrality
high_cent <- as.numeric(quantile(colMeans(abundance_cent[ , sp_names]),
                                    probs = 0.75))
keystone_sp <- names(abundance_cent[, sp_names][which(colMeans(abundance_cent[ , sp_names])
                                                      >= high_cent)])

## Gather centrality observations into tidy long format
cent_mod_data = abundance_cent %>%
  tidyr::gather('species' = sp_names, species, centrality) %>%
  dplyr::filter(centrality > 0) %>%
  dplyr::filter(species %in% keystone_sp) %>%

  #logit transform for skewed proportional data
  dplyr::mutate(centrality = log(centrality / (1 - centrality))) %>%
  dplyr::select(-beta_os) %>%
  dplyr::filter(is.finite(centrality))

## Load species' traits informationn to be used as predictors of centrality
data("Bird.traits")

cent_mod_data = cent_mod_data %>%
  dplyr::left_join(Bird.traits)

## Use the covariate names to build a linear mixed model formula
fixed_effects <- paste(colnames(cent_mod_data[, !names(cent_mod_data) %in% c('species','centrality')]),
                       collapse = '+')

full_formula <- paste(paste('centrality', '~', fixed_effects),
                      '(1 | species)', sep = " + ")

#### Perform 100 rounds of 4-fold cv to estimate coefficient quantiles and rsquared ####
all_cvs <- parallel::mclapply(seq_len(5), function(j){

fit_cvs <- lapply(seq_len(4), function(k){

## extract 75% of observations as training data
row_indices <- seq_len(nrow(cent_mod_data))
in_train <- sample(row_indices, floor(nrow(cent_mod_data) * .75), F)

## run the model on training data
train_mod <- lmerTest::lmer(as.formula(full_formula),
                        data = cent_mod_data[in_train, ], REML = FALSE)

## use backward step-wise selection to choose best predictors
step_test <- lmerTest::step(train_mod)

best_mod <- lmerTest::get_model(step_test)

## calculate rsquared by predicting with-held (test) data
test_data <- cent_mod_data[-in_train, ]
test_data$predict = predict(best_mod, newdata = test_data)
rsquared <- cor.test(test_data$centrality, test_data$predict)[[4]]

## extract coefficients
mod_results <- data.frame(coef(best_mod)$species)
intercepts <- data.frame(Species = rownames(mod_results),
                         Intercept = mod_results[,1])

mod_coefs <- data.frame(Parameter = rownames(data.frame(summary(best_mod)$coefficients[-1, ])),
                        Coefficient = summary(best_mod)$coefficients[-1, 1])

list(coefficients = mod_coefs,
     intercepts = intercepts,
     rsquared = rsquared)
})

list(coef.fold.summary = do.call(rbind, purrr::map(fit_cvs, 'coefficients')),
     intercept.summary = do.call(rbind, purrr::map(fit_cvs, 'intercepts')),
     rsquared.summary = do.call(rbind, purrr::map(fit_cvs, 'rsquared')))

}, mc.cores = 3)

#### Extract summary statistics of model coefficients ####
coef.summary <- do.call(rbind, purrr::map(all_cvs, 'coef.fold.summary')) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(Lower = round(min(Coefficient), 4),
                   Median = round(quantile(Coefficient, probs = 0.5), 4),
                   Upper = round(max(Coefficient), 4)) %>%
  dplyr::arrange(-abs(Median))

rsquared.summary <- data.frame(cor = do.call(rbind, purrr::map(all_cvs, 'rsquared.summary'))) %>%
  dplyr::mutate(Group = 1) %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(Lower = round(min(cor), 4),
                   Median = round(quantile(cor, probs = 0.5), 4),
                   Upper = round(max(cor), 4)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-Group)

intercept.summary <- do.call(rbind, purrr::map(all_cvs, 'intercept.summary')) %>%
  dplyr::group_by(Species) %>%
  dplyr::mutate(Intercept = boot::inv.logit(Intercept)) %>%
  dplyr::summarise(Lower = round(min(Intercept), 4),
                   Median = round(quantile(Intercept, probs = 0.5), 4),
                   Upper = round(max(Intercept), 4)) %>%
  dplyr::arrange(-abs(Median))


#### Estimate predictors of interaction beta diversity ####
beta_mod_data <- abundance_cent[,!(names(abundance_cent) %in% sp_names)] %>%
  dplyr::mutate(beta_os = as.vector(scale(beta_os)))

## Build a linear model formula
fixed_effects <- paste(colnames(beta_mod_data[, !names(beta_mod_data) %in% c('beta_os')]),
                       collapse = '+')

full_formula <- paste('beta_os', '~', fixed_effects)

## extract 75% of observations as training data
row_indices <- seq_len(nrow(beta_mod_data))
in_train <- sample(row_indices, floor(nrow(beta_mod_data) * .75), F)

train_mod <- lm(as.formula(full_formula),
                            data = beta_mod_data[in_train, ])

## use backward step-wise selection to choose best predictors
step_test <- step(train_mod, trace = 0)

## calculate rsquared by predicting with-held (test) data
test_data <- beta_mod_data[-in_train, ]
test_data$predict = predict(step_test, newdata = test_data)
rsquared <- cor.test(test_data$beta_os, test_data$predict)[[4]]

summary(step_test)

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
