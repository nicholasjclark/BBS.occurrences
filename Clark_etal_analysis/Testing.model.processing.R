network_mods$sp_cent_coefs
network_mods$betadiv_coefs
network_mods$betadiv_intercepts
network_mods$sp_cent_rsquared
network_mods$betadiv_rsquared
network_mods$sp_cent_intercepts


beta_fixed_effects <- paste(network_mods$betadiv_coefs$Parameter,
                            collapse = '+')

beta_full_formula <- paste('beta_os', '~', beta_fixed_effects)

mod <- lm(as.formula(beta_full_formula),
          data = network_mods$betadiv_mod_data)

## Predict outcomes using 99% confidence and 68% prediction intervals
mod.predict <- cbind(network_mods$betadiv_mod_data,
                     predict(mod, interval = 'confidence',
                             level = 0.9999999999))
predictions <- predict(mod, interval = 'prediction',
                       level = 0.68)

colnames(predictions) <- c('pred.fit','pred.lwr','pred.upr')
mod.predict <- cbind(mod.predict,predictions)

geom.text.size <- 4
theme.size <- (14/5) * geom.text.size

## Plot the predicted regression line and confidence interval for
# habitat niche distance, smooth prediction intervals using loess
g1 <- ggplot2::ggplot(mod.predict,
                      ggplot2::aes(Prop.forest)) +
  ggplot2::stat_smooth(ggplot2::aes(y = lwr), size = 0.1,
                       colour = "white", n = 500, span = 1, se = FALSE) +
  ggplot2::stat_smooth(ggplot2::aes(y = upr), size = 0.1,
                       colour = "white", n = 500, span = 1, se = FALSE)
gg1 <- ggplot2::ggplot_build(g1)
df2 <- data.frame(x = gg1$data[[1]]$x,
                  ymin = gg1$data[[1]]$y,
                  ymax = gg1$data[[2]]$y)

g2 <- ggplot2::ggplot(mod.predict,
                     ggplot2::aes(Prop.forest)) +
  ggplot2::stat_smooth(ggplot2::aes(y = pred.lwr), size = 0.1,
                       colour = "white", n = 500, span = 1, se = FALSE) +
  ggplot2::stat_smooth(ggplot2::aes(y = pred.upr), size = 0.1,
                       colour = "white", n = 500, span = 1, se = FALSE)
gg2 <- ggplot2::ggplot_build(g2)
df3 <- data.frame(x = gg2$data[[1]]$x,
                  ymin = gg2$data[[1]]$y,
                  ymax = gg2$data[[2]]$y)

 forest.plot <- g1 +
   geom_point(data = mod.predict, aes(x = Prop.forest, y = beta_os),
              alpha = 0.1, size = 0.15, colour = 'grey60') +
  ggplot2::geom_ribbon(data = df3,
                        ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "red", alpha = 0.3) +
  ggplot2::geom_ribbon(data = df2,
                       ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
                       fill = "red", alpha = 0.85) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1),
                 axis.text = ggplot2::element_text(size = 9),
                 legend.title = ggplot2::element_blank())+
  ggplot2::ylab(expression(paste("Interaction dissimilarity (", beta["os"],")"))) +
  ggplot2::xlab('Percent forest cover') +
  ggplot2::theme(axis.text = ggplot2::element_text(size = theme.size - 1, colour = "black"),
                 strip.text.x = ggplot2::element_text(size = theme.size, colour = "black"),
                 axis.title = ggplot2::element_text(size = theme.size, colour = 'black')) +
  ggplot2::theme(plot.title = ggplot2::element_text(color="black", face="bold", hjust = 0.5,
                                                    size = theme.size))

pdf('forest.plot.pdf',width = 3.1,height = 2.8)
forest.plot
dev.off()

Miss.bin.mod1 <- region_mods$`Mississippi Flyway1`$occurrence_mod
View(Miss.bin.mod1$coefficients)
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
