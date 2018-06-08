Miss.bin.mod1 <- region_mods$`Mississippi Flyway1`$occurrence_mod
binplot1 <- plotBinomial(lassoBinomial = Miss.bin.mod1, cutoff = 0.1)

Miss.abund.mod1 <- region_mods$`Mississippi Flyway1`$abundance_mod
abundplot1 <- plotGaussian(lassoAbund = Miss.abund.mod1, cutoff = 0.1)
gridExtra::grid.arrange(binplot1, abundplot1, ncol = 2)

###2
Miss.bin.mod2 <- region_mods$`Mississippi Flyway2`$occurrence_mod
binplot2 <- plotBinomial(lassoBinomial = Miss.bin.mod2, cutoff = 0)

Miss.abund.mod2 <- region_mods$`Mississippi Flyway2`$abundance_mod
abundplot2 <- plotGaussian(lassoAbund = Miss.abund.mod2, cutoff = 0)
gridExtra::grid.arrange(binplot2, abundplot2, ncol = 2)

###3
Miss.bin.mod3 <- region_mods$`Mississippi Flyway3`$occurrence_mod
binplot3 <- plotBinomial(lassoBinomial = Miss.bin.mod3, cutoff = 0)

Miss.abund.mod3 <- region_mods$`Mississippi Flyway3`$abundance_mod
abundplot3 <- plotGaussian(lassoAbund = Miss.abund.mod3, cutoff = 0)
gridExtra::grid.arrange(binplot3, abundplot3, ncol = 2)

###4
Miss.bin.mod4 <- region_mods$`Mississippi Flyway4`$occurrence_mod
binplot4 <- plotBinomial(lassoBinomial = Miss.bin.mod4, cutoff = 0)

Miss.abund.mod4 <- region_mods$`Mississippi Flyway4`$abundance_mod
abundplot4 <- plotGaussian(lassoAbund = Miss.abund.mod4, cutoff = 0)
gridExtra::grid.arrange(binplot4, abundplot4, ncol = 2)


summariseAbund_res(Miss.abund.mod1)
summariseAbund_res(Miss.abund.mod2)
summariseAbund_res(Miss.abund.mod3)
summariseAbund_res(Miss.abund.mod4)

#### Summarise covariates function ####
summariseAbund_res = function(lassoAbund){

## Summarise interaction coefficients
library(dplyr)
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

cov.indirect.overview <- cov.indirect.overview[order(-cov.indirect.overview$Number.nonzero.coefs),]
return(list(tot_pos_interactions = coefs.pos,
            tot_neg_interactions = coefs.neg,
            direct_effects = cov.overview,
            indirect_effects = cov.indirect.overview))
}
