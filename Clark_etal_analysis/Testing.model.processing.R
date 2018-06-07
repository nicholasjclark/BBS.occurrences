plotBinomial(lassoBinomial = Miss.bin.mod, cutoff = 0.05)

## Repeat for abundances
abund.coefs <- Miss.crf.mod$direct_coefs[, sp.names]
abund.influence.vec <- abund.coefs[upper.tri(abund.coefs)]
abund.nonzeros <- which(abund.influence.vec != 0)
abund.results <- data.frame(coef(summary(lm(abund.influence.vec[abund.nonzeros] ~
                                              phy.dist.vec[abund.nonzeros]^2 +
                                              ecol.dist.vec[abund.nonzeros]^2))))
abund.results <- round(abund.results, 4)

## Which covariate has the strongest influence on abundance covariance?
cov.influence <- which.max(unlist(lapply(Miss.crf.mod$indirect_coefs, function(x){
  sum(abs(data.frame(x)))
})))

View(Miss.crf.mod$indirect_coefs[[cov.influence]])

inf.cov.coefs <- data.frame(Miss.crf.mod$indirect_coefs[[cov.influence]])[, sp.names]
inf.cov.influence.vec <- inf.cov.coefs[upper.tri(inf.cov.coefs)]
inf.cov.nonzeros <- which(inf.cov.influence.vec != 0)
inf.cov.results <- data.frame(coef(summary(lm(inf.cov.influence.vec[inf.cov.nonzeros] ~
                                                phy.dist.vec[inf.cov.nonzeros] +
                                                ecol.dist.vec[inf.cov.nonzeros]))))
inf.cov.results <- round(inf.cov.results, 4)
