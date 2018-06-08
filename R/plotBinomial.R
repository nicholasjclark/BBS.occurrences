#'Construct regression plots of species' pairwise influences on occurrence probabilities
#'using model output from \code{lassoBinomial}
#'
#'
#'
#'@param lassoBinomial A model object returned from \code{lassoBinomial_comm}
#'@param cutoff Positive numeric value indicating the minimum quantile of biotic coefficients'
#'absolute values (excluding \code{0}'s) that are considered as reflective of an
#'important influence on occurrence probability. Default is \code{0.05}, meaning that the top 95
#'percent of biotic coefficients are retained when performing the polynomial regression
#'@details Species' pairwise influences on occurrence probabilities are extracted from
#'\code{lassoBinomial$Coefficients} and predicted
#'by their respective habitat niche and phylogenetic distances (extracted from
#'\code{data("Bird.ecol.distances")} and \code{data("Bird.phy.distances")}) using
#'a second-order polynomial regression of the form
#'\code{lm(Influence ~ poly(Phy.dist,2) + poly(Ecol.dist,2))}. Only pairwise
#'influences whose absolute magnitudes are above the specified \code{cutoff} quantile
#'are used in the regression
#'@return A \code{ggplot2} object visualising predicted habitat niche and phylogenetic distance
#'regression lines and 99% confidence intervals, along with the proportion of explained variance
#'captured by each predictor
#'
#'@export
plotBinomial = function(lassoBinomial, cutoff){

#### Import phylogenetic and ecological distance datasets ####
data("Bird.ecol.distances")
data("Bird.phy.distances")

#### Extract ecological and phylogenetic pairwise distances for the assessed species ####
sp.names <- rownames(lassoBinomial$Coefficients)
Bird.ecol.distances <- Bird.ecol.distances[sp.names, sp.names]
ecol.dist.vec <- Bird.ecol.distances[upper.tri(Bird.ecol.distances)]

Bird.phy.distances <- Bird.phy.distances[sp.names, sp.names]
phy.dist.vec <- Bird.phy.distances[upper.tri(Bird.phy.distances)]

## Extract vector of occurrence influence coefficients
occur.coefs <- lassoBinomial$Coefficients[, sp.names]
influence.vec <- occur.coefs[upper.tri(occur.coefs)]

## Calculate quantiles of non-zero influence coefficients
influence.vec <- influence.vec[which(influence.vec != 0)]

if(missing(cutoff)){
  cutoff <- as.numeric(quantile(abs(influence.vec),
                                                 c(0.05)))
} else {
  cutoff <- as.numeric(quantile(abs(influence.vec),
                                c(cutoff)))
}

## Keep only the coefficients that are above the specified cutoff magnitude
occur.nonzeros <- which(abs(influence.vec) > cutoff)

#### Run a simple second-order polynomial model for occurrences ####
occur.data <- data.frame(Influence = influence.vec[occur.nonzeros],
                         Phy.dist = phy.dist.vec[occur.nonzeros],
                         Ecol.dist = ecol.dist.vec[occur.nonzeros])

mod <- lm(Influence ~ poly(Phy.dist, 2) +
            poly(Ecol.dist, 2), data = occur.data)

## Extract proportions of explained variance for each distance predictor
Betas <- (mod$coefficients[2:5]^2 / sum(mod$coefficients[2:5]^2))
Var.exp <- Betas * as.numeric(summary(mod)[8])

## Predict outcomes using 99% confidence intervals
mod.predict <- cbind(occur.data,
                     predict(mod, interval = 'confidence',
                             level = 0.99))

## Set y limits for plotting
ceiling_dec <- function(x, level = 1) round(x + 5*10 ^ (-level - 1), level)
y.upper <- ceiling_dec(max(mod.predict$upr), 1)
y.lower <- -1 * ceiling_dec(max(abs(mod.predict$lwr)), 1)

#### Plotting prediction lines ####
## Set the default text sizes
geom.text.size <- 4
theme.size <- (14/5) * geom.text.size

## Plot the predicted regression line and confidence interval for
# habitat niche distance, smooth prediction intervals using loess
g1 <- ggplot2::ggplot(mod.predict,
                      ggplot2::aes(Ecol.dist)) +
  ggplot2::stat_smooth(ggplot2::aes(y = lwr), size = 0.1,
                       colour = "black", n = 200, span = 1,
                       method = "loess", se = FALSE) +
  ggplot2::stat_smooth(ggplot2::aes(y = upr), size = 0.1,
                       colour = "black", n = 200, span = 1,
                       method = "loess", se = FALSE) +
  ggplot2::stat_smooth(ggplot2::aes(y = fit), size = 0.6,
                       colour = "black", n = 200, span = 1,
                       method = "loess", se = FALSE)
gg1 <- ggplot2::ggplot_build(g1)
df2 <- data.frame(x = gg1$data[[1]]$x,
                  ymin = gg1$data[[1]]$y,
                  ymax = gg1$data[[2]]$y)

## Create the text for visualising variance explained
hab_text <- substitute(italic(Var[exp])~"="~r2,
                       list(r2 = format(round(sum(Var.exp[3:4]), 3),
                                        digits = 3)))
dftext <- data.frame(Ecol.dist = -Inf, y = Inf,
                     hab_text = as.character(as.expression(hab_text)))

## Build the habitat niche plot
hab.plot <- g1 +
  ggplot2::geom_ribbon(data = df2,
                       ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
                       fill = "blue", alpha = 0.4) +
  ggplot2::ylim(y.lower, y.upper) +
  ggplot2::xlim(0, 1) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1),
                 axis.text = ggplot2::element_text(size = 9),
                 legend.title = ggplot2::element_blank())+
  ggplot2::labs(y = "", x = "") +
  ggplot2::geom_text(ggplot2::aes(label = hab_text, y = y, x = Ecol.dist),
                     data = dftext, parse = TRUE,
                     vjust = 1.5, hjust = -0.05) +
  ggplot2::theme(axis.text = ggplot2::element_text(size = theme.size - 1, colour = "black"),
                 strip.text.x = ggplot2::element_text(size = theme.size, colour = "black"),
                 axis.title = ggplot2::element_text(size = theme.size, colour = 'black')) +
  ggplot2::ggtitle("Habitat niche") +
  ggplot2::theme(plot.title = ggplot2::element_text(color="black", face="bold", hjust = 0.5,
                                                    size = theme.size))

## Repeat using the phylogenetic distance predictor
g1 <- ggplot2::ggplot(mod.predict, ggplot2::aes(Phy.dist)) +
  ggplot2::stat_smooth(ggplot2::aes(y = lwr), size = 0.1,
                       colour = "black", n = 200, span = 1,
                       method = "loess", se = FALSE) +
  ggplot2::stat_smooth(ggplot2::aes(y = upr), size = 0.1,
                       colour = "black", n = 200, span = 1,
                       method = "loess", se = FALSE) +
  ggplot2::stat_smooth(ggplot2::aes(y = fit), size = 0.6,
                       colour = "black", n = 200, span = 1,
                       method = "loess", se = FALSE)
gg1 <- ggplot2::ggplot_build(g1)
df2 <- data.frame(x = gg1$data[[1]]$x,
                  ymin = gg1$data[[1]]$y,
                  ymax = gg1$data[[2]]$y)
phy_text <- substitute(italic(Var[exp])~"="~r2,
                       list(r2 = format(round(sum(Var.exp[1:2]), 3),
                                        digits = 3)))
dftext <- data.frame(Phy.dist = -Inf, y = Inf,
                     phy_text = as.character(as.expression(phy_text)))

phy.plot <- g1 + ggplot2::geom_ribbon(data = df2,
                                      ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
                                      fill = "red", alpha = 0.4) +
  ggplot2::ylim(y.lower, y.upper) +
  ggplot2::xlim(0, 1) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1),
                 axis.text = ggplot2::element_text(size = 9),
                 legend.title = ggplot2::element_blank())+
  ggplot2::labs(y = "", x = substitute("Species' pairwise dissimilarity (" * {italic(N) == n} * ")",
                                       list(n = nrow(occur.data)))) +
  ggplot2::geom_text(ggplot2::aes(label = phy_text, y = y, x = Phy.dist),
                     data = dftext, parse = TRUE,
                     vjust = 1.5, hjust = -0.05) +
  ggplot2::theme(axis.text = ggplot2::element_text(size = theme.size - 1, colour = "black"),
                 strip.text.x = ggplot2::element_text(size = theme.size, colour = "black"),
                 axis.title = ggplot2::element_text(size = theme.size, colour = 'black')) +
  ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 7.5, r = 0,
                                                                      b = 0, l = 0))) +
  ggplot2::ggtitle("Evolutionary history") +
  ggplot2::theme(plot.title = ggplot2::element_text(color="black", face="bold", hjust = 0.5,
                                                    size = theme.size))

#### Combine the two plots as a grid object and return ####
return(gridExtra::grid.arrange(hab.plot, phy.plot, ncol = 1,
                               heights = c(1, 1),
                               left = grid::textGrob("Species' pairwise influence on occurrence",
                                                     rot = 90, vjust = 1.5)))

}

