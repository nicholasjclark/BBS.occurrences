library(nlme)
library(dplyr)
beta_data <- all_regions %>%
  dplyr::select(Latitude, Longitude, Site.group,Prop.forest,Tot.precip.yr,Year) %>%
  dplyr::mutate(x = (111.13 * Longitude * cos(34.021 * 3.14/180)),
                y = 111.13 * Latitude,
                prop.for.sc = as.vector(scale(Prop.forest)),
                precip.sc = as.vector(scale(Tot.precip.yr))) %>%
  dplyr::filter(Site.group == 'Mississippi Flyway')
beta_data <- beta_data[complete.cases(beta_data),]

# linear model, no autocorrelation
lin.mod <- gls(prop.for.sc ~ precip.sc, data = beta_data)
summary(lin.mod)
vario2 <- Variogram(lin.mod, form = ~x + y, resType = "pearson")
plot(vario2, smooth = TRUE)

# exponential autocorrelation, no nugget
# grouped by year (different years are assumed to be uncorrelated) as
# it is difficult to detect temporal autocorrelation without at least 20 - 25 years worth of data
# see (Teller et al., 2016) at
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12486)
exp.mod <- gls(prop.for.sc ~ precip.sc,
               correlation = corExp(form = ~x + y | Year), data = beta_data)
summary(exp.mod)
anova(lin.mod, exp.mod)

# exponential autocorrelation with nugget
nugg.mod <- gls(prop.for.sc ~ precip.sc,
               correlation = corExp(form = ~x + y | Year, nugget = TRUE), data = beta_data)
anova(lin.mod, exp.mod, nugg.mod)

# nugget not supported, examine the variogram from exp.mod
vario.exp <- Variogram(exp.mod, form = ~x + y, resType = "pearson")
plot(vario.exp, smooth = FALSE)

# plot normalized residual variogram, shouldn't see any pattern if autocorrelation
# is properly captured (use the asymptotic distance as the max)
vario.exp.norm <- Variogram(exp.mod, form = ~x + y, resType = "normalized", maxDist = 500)
plot(vario.exp.norm, smooth = FALSE)

resp <- as.matrix(cbind(beta_data[,c('Longitude','Latitude')], as.numeric(predict(exp.mod))))
e <- raster::extent(resp[,1:2])

r <- raster::raster(e, ncol=50, nrow=50)

# you need to provide a function 'fun' for when there are multiple points per cell
x <- raster::rasterize(resp[, 1:2], r, resp[,3], fun=mean)
library(raster)
plot(x, main = 'Predicted forest cover')

## Try with a spatial smoother (GAM)
library(mgcv)
m1_ml <- gls(yield ~ variety - 1, data = d, method = "ML")

m_gam1 <- gam(prop.for.sc ~ precip.sc, data = beta_data)
summary(m_gam1)

m_gam2 <- gam(prop.for.sc ~ precip.sc + te(x, y), data = beta_data)
summary(m_gam2)
AIC(m_gam1,m_gam2)
plot(m_gam2, pers = TRUE)

m_gam3 <- gam(prop.for.sc ~ precip.sc + s(x, y), data = beta_data)
summary(m_gam3)
plot(m_gam3, pers = TRUE)

m_gam4 <- gamm(prop.for.sc ~ precip.sc + te(x, y), random = list(Year = ~1),
               data = beta_data)
summary(m_gam4$lme)
AIC(m_gam1,m_gam2, m_gam3, m_gam4$lme)

## s function fits best, based on AIC metrics
vis.gam(m_gam3, plot.type="contour",
        view=c("x","y"), asp=1, type="response", contour.col="black", n.grid=100)
gam.check(m_gam3) #k is too low here
preds <- predict(m_gam3)
cor.test(beta_data$prop.for.sc, preds)

#because k was too low, may need to optimise
m_gam5 <- gam(prop.for.sc ~ precip.sc + s(x, y,k = 10), data = beta_data)
AIC(m_gam5)
vis.gam(m_gam5, plot.type="contour",color = 'heat',
        view=c("x","y"), asp=1, type="response", contour.col="black", n.grid=100)

# try with lat and long and try parallel computing
# lat and long results in lower AIC, different scales?
m_gam6 <- gam(prop.for.sc ~ precip.sc + s(Longitude, Latitude,k = 10),
              data = beta_data,
              control = gam.control(nthreads = 3))
plot(m_gam6, residuals=TRUE)
AIC(m_gam6)
vis.gam(m_gam6, plot.type="contour",color = 'heat',
        view=c("Longitude","Latitude"), asp=1, type="response",
        contour.col="black", n.grid=100)

gam.check(m_gam6)
preds.2 <- predict(m_gam5)
cor.test(beta_data$prop.for.sc, preds.2)

