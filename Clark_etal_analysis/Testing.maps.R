install.packages("fiftystater")
library(fiftystater)
library(dplyr)
points <- Site.descriptors %>%
  dplyr::ungroup() %>%
  dplyr::filter(Year == 2008) %>%
  dplyr::select(Latitude, Longitude, Prop.anthro.urban) %>%
  dplyr::distinct()

points.spdf <- sp::SpatialPointsDataFrame(coords = points[, c('Longitude','Latitude')],
                                          data = points,
                                        proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))
fortify.points <- data.frame(points.spdf)

map.usa_states <- map_data("world", c("USA","hawaii"), xlim=c(-180,-65), ylim=c(19,72),interior = FALSE)


ggplot() +
  geom_polygon(data = map.usa_states,
               aes(x = long, y = lat, group = group), fill = "#484848") +
  geom_point(data = fortify.points, aes(x = Longitude, y = Latitude,
                                        colour = Prop.anthro.urban))
