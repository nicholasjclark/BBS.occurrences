library(dplyr)
library(ggplot2)
map.usa_states <- map_data("world", c("USA","hawaii"), xlim=c(-180,-65), ylim=c(19,72),interior = FALSE)

all_regions %>%
  group_by(Site.group) %>%
  summarise(n.sites = n_distinct(Flyway.region))

all_regions_distinct = all_regions %>%
  group_by(Site.group) %>%
  dplyr::select(Flyway.region:Longitude) %>%
  ungroup() %>%
  distinct()

geom.text.size <- 4
theme.size <- (14/5) * geom.text.size
ggplot(data = all_regions_distinct) +
  geom_polygon(data = map.usa_states,
               aes(x = long, y = lat, group = group), fill = "#484848") +
  facet_wrap(~ Site.group) +
  geom_point(data = all_regions, aes(x = Longitude, y = Latitude,
                                        colour = as.factor(Flyway.region)),
             size = 0.4, alpha = 0.6)  +
  theme_bw() +
  theme(legend.position = "") +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1),
                 axis.text = ggplot2::element_text(size = 9),
                 legend.title = ggplot2::element_blank()) +
  ggplot2::coord_map(xlim = c(-169 , -65),
                     ylim = c(22, 71.5)) +
  ggplot2::labs(y = "Latitude", x = "Longitude") +
  ggplot2::theme(axis.text = ggplot2::element_text(size = theme.size - 1, colour = "black"),
                                              strip.text.x = ggplot2::element_text(size = theme.size, colour = "black"),
                                              axis.title = ggplot2::element_text(size = theme.size, colour = 'black'))
