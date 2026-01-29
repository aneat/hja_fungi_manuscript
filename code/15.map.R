##load packages
library("terra")
library("tidyterra")
library("tidyverse")
library("RColorBrewer")
library("sf")
library("usmap")
library("cowplot")
library("png")
library("ggpubr")
library("jpeg")
library("ggspatial")

##make a data frame from elevation raster to be used in ggplot
elevation <- terra::rast("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\hja map\\microclimate_rasters\\hja_dem_2016_1m_corrected.tif")
elevation
plot(elevation)
plot(sites, add = TRUE)

##create points representing each of the ten sites
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil.rds")
sites <- subset(env1, select = c(5, 10, 11, 12))
sites <- sites %>% mutate_at(vars(matches("longitude")), `-`)
sites <- sites[, c(1,2,4,3)]

hja_sites <- unique(sites) %>%
  st_as_sf(coords = c(3,4), crs = st_crs(4326)) %>%
  st_transform(crs = 26910)

hja_sites <- hja_sites %>%
  arrange(elevation) %>%
  mutate(numbers = 1:10)

hja_sites$numbers <- gsub("3", "3   ", hja_sites$numbers)

sites <- vect(hja_sites)

##initial map(s) of hja
site_cols <- brewer.pal(11, "PRGn")[1:6]
one <- ggplot(sites) +
 geom_spatraster(data = elevation) + geom_spatvector(aes(color = sites$elevation, stroke = 2, size = 8)) + theme_classic() + scale_fill_gradient(low = "grey82", high = "black",na.value = "transparent", name = "terrain elevation (m)") + scale_color_gradientn(colors = site_cols, name = "site elevation (m)") + theme(text = element_text(size = 23), legend.position = c(.91, .79)) + scale_size(guide = "none") + annotation_scale(plot_unit = "m", text_cex = 1.5) + annotation_north_arrow(location = "bl",which_north = "true",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  pad_x = unit(0.02, "cm"), pad_y = unit(0.8, "cm"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  style = north_arrow_nautical,width = unit(3, "cm"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  height = unit(3.5, "cm")) +   geom_spatvector_text(aes(label = sites$numbers), color = "white",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   fontface = "bold", nudge_x = 50, nudge_y = 400, size = 9) + theme(axis.title = element_blank())
one

##make inset map
hja_site <- subset(hja_sites, site == "RS08")
hja_site$site_name <- c("HJA Forest")
hja_site <- usmap_transform(hja_site)
inset <- (plot_usmap(include = c("OR", "WA"), labels = TRUE) +theme(panel.background = element_rect(color = "black")) + geom_sf(data = hja_site, color = "#40004B", alpha = 0.7, size = 7) +
            geom_label(data = hja_site, aes(label = site_name, geometry = geometry), size = 6, alpha = 0.8, stat = "sf_coordinates", nudge_y = 90000, nudge_x = 50000))
inset
map_final <-
  ggdraw() +
  draw_plot(one) +
  draw_plot(inset, x = 0.05, y = .55, width = .35, height = .4)
map_final

ggsave("map1.png", map_final, units = "in", width = 15, height = 10, dpi = 600)

