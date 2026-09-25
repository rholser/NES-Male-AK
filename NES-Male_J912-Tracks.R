##### Load Libraries
# Programming Language
library(tidyverse)

# Analytical and Spatial Toolboxes
library(sf)
library(raster)
library(marmap)
library(tidyterra)
library(terra)

# Plotting and Mapping
library(NatParksPalettes)
library(rnaturalearth)  # basemap
library(viridis)
library(maptiles)
library(patchwork)

##### identify files and load ######
tracks_data<-read_csv('Data/All_Tracks.csv')

keep_ids <- c("2025035","2026013")

tracks_data <- tracks_data %>%
  rename(TOPPID = id) %>%
  mutate(TOPPID = as.factor(TOPPID),
         lon360 = (lon + 360) %% 360) %>%
  group_by(TOPPID) %>%
  mutate(Year = as.factor(year(first(date)))) %>%
  ungroup() %>%
  filter((TOPPID %in% keep_ids)) %>%
  droplevels()

# Convert to SpatialPointsDataFrame
laea <- st_crs("+proj=laea +lon_0=-170 +lat_0=55 +units=m")
prj <- "+proj=laea +lon_0=-170 +lat_0=55 +datum=WGS84 +units=m +no_defs"

tracks_sf <- tracks_data |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(laea)

# Create lines
tracks_lines <- tracks_sf |>
  arrange(TOPPID, date) |>  # Arrange before summarizing
  group_by(TOPPID) |> #Without this the points are connected in a weird order this keeps it so that st_cast mess up the order of the points
  summarize(do_union = FALSE, .groups = "drop") |>
  st_cast("LINESTRING") 

tracks_prj <- tracks_lines |>
  st_segmentize(units::set_units(10, "km", mode = "standard")) |>
  st_transform(prj)


# --- Define zoom bounding boxes --- #
# Zoom window
lon_rng <- c(-172, -167)
lat_rng <- c(51.5, 54)
zoom_box <- st_bbox(c(xmin = min(lon_rng), xmax = max(lon_rng),
                      ymin = min(lat_rng), ymax = max(lat_rng)),
                    crs = 4326) |>
  st_as_sfc() |>
  st_segmentize(units::set_units(5, "km", mode = "standard")) |>
  st_transform(prj)

zb <- st_bbox(zoom_box)

# Define wide bounding box
lon_wide <- c(-178, -120)
lat_wide <- c(34, 65)
wb <- st_bbox(c(xmin = lon_wide[1], xmax = lon_wide[2],
                ymin = lat_wide[1], ymax = lat_wide[2]), crs = 4326) |>
  st_as_sfc() |>
  st_segmentize(units::set_units(25, "km", mode = "standard")) |>
  st_transform(prj) |>
  st_bbox()

track_cols <- natparks.pals("Denali",4)
track_cols <- track_cols[c(1,3)]

# --- Plots --- #
# Zoomed plot
aoi <- st_as_sfc(st_bbox(c(xmin = -173, xmax = -166, ymin = 50, ymax = 55),
                         crs = 4326))
tiles <- get_tiles(aoi, provider = "Esri.OceanBasemap", zoom = 9,
                   crop = TRUE, cachedir = "Data/tiles")

zoom<-ggplot() +
  geom_spatraster_rgb(data = tiles, maxcell = Inf) +
  geom_sf(data = tracks_prj, aes(colour = TOPPID), linewidth = 0.8, alpha = 0.8) +
  scale_colour_manual(
    values = track_cols,
    labels = c("2025035" = "PM 2025",
               "2026013" = "PB 2026"),
    name   = "Trip")+
  coord_sf(crs = prj, 
           xlim =  zb[c("xmin", "xmax")], 
           ylim = zb[c("ymin", "ymax")], 
           expand = FALSE) +
  #labs(caption = get_credit("Esri.OceanBasemap")) +
  theme(legend.position  = "none",
        axis.title       = element_blank(),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.ticks.length = unit(0, "pt"))

ggsave("Figures/J912_Tracks_Zoom_HiRes.png", 
       width = 8, height = 6, units = "in", bg = "white", dpi = 600)

# Wide plot
aoi_w <- st_as_sfc(st_bbox(c(xmin = -176, xmax = -114, ymin = 25, ymax = 64),
                           crs = 4326))
tiles_w <- get_tiles(aoi_w, provider = "Esri.OceanBasemap", zoom = 6,
                     crop = TRUE, cachedir = "Data/tiles")
lims <- st_sfc(
  st_polygon(list(cbind(
    c(-174, -124, -124, -174, -174),   # lon, closed ring
    c(  28,   28,   56,   56,   28)))),    # lat
  crs = 4326) |>
  st_segmentize(units::set_units(20, km)) |>
  st_transform(prj)

wb2 <- st_bbox(lims)

wide<-ggplot() +
  geom_spatraster_rgb(data = tiles_w, maxcell = 5e6) +
  geom_sf(data = tracks_prj, aes(colour = TOPPID), linewidth = 0.8) +
  geom_sf(data = zoom_box, colour = "black", fill = NA, linewidth = 0.8) +
  scale_colour_manual(values = track_cols,
                      labels = c("2025035" = "PM 2025",
                                 "2026013" = "PB 2026"),
                      name   = "Trip")+
  scale_x_continuous(breaks = seq(-174, -120, by = 10)) +
  scale_y_continuous(breaks = seq(32, 60, by = 10)) +
  coord_sf(crs = prj,
           xlim = wb2[c("xmin", "xmax")],
           ylim = wb2[c("ymin", "ymax")],
           expand = FALSE,
           datum = sf::st_crs(4326),
           label_axes = list(bottom = "E", left = "N"))+
  theme_bw(base_size = 10) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = alpha("grey40", 0.4),
                                        linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.title       = element_blank())

ggsave("Figures/J912_Tracks_HiRes.png", 
       width = 8, height = 6, units = "in", bg = "white", dpi = 600)

# Combined plot
combined <- wide + zoom +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

inset <- wide +
  inset_element(zoom + theme(plot.background = element_rect(colour = "black",
                                                            linewidth = 0.6)),
                left = 0, bottom = 0, right = 0.45, top = 0.425)

ggsave("Figures/J912_track_map_inset.png", inset, width = 8, height = 6.5, dpi = 600)
