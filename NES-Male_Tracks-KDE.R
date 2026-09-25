##### Load Libraries
# Programming Language
library(tidyverse)

# Analytical and Spatial Toolboxes
library(adehabitatHR)   # KDE for space use (UD)
library(MASS)           # kde2d for weighted KDE
library(KernSmooth)
library(sf)
library(raster)
library(marmap)
library(tidyterra)
library(terra)

# Plotting and Mapping
library(NatParksPalettes)
library(rnaturalearth)  # basemap
library(viridis)

##### identify files and load ######
tracks_data<-read_csv('Data/All_Tracks.csv')

# IDs to exclude - too short, low quality
skip_ids <- c("2006064","2006065","2007021","2008051","2008053","2008058",
              "2010042","2010084","2010085","2011044","2011045","2011046",
              "2011047","2011048","2015021","2015022","2015058","2026013")

# Format data
tracks_data <- tracks_data |>
  rename(TOPPID = id) |>
  mutate(TOPPID = as.factor(TOPPID),
         lon360 = (lon + 360) %% 360) |>
  group_by(TOPPID) |>
  mutate(Year = as.factor(year(first(date))),
         DayofTrip = ceiling(as.numeric(difftime(date, first(date), units = "days"))))  |>
  ungroup() |>
  filter(!(TOPPID %in% skip_ids)) |>
  droplevels()


####### Tracks KDE ########
# --- Convert to SpatialPointsDataFrame ---
tracks_sf<-  tracks_data |>
  st_as_sf(coords = c("lon360", "lat"), crs = 4326) |>
  st_transform(crs = "+proj=laea +lon_0=-170 +lat_0=55 +units=m")

tracks_proj <- tracks_sf |>
  as("Spatial")
tracks_proj$group <- "all"

# --- Group-level bandwidth (all animals pooled) ---
all_coords <- coordinates(tracks_proj)
hx_all <- dpik(all_coords[, 1])
hy_all <- dpik(all_coords[, 2])
h_group <- sqrt(hx_all * hy_all)

kde_combined <- kernelUD(tracks_proj[, "group"], h = h_group, grid = 500, extent = 0.25)

kde_combined_df <- as.data.frame(
  raster(as(kde_combined$all, "SpatialPixelsDataFrame")), xy = TRUE) %>%
  rename(density = 3) %>%
  filter(!is.na(density))


####### Download Bathymetry #########
# Contours based on ETOPO2
bathyetopo<- marmap::getNOAA.bathy(lon1 = 170, lon2 = -120, lat1 = 36, lat2 = 65, 
                                   resolution = 0.25, keep=T, antimeridian = T)

# Define 200 and 1000m contours
etopoBreaks<-bathyetopo |>
  marmap::as.raster()|> 
  rast()|>
  as.contour(levels=c(-200,-1000)) |>
  tidyterra::as_sf() |>
  st_transform(crs = "+proj=laea +lon_0=-170 +lat_0=55 +units=m")

# Create lines
tracks_lines <- tracks_sf |>
  arrange(TOPPID, date) |>  # Arrange before summarizing
  group_by(TOPPID) |> #Without this the points are connected in a weird order this keeps it so that st_cast mess up the order of the points
  summarize(do_union = FALSE, .groups = "drop") |>
  st_cast("LINESTRING")


####### Plotting #######
# Basemap
world <- ne_countries(scale = "medium", returnclass = "sf") |>
  # match projection
  st_transform(crs = st_crs("+proj=laea +lon_0=-170 +lat_0=55 +units=m"))  

world_clipped <-  world |>
  st_make_valid() |>
  st_crop(world, kde_combined_df)

# Build the graticule matching your plot extent
grat <- sf::st_graticule(
  lat = seq(floor(min(kde_combined_df$y)), ceiling(max(kde_combined_df$y)), by = 5),
  lon = seq(floor(min(kde_combined_df$x)), ceiling(max(kde_combined_df$x)), by = 20))

# Bounding boxes for three animals
boxes <- data.frame(
  name = c("G841", "J914", "J916"),
  xmin = c(178, 185, 205),
  xmax = c(195, 195, 210),
  ymin = c(51, 51, 55),
  ymax = c(52.5, 53, 57.5))

# Helper function to convert one row into an sf polygon
make_bbox_sf <- function(xmin, ymin, xmax, ymax, name, crs = 4326) {
  st_as_sfc(st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
                    crs = crs)) |>
    st_sf(name = name, geometry = _)}

# Build the sf object from the data frame
bboxes <- pmap(boxes, make_bbox_sf) |>
  bind_rows() 
#  sf::st_segmentize(dfMaxLength = 0.1)   # add a vertex every 0.1 degrees

colors <- natparks.pals("Cuyahoga",3)
seal_colors <- colors[c(3,1,2)]

colors <- natparks.pals("Torres",7)
seal_colors <- colors[c(2,3,4)]

# Main Map
map<-ggplot() +
  geom_raster(data = kde_combined_df, aes(x = x, y = y, fill = density), 
              interpolate = TRUE) +
  scale_fill_viridis(option = "viridis", name = "Relative Density",
                     na.value = "transparent",
                     guide = guide_colorbar(label = FALSE)) +
  geom_sf(data = world_clipped, fill = "gray40", color = "gray60", 
          linewidth = 0.2) +
  geom_sf(data = etopoBreaks, aes(linetype = as.factor(level)), 
          color = "white", linewidth = 0.2, alpha = 0.8) +
  geom_sf(data = grat, color = "grey70", linewidth = 0.2, 
          linetype = "dashed", alpha = 0.5) +
  scale_linetype_manual(values = c('solid', "dotted"), 
                        labels = c("1000m", "200m"),
                        name = "Depth") +
  geom_sf(data = tracks_lines, color = "slateblue" ,linewidth = 0.2, alpha = 0.4)+
  geom_sf(data = bboxes, fill = NA, aes(color = name), linewidth = 0.6) +
  scale_color_manual("Seal ID", values=seal_colors) +
  guides(color = "none") +
  coord_sf(xlim = c(min(kde_combined_df$x), max(kde_combined_df$x)),
           ylim = c(min(kde_combined_df$y), max(kde_combined_df$y)),
           expand = FALSE,
           datum = sf::st_crs(4326)) +
  scale_x_continuous(breaks = seq(160, 240, by = 20), position = "bottom") +
  scale_y_continuous(breaks = seq(30, 65, by = 10)) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.085),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "bottom",   # aligns legends to bottom of the box
        legend.background = element_rect(fill = alpha("black", 0.5), color = NA),
        legend.text = element_text(color = "white", size = 11),
        legend.title = element_text(color = "white", size = 12),
        legend.key.size = unit(0.4, "cm"),
        panel.spacing = unit(0.5, "lines"),
        axis.title = element_blank(),
        axis.text = element_text(size = 14))
map
ggsave("Figures/Track_KDE_HiRes.png", 
       width = 14, height = 6, units = "in", bg = "white", dpi = 600)

