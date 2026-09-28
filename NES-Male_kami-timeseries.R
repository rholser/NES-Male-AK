# Load libraries
# Programming
library(tidyverse)
library(data.table)
library(geosphere)

# Spatial tools
library(sf)
library(tidyterra)
library(terra)
library(marmap)
library(rnaturalearth)  # basemap
library(rnaturalearthdata)

# Plotting tools
library(ggspatial)
library(NatParksPalettes)
library(viridis)
library(ggside)
library(hexbin)
library(scico)
library(ggnewscale)

#####------identify files and load--------##############
# load and format
combined_data <- read_csv("Data/All_Kami_data.csv")

combined_data <- combined_data |>
  mutate(TOPPID = as.factor(TOPPID),
         Time = as.POSIXct(Time, format = '%d-%b-%Y %H:%M:%S'))

# load and format
combined_divestat <- read_csv("Data/All_Kami_DiveStat.csv")

combined_divestat <- combined_divestat |>
  mutate(TOPPID = as.factor(TOPPID),
         DateTime = as.POSIXct(DateTime, format = '%d-%b-%Y %H:%M:%S')) |>
  group_by(TOPPID) |>
  mutate(DayOfTrip = JulDate - first(JulDate)) |>
  arrange(DateTime, .by_group = TRUE) |>
  # calculate dive-to-dive transit rate:
  mutate(dist_m = c(NA, distHaversine(
    cbind(Lon[-n()], Lat[-n()]),
    cbind(Lon[-1],   Lat[-1]))),
    dt_sec = c(NA, diff(as.numeric(DateTime))),
    speed_m_s = dist_m / dt_sec) |>
  ungroup() |>
  # convert transit rate to km/hr
  mutate(speed_km_hr = speed_m_s*(60*60)/1000,
         Lon360 = case_when(Lon < 0 ~ Lon + 360,
                            TRUE ~ Lon))

#####------Download Bathymetry Data--------###############
options(timeout = 600)

# GEBCO 15 arc-second bathymetry data
bathyetopo<- marmap::getNOAA.bathy(lon1 = 170, lon2 = -120, lat1 = 36, lat2 = 65, 
                                   resolution = 0.25, keep=TRUE, antimeridian = TRUE)

# Extract bathymetry for each dive
combined_divestat$Bathy <- marmap::get.depth(
  bathyetopo, x = combined_divestat[, c("Lon360", "Lat")],  # must be Lon, Lat order
  locator = FALSE)$depth

# Change bathymetry to positive depth values
combined_divestat$Bathy <- -(combined_divestat$Bathy)

# Calculate difference between max dive depth and bathymetry
combined_divestat$ExcessDepth<-combined_divestat$Maxdepth - combined_divestat$Bathy

# Create adjusted baythmetry. If dive depth > bathy, use max depth as bathymetry
combined_divestat <- combined_divestat |>
  mutate(Bathy_adj = case_when(ExcessDepth > 0 ~ Maxdepth,
                               TRUE ~ Bathy),
         DayOfTripR = round(DayOfTrip),
         ToD = case_when( SolarEl >= -6 ~ "Day",
                          SolarEl < -6 ~ "Night")) %>%
  mutate(ToD = as.factor(ToD))

#####-------Map DiveStat metrics onto kamikami timeseries----------########

# Convert both data sets to data.table
setDT(combined_data)
setDT(combined_divestat)

# Set keys for join — must key on TOPPID + datetime to join within individual
setkey(combined_data, TOPPID, Time)
setkey(combined_divestat, TOPPID, DateTime)

# Forward-fill solar elevation, bathymetry, etc. per individual
combined_kami <- combined_divestat[, .(TOPPID, DateTime, SolarEl, DiveNumber, PDI, DayOfTrip, speed_km_hr, Bathy, Bathy_adj, Lat, Lon360)][
  combined_data, roll = Inf, on = .(TOPPID, DateTime = Time)]

result <- combined_kami |>
  filter(CorrectedDepth>25) |>
  filter(!is.na(SolarEl), KAMI_L>0) |>
  mutate(ToD = case_when( SolarEl >= -6 ~ "Day",
                          SolarEl < -6 ~ "Night")) |>
  mutate(ToD = as.factor(ToD))

# Day of trip statistics
DayofTrip <- combined_divestat |>
  group_by(TOPPID, DayOfTripR) |>
  summarize(KamiEvents = sum(KamiEventSumDive),
            Lat = mean(Lat, na.rm = TRUE),
            Lon360 = mean(Lon360, na.rm = TRUE),
            Speed = mean(speed_km_hr, na.rm = TRUE),
            .groups = "drop")

# Proportion of dives with kami events  
Kami_prop <- combined_divestat |>
  group_by(TOPPID) |>
  summarize(nDives      = sum(!is.na(KamiNoSumDive100m)),
            KamiDives   = sum(KamiNoSumDive100m > 0, na.rm = TRUE),
            noKamiDives = nDives - KamiDives,
            propKami    = KamiDives / nDives,
            .groups = "drop")

#######-------------Bathy-Kami Plots----------##################
shade_poly <- data.frame(
  x = c(0, 1250, 0),
  y = c(0, 1250, 1250))

gradient_df <- data.frame(
  x = seq(-65, 65, length.out = 500)) |>
  dplyr::mutate(shade = dplyr::case_when(
    x < -18 ~ 0,           # astronomical night
    x < 0   ~ (x + 18)/18, # twilight ramp 0→1
    TRUE    ~ 1))            # day

colors<-scico(n = 4, palette = 'nuuk', categorical = TRUE)
ToD_colors<-colors[c(2,1)]

################-------J916-2025037-------##################
#Subset data and summarize
J916_plot<-result |>
  filter(TOPPID=="2025037") |>
  filter(Lon360<=210)

J916_summary <- J916_plot |>
  group_by(Hour) |>
  summarise(n = n(),                              # foraging event count per hour
            mean_SolarEl = mean(SolarEl, na.rm = TRUE),
            .groups = "drop")

# Identify dark hours and find their boundaries
dark_hours <- J916_summary |>
  filter(mean_SolarEl < -6) |>
  pull(Hour)

# Build segments at the edges where dark transitions to light (or vice versa)
# An edge exists between hour h and h+1 if exactly one of them is "dark"
all_hours <- 0:23
edge_positions <- sapply(all_hours, function(h) {
  this_dark <- h %in% dark_hours
  next_dark <- ((h + 1) %% 24) %in% dark_hours
  if (this_dark != next_dark) (h + 0.5) else NA_real_}) |> 
  na.omit()

# y-range for the segments — match your bar heights
y_max <- max(J916_summary$n)

####Plots
# ----- Figure 2 J916 Subplot ----- #
ggplot(data = J916_plot) +
  geom_polygon(data = shade_poly, aes(x = x, y = y),
               fill = "grey40", alpha = 0.5) +
  geom_hex(aes(x=Bathy_adj,y=CorrectedDepth), bins = 100)+
  scale_fill_scico("Feeding\nEvents", palette = 'lapaz') +
  # new_scale_fill() +                           # reset fill
  # geom_ysidedensity(aes(y = CorrectedDepth, group=ToD, fill = ToD), alpha = 1, linewidth = 1, position = "stack") +  # marginal density on y-axis
  # scale_fill_manual("Time of Day", values=ToD_colors) +
  xlim(2000,0)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Foraging Depth (m)", x = "Sea Floor Depth (m)", title = "J916 Feeding Events") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank(),
        legend.position = "")

ggsave('Figures/J916_Kami-Bathymetry_adj.png', width=10, height=7, dpi=600)

# ----- Figure 5 J916 Subplots ----- #
ggplot(data = J916_plot) +
  # Background gradient
  geom_rect(data = gradient_df,
            aes(xmin = x, xmax = dplyr::lead(x), 
                ymin = -Inf, ymax = Inf, fill = x),
            show.legend = FALSE) +
  scico::scale_fill_scico(palette = "nuuk", guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_hex(aes(x=SolarEl,y=CorrectedDepth), bins = 50)+
  scale_fill_scico("Feeding\nEvents", palette = 'lapaz') +
  xlim(-65,65)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Foraging Depth (m)", x = "Solar Elevataion", title = "J916 Feeding Events") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank(),
        legend.position = c(0.92, 0.5),
        legend.justification = c("right","center"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

ggsave('Figures/J916_Kami-SolarEl.png', width=8, height=6, dpi=600)

ggplot(J916_summary, aes(x = Hour, y = n, fill = mean_SolarEl)) +
  geom_col(width = 1, color = NA) +
  geom_segment(data = data.frame(x = edge_positions),
               aes(x = x, xend = x, y = 0, yend = y_max * 1.05),
               inherit.aes = FALSE,
               color = "grey20", linewidth = 1.5) +
  coord_polar(start = 0) +
  scale_x_continuous(breaks = 0:23, limits = c(-0.5, 23.5)) +
  scale_fill_scico(name = "Solar\nElevation",
                   palette = "nuuk",
                   limits = c(-65, 65)) +
  ggthemes::theme_few() +
  labs(x = "Hour of the Day (UTC)", y = "") +
  theme(panel.grid.major = element_line(size = 0.3, colour = "grey", linetype = "dashed"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

ggsave('Figures/J916_Kami-Hour-SolarEl.png', width=7, height=6, dpi=600)

################-----G841-2025038-----##################
#Subset data and summarize
G841_plot<-result |>
  filter(TOPPID=="2025038") |>
  filter(Lon360<=195 & Lat >51)

G841_summary <- G841_plot |>
  group_by(Hour) |>
  summarise(n = n(),                              # foraging event count per hour
            mean_SolarEl = mean(SolarEl, na.rm = TRUE),
            .groups = "drop")

# Identify dark hours and find their boundaries
dark_hours <- G841_summary |>
  filter(mean_SolarEl < -6) |>
  pull(Hour)

# Build segments at the edges where dark transitions to light (or vice versa)
# An edge exists between hour h and h+1 if exactly one of them is "dark"
all_hours <- 0:23
edge_positions <- sapply(all_hours, function(h) {
  this_dark <- h %in% dark_hours
  next_dark <- ((h + 1) %% 24) %in% dark_hours
  if (this_dark != next_dark) (h + 0.5) else NA_real_
}) |> na.omit()

# y-range for the segments — match your bar heights
y_max <- max(G841_summary$n)

####Plots
# ----- Figure 2 G841 Subplot ----- #
ggplot(data = G841_plot) +
  geom_polygon(data = shade_poly, aes(x = x, y = y),
               fill = "grey40", alpha = 0.5) +
  geom_hex(aes(x=Bathy_adj,y=CorrectedDepth), bins = 100)+
  scale_fill_scico("Feeding\nEvents", palette = 'lapaz') +
  xlim(2000,0)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Foraging Depth (m)", x = "Sea Floor Depth (m)", title = "G841 Feeding Events") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank(),
        legend.position = "")

ggsave('Figures/G841_Kami-Bathymetry_adj.png', width=10, height=7, dpi=600)

# ----- Figure 5 G841 Subplots ----- #
ggplot(data = G841_plot) +
  # Background gradient
  geom_rect(data = gradient_df,
            aes(xmin = x, xmax = dplyr::lead(x), 
                ymin = -Inf, ymax = Inf, fill = x),
            show.legend = FALSE) +
  scico::scale_fill_scico(palette = "nuuk", guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_hex(aes(x=SolarEl,y=CorrectedDepth), bins = 50)+
  scale_fill_scico("Feeding\nEvents", palette = 'lapaz') +
  xlim(-65,65)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Foraging Depth (m)", x = "Solar Elevataion", title = "G841 Feeding Events") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank(),
        legend.position = c(0.92, 0.5),
        legend.justification = c("right","center"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

ggsave('Figures/G841_Kami-SolarEl.png', width=8, height=6, dpi=600)

ggplot(G841_summary, aes(x = Hour, y = n, fill = mean_SolarEl)) +
  geom_col(width = 1, color = NA) +
  geom_segment(data = data.frame(x = edge_positions),
               aes(x = x, xend = x, y = 0, yend = y_max * 1.05),
               inherit.aes = FALSE,
               color = "grey20", linewidth = 1.5) +
  coord_polar(start = 0) +
  scale_x_continuous(breaks = 0:23, limits = c(-0.5, 23.5)) +
  scale_fill_scico(name = "Solar\nElevation",
                   palette = "nuuk",
                   limits = c(-65, 65)) +
  ggthemes::theme_few() +
  labs(x = "Hour of the Day (UTC)", y = "") +
  theme(panel.grid.major = element_line(size = 0.3, colour = "grey", linetype = "dashed"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

ggsave('Figures/G841_Kami-Hour-SolarEl.png', width=7, height=6, dpi=600)

# ----- Figure S1 ----- #
plot1 <- combined_kami |>
  filter(TOPPID == "2025038") |>
  filter(DateTime > as.POSIXct("2025-11-05") & DateTime < as.POSIXct("2025-11-10"))


ggplot(data=subset(plot1, DateTime < as.POSIXct("2025-11-08")), aes(x=DateTime, y=CorrectedDepth)) + 
  geom_line(aes(color = SolarEl), alpha = 0.7)+
  geom_line(aes(x = DateTime, y = Bathy), color = "black")+
  #geom_point(shape=16, size = 0.6, alpha = 1) +
  geom_point(data = subset(plot1, KAMI_L>0 & DateTime < as.POSIXct("2025-11-08")), 
             color = "darkred", shape = 16, size = 1.5)+
  ggthemes::theme_few()+
  ylim(600, 0)+
  scale_color_viridis(option = "magma", name = "KAMI_L",
                      na.value = "transparent") +
  labs( y = "Depth (m)", x = "Time", title = "G841 Benthic Kami Events")+
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", 
                                        linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

ggsave('Figures/G841_benthic_timeseries2_v2.png',width=18, height = 9, dpi = 300)


################---J914-2025036---##################
#Subset data and summarize
J914_plot<-result |>
  filter(TOPPID=="2025036") |>
  filter(Lon360<=195)

J914_summary <- J914_plot |>
  group_by(Hour) |>
  summarise(n = n(),                              # foraging event count per hour
            mean_SolarEl = mean(SolarEl, na.rm = TRUE),
            .groups = "drop")

# Identify dark hours and find their boundaries
dark_hours <- J914_summary |>
  filter(mean_SolarEl < -6) |>
  pull(Hour)

# Build segments at the edges where dark transitions to light (or vice versa)
# An edge exists between hour h and h+1 if exactly one of them is "dark"
all_hours <- 0:23
edge_positions <- sapply(all_hours, function(h) {
  this_dark <- h %in% dark_hours
  next_dark <- ((h + 1) %% 24) %in% dark_hours
  if (this_dark != next_dark) (h + 0.5) else NA_real_}) |> na.omit()

# y-range for the segments — match your bar heights
y_max <- max(J914_summary$n)

####Plots
# ----- Figure 2 J914 Subplot ----- #
ggplot(data = J914_plot) +
  geom_polygon(data = shade_poly, aes(x = x, y = y),
               fill = "grey40", alpha = 0.5) +
  geom_hex(aes(x=Bathy_adj,y=CorrectedDepth), bins = 100)+
  scale_fill_scico("Feeding\nEvents", palette = 'lapaz') +
  xlim(2000,0)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Foraging Depth (m)", x = "Sea Floor Depth (m)", title = "J914 Feeding Events") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank(),
        legend.position = "")

ggsave('Figures/J914_Kami-Bathymetry_adj.png', width=10, height=7, dpi=600)

# ----- Figure 5 J914 Subplots ----- #
ggplot(data = J914_plot) +
  # Background gradient
  geom_rect(data = gradient_df,
            aes(xmin = x, xmax = dplyr::lead(x), 
                ymin = -Inf, ymax = Inf, fill = x),
            show.legend = FALSE) +
  scico::scale_fill_scico(palette = "nuuk", guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_hex(aes(x=SolarEl,y=CorrectedDepth), bins = 50)+
  scale_fill_scico("Feeding\nEvents", palette = 'lapaz') +
  xlim(-65,65)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Foraging Depth (m)", x = "Solar Elevataion", title = "J914 Feeding Events") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank(),
        legend.position = c(0.92, 0.5),
        legend.justification = c("right","center"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

ggsave('Figures/J914_Kami-SolarEl.png', width=8, height=6, dpi=600)

ggplot(J914_summary, aes(x = Hour, y = n, fill = mean_SolarEl)) +
  geom_col(width = 1, color = NA) +
  geom_segment(data = data.frame(x = edge_positions),
               aes(x = x, xend = x, y = 0, yend = y_max * 1.05),
               inherit.aes = FALSE,
               color = "grey20", linewidth = 1.5) +
  coord_polar(start = 0) +
  scale_x_continuous(breaks = 0:23, limits = c(-0.5, 23.5)) +
  scale_fill_scico(name = "Solar\nElevation",
                   palette = "nuuk",
                   limits = c(-65, 65)) +
  ggthemes::theme_few() +
  labs(x = "Hour of the Day (UTC)", y = "") +
  theme(panel.grid.major = element_line(size = 0.3, colour = "grey", linetype = "dashed"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

ggsave('Figures/J914_Kami-Hour-SolarEl.png', width=7, height=6, dpi=600)

#############---All Seals, all regions---############

# ----- Figure 4 ----- #
ggplot(data = result) +
  geom_polygon(data = shade_poly, aes(x = x, y = y),
               fill = "grey40", alpha = 0.5) +
  geom_hex(aes(x=Bathy_adj,y=CorrectedDepth), bins = 100)+
  scale_fill_scico("Feeding Events", palette = 'lapaz') +
  new_scale_fill() +                           # reset fill
  geom_ysidedensity(aes(y = CorrectedDepth, group=ToD, fill = ToD), alpha = 1, linewidth = 1, position = "stack") +  # marginal density on y-axis
  scale_fill_manual("Time of Day", values=ToD_colors) +
  xlim(6000,0)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Depth (m)", x = "Sea Floor Depth (m)") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank())

ggsave('Figures/All_Kami-Bathymetry_adj.png', width=12, height=6, dpi=300)

# ----- Figure S2 Subplots ----- #
ggplot(data = subset(combined_divestat, KamiNoSumDive == 0)) +                        # reset fill
  geom_density(aes(y = Maxdepth, group=ToD, fill = ToD), alpha = 1, linewidth = 1, position = "stack") +  # marginal density on y-axis
  scale_fill_manual("Time of Day", values=ToD_colors) +
  #xlim(6000,0)+
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Depth (m)", x = "Density", title = "Dives without Jaw Motion") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank())

ggsave('Figures/All_noKami_Density.png', width=4, height=6, dpi=300)

ggplot(data = subset(combined_divestat, KamiNoSumDive > 0)) +
  geom_density(aes(y = Maxdepth, group=ToD, fill = ToD), alpha = 1, linewidth = 1, position = "stack") +  # marginal density on y-axis
  scale_fill_manual("Time of Day", values=ToD_colors) +
  ylim(1250,0)+
  ggthemes::theme_few()+  
  labs(y = "Depth (m)", x = "Density", title = "Dives with Jaw Motion") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank())

ggsave('Figures/All_Kami_Density.png', width=4, height=6, dpi=300)

#######---Day of Trip ~ Kami Kami Plot---###########
#Start and end of transit defined as daily movement speed above/below 3km/hr
transit <- data.frame(SealID    = c("J914", "J916", "G841"),
                      TOPPID    = c("2025036", "2025037", "2025038"),
                      #out_end   = c(45, 30, 49), # end from speed
                      out_end   = c(35, 30, 31), # end from speed + latitude
                      back_start = c(88, 86, 83))    # start of return transit

# Define colors to match map
colors <- natparks.pals("Torres",7)
seal_colors <- colors[c(3,4,2)]

DayofTrip <- DayofTrip |>
  left_join(transit |> dplyr::select(TOPPID, SealID), by = "TOPPID")

# Remove one unusually high kami event day from 2025036 for scaling
test <- DayofTrip |>
  filter(!(TOPPID == "2025036" & KamiEvents >600))

# Calculate cumulative distributions
ecdf_data <- DayofTrip |>
  group_by(SealID, TOPPID) |>
  arrange(DayOfTripR, .by_group = TRUE) |>
  mutate(ecdf_val = cumsum(KamiEvents) / sum(KamiEvents, na.rm = TRUE)) |>
  group_by(SealID) |>
  mutate(ecdf_scaled = ecdf_val * max(KamiEvents, na.rm = TRUE)) |>
  ungroup()

# ----- Figure 3 ----- #
ggplot(data = test, aes(x = DayOfTripR, y=KamiEvents, group = TOPPID, color = TOPPID))+
  geom_rect(data = transit, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = out_end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4) +
  geom_rect(data = transit, inherit.aes = FALSE,
            aes(xmin = back_start, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.4) +
  geom_line(data = ecdf_data, aes(x = DayOfTripR, y = ecdf_scaled), color = "grey30",
            linewidth = 1.5, alpha=0.8) +  
  geom_point(size=2)+
  geom_smooth(method = "gam") +
  scale_color_manual("Seal ID", values = seal_colors) +
  scale_y_continuous(name = "Total Feeding Events",
    sec.axis = sec_axis(~ . , name = "Cumulative Proportion", labels = NULL)) +
  ggthemes::theme_few()+  
  labs(y = "Total Feeding Events", x = "Day of Trip") +
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        axis.title.y.right = element_text(margin = margin(l = 35)),
        ggside.axis.text.x = element_blank(),
        ggside.axis.ticks.x = element_blank(),
        legend.position = "")+
  facet_wrap(~SealID, ncol = 1, scale = "free_y")

ggsave('Figures/All_Kami-DayOfTrip.png', width=9, height=8, dpi=600)

