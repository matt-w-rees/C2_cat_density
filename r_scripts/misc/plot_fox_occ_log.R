
library(dplyr)
library(mgcv)
library(sp)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(terra)


# LOAD  -------------------------------------------------------------
# load fox presence absence records
records <- read.csv("raw_data/spp_records_pa.csv")
# split by region
records_glenelg <- filter(records, region == "glenelg")
records_otways <- filter(records, region == "otways")

# lower glenelg shapefiles
LGNPN_traps <- read.table("raw_data/LGNPN_traps.txt")
LGNPS_traps <- read.table("raw_data/LGNPS_traps.txt")
LGNP_traps <- rbind(LGNPN_traps, LGNPS_traps)
names(LGNP_traps) <- c("station", "x", "y")
head(LGNP_traps)

# dissolved buffer zone
LGNPN_dissolved <- st_read(dsn='raw_data/lower_glenelg_shp/LGNP_dissolved_32754.gpkg') 

rivers <- st_read(dsn='raw_data/rivers/rivers_buffer_glenelg_32754.gpkg') 
glenelg_river <- filter(rivers, NAME == "GLENELG RIVER")
glenelg_river <- st_union(glenelg_river)
glenelg_river <- st_as_sf(glenelg_river)
plot(glenelg_river)

# load models
gam_g_fox <- readRDS("models/gam_g_fox.RData")
gam_o_fox <- readRDS("models/gam_o_fox.RData")



# GLENELG FOX PLOT -----------------------------------------------------------
## STEP 1) Make buffer zone around cameras to restrict plotting - optional 
#take just the cams
records_glenelg_cams <- distinct(records_glenelg, station, .keep_all = TRUE)
records_glenelg_cams <- subset(records_glenelg_cams, select = c(station, x, y))
# add lower glenelg cams - for plotting
records_glenelg_cams_lg <- rbind(records_glenelg_cams, LGNP_traps)
records_glenelg_cams_lg_sf <- st_as_sf(records_glenelg_cams_lg, coords = c("x", "y"), crs = 32754) 
plot(records_glenelg_cams_lg_sf)

# change it to sf class (not lg)
records_glenelg_cams_sf <- st_as_sf(records_glenelg_cams, coords = c("x", "y"), crs = 32754) 
plot(records_glenelg_cams_sf)

# make a 4km dissolved buffer around each camera
glenelg_cams_buffer <- st_as_sf(st_union(st_buffer(records_glenelg_cams_sf, 4000)))
plot(glenelg_cams_buffer)

# add in lower glenelg buffers
glenelg_cams_buffer_union <- st_union(glenelg_cams_buffer, LGNPN_dissolved)
plot(glenelg_cams_buffer_union)

# convert to SpatVector
glenelg_cams_buffer_union_spv <- terra::vect(glenelg_cams_buffer_union)
plot(glenelg_cams_buffer_union_spv)

# extract coordinates to dataframe
glenelg_buffer_df <- as.data.frame(geom(glenelg_cams_buffer_union_spv))
head(glenelg_buffer_df)
summary(glenelg_buffer_df)

# split by each group (seperate polygon for each grid)
buffer_df1 <- glenelg_buffer_df[which(glenelg_buffer_df$part == "1"),]
buffer_df2 <- glenelg_buffer_df[which(glenelg_buffer_df$part == "2"),]
buffer_df3 <- glenelg_buffer_df[which(glenelg_buffer_df$part == "3"),]
buffer_df4 <- glenelg_buffer_df[which(glenelg_buffer_df$part == "4"),]
buffer_df5 <- glenelg_buffer_df[which(glenelg_buffer_df$part == "5"),]

## STEP 2) MAKE A GRID TO PREDICT DATA INTO
data_glenelg_plot = expand.grid(
  x = seq(min(glenelg_buffer_df$x), 
          max(glenelg_buffer_df$x),
          length=500),
  y = seq(min(glenelg_buffer_df$y),
          max(glenelg_buffer_df$y),
          length=500),
  survey_duration = 60)

# subset to just locations within  of each camera 
data_glenelg_plot1 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df1, x, y)),]
data_glenelg_plot2 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df2, x, y)),]
data_glenelg_plot3 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df3, x, y)),]
data_glenelg_plot4 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df4, x, y)),]
data_glenelg_plot5 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df5, x, y)),]
data_glenelg_plot <- rbind(data_glenelg_plot1, data_glenelg_plot2, data_glenelg_plot3, data_glenelg_plot4, data_glenelg_plot5)

# predict model results into dataframe
data_glenelg_plot <- cbind(data_glenelg_plot, predict.gam(gam_g_fox, newdata = data_glenelg_plot, se.fit = TRUE, type = "link", exclude = c("s(survey_duration)")))
data_glenelg_plot <- rename(data_glenelg_plot, fox_predicted = fit,  fox_predicted_se = se.fit) # rename
# make spatial for plotting
data_glenelg_plot_sf <- st_as_sf(data_glenelg_plot, coords = c("x","y"), crs = 32754)

# plot
g_fox_plot <- ggplot() + 
  geom_sf(data = data_glenelg_plot_sf, aes(colour = fox_predicted), size = 3) + 
  scale_colour_viridis("log(fox occurrence)", option = "viridis", limits = c(-2.2, 0.91)) + 
  geom_sf(data = glenelg_river, size = 1.2) +
  geom_sf(data = records_glenelg_cams_lg_sf, colour = "white", fill = "white", alpha = 1, size = 1, shape = 16) +
  theme_bw(10) + 
  theme(axis.title = element_blank(),
        legend.position = "none") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 8500, y = mean(data_glenelg_plot$y) + 24000, label = "Hotpsur (R2; NI)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) + 23000, y = mean(data_glenelg_plot$y) + 13000, label = "Annya (R1; NI)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 14500, y = mean(data_glenelg_plot$y) - 14500, label = "Cobboboonee (R1; I)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) + 9500, y = mean(data_glenelg_plot$y) - 23000, label = "Mt Clay (R2; I)") +
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 24500, y = mean(data_glenelg_plot$y) + 13000, label = "LGNP north (R3; NI)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 33000, y = mean(data_glenelg_plot$y) - 8000, label = "LGNP south (R3; I)")
g_fox_plot

# plot uncertainty
# plot
g_fox_plot_se <- ggplot() + 
  geom_sf(data = data_glenelg_plot_sf, aes(colour = fox_predicted_se), size = 3) + 
  scale_colour_viridis("Standard error", option = "viridis") + 
  geom_sf(data = glenelg_river, size = 1.2) +
  geom_sf(data = records_glenelg_cams_lg_sf, colour = "white", fill = "white", alpha = 1, size = 1, shape = 16) +
  theme_bw(10) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 8500, y = mean(data_glenelg_plot$y) + 24000, label = "Hotpsur (R2; NI)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) + 23000, y = mean(data_glenelg_plot$y) + 13000, label = "Annya (R1; NI)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 14500, y = mean(data_glenelg_plot$y) - 14500, label = "Cobboboonee (R1; I)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) + 9500, y = mean(data_glenelg_plot$y) - 23000, label = "Mt Clay (R2; I)") +
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 24500, y = mean(data_glenelg_plot$y) + 13000, label = "LGNP north (R3; NI)") + 
  annotate("text", size = 3.5, x = mean(data_glenelg_plot$x) - 33000, y = mean(data_glenelg_plot$y) - 8000, label = "LGNP south (R3; I)")
g_fox_plot_se


# OTWAY PLOTS -----------------------------------------------------------
## STEP 1) Make buffer zone around cameras to restrict plotting - optional 
# change it to sf class
records_otways_cams <- st_as_sf(records_otways, coords = c("x", "y"), crs = 32754) 
# make a 4km buffer around each camera
otways_cams_buffer = st_buffer(records_otways_cams, 4000)
# dissolve the buffer
otways_cams_buffer = st_union(otways_cams_buffer)
# convert back to dataframe
otways_buffer_df <- fortify(as_Spatial(otways_cams_buffer))%>%
  transmute(x = long, y = lat, order = order, group = group)


## STEP 2) MAKE A GRID TO PREDICT DATA INTO
data_otways_plot = expand.grid(
  x = seq(min(otways_buffer_df$x), 
          max(otways_buffer_df$x),
          length=100),
  y = seq(min(otways_buffer_df$y),
          max(otways_buffer_df$y),
          length=100),
  station = "T053", 
  station_year = "T053x2017",
  year = seq(2017,2019,by=1),
  survey_duration = 60)

# subset to just locations within buffer zone
data_otways_plot = data_otways_plot[with(data_otways_plot, inSide(otways_buffer_df, x, y)),]
# predict
data_otways_plot <- cbind(data_otways_plot, predict.gam(gam_o_fox, newdata = data_otways_plot, se.fit = TRUE, type = "link", exclude = c("s(station)", "s(survey_duration)")))
data_otways_plot <- rename(data_otways_plot, fox_predicted = fit,  fox_predicted_se = se.fit) # rename
# make spatial for plotting
data_otways_plot_sf <- st_as_sf(data_otways_plot, coords = c("x","y"), crs = 32754)

## STEP 3) PLOT
# fox plot
o_fox_plot <- ggplot() + 
  geom_sf(data = data_otways_plot_sf, aes(colour = fox_predicted), size = 3) + 
  scale_colour_viridis("log(fox occurrence)", option = "viridis", limits = c(-2.2, 0.91)) + 
  facet_wrap(~year, nrow = 1) +
  geom_sf(data = records_otways_cams, colour = "white", fill = "white", alpha = 1, size = 1, shape = 16) +
  theme_bw(10) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom") + 
  annotate("text", x = mean(data_otways_plot$x) + 7500, y = mean(data_otways_plot$y) + 7500, label = "NI") + 
  annotate("text", x = mean(data_otways_plot$x) - 7500, y = mean(data_otways_plot$y) - 7500, label = "I") 
o_fox_plot

# fox plot - uncertainty
o_fox_plot_se <- ggplot() + 
  geom_sf(data = data_otways_plot_sf, aes(colour = fox_predicted_se), size = 3) + 
  scale_colour_viridis("Standard error", option = "viridis") + 
  facet_wrap(~year, nrow = 1) +
  geom_sf(data = records_otways_cams, colour = "white", fill = "white", alpha = 1, size = 1, shape = 16) +
  theme_bw(10) + 
  theme(axis.title = element_blank(),
        legend.position = "bottom") + 
  annotate("text", x = mean(data_otways_plot$x) + 7500, y = mean(data_otways_plot$y) + 7500, label = "NI") + 
  annotate("text", x = mean(data_otways_plot$x) - 7500, y = mean(data_otways_plot$y) - 7500, label = "I") 
o_fox_plot_se


# COMBINE AND PLOT --------------------------------------------------------

png("figs/fig2A_600dpi.png", width = 8.5, height = 5, res = 600, units = "in")
g_fox_plot + plot_annotation(title = '(a) Glenelg region') 
dev.off()

png("figs/fig2B_600dpi.png", width = 8.5, height = 4.25, res = 600, units = "in")
o_fox_plot + plot_annotation(title = '(b) Otway region') 
dev.off()

# to merge, using imagemagick - type in the terminal:
# convert figs/fig2A_600dpi.png figs/fig2B_600dpi.png -gravity center -append figs/fig2_600dpi.png


# supp uncertainty figs

png("figs/fox_occ_se_glenelg_600dpi.png", width = 8.5, height = 5, res = 600, units = "in")
g_fox_plot_se
dev.off()

png("figs/fox_occ_se_otways_600dpi.png", width = 8.5, height = 4.25, res = 600, units = "in")
o_fox_plot_se
dev.off()

