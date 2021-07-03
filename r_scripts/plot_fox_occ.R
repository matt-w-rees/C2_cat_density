
library(dplyr)
library(mgcv)
library(sp)
library(sf)
library(ggplot2)
library(viridis)


# LOAD  -------------------------------------------------------------
# load fox presence absence records
records <- read.csv("raw_data/spp_records_pa.csv")
# split by region
records_glenelg <- filter(records, region == "glenelg")
records_otways <- filter(records, region == "otways")

# load models
gam_g_fox <- readRDS("models/gam_g_fox.RData")
gam_o_fox <- readRDS("models/gam_o_fox.RData")


# GLENELG FOX PLOT -----------------------------------------------------------
## STEP 1) Make buffer zone around cameras to restrict plotting - optional 
#take just the cams
records_glenelg_cams <- distinct(records_glenelg, station, .keep_all = TRUE)
# change it to sf class
records_glenelg_cams <- st_as_sf(records_glenelg_cams, coords = c("x", "y"), crs = 32754) 
# make a 4km buffer around each camera
glenelg_cams_buffer = st_buffer(records_glenelg_cams, 4000)
# dissolve the buffer
glenelg_cams_buffer = st_union(glenelg_cams_buffer)
# convert back to dataframe
glenelg_buffer_df <- fortify(as_Spatial(glenelg_cams_buffer))%>%
  transmute(x = long, y = lat, order = order, group = group)
# split by each group (seperate polygon for each grid)
buffer_df1 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.1"),]
buffer_df2 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.2"),]
buffer_df3 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.3"),]
buffer_df4 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.4"),]

## STEP 2) MAKE A GRID TO PREDICT DATA INTO
data_glenelg_plot = expand.grid(
  x = seq(min(glenelg_buffer_df$x), 
          max(glenelg_buffer_df$x),
          length=50),
  y = seq(min(glenelg_buffer_df$y),
          max(glenelg_buffer_df$y),
          length=50),
  survey_duration = 60)

# subset to just locations within  of each camera 
data_glenelg_plot1 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df1, x, y)),]
data_glenelg_plot2 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df2, x, y)),]
data_glenelg_plot3 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df3, x, y)),]
data_glenelg_plot4 = data_glenelg_plot[with(data_glenelg_plot, inSide(buffer_df4, x, y)),]
data_glenelg_plot <- rbind(data_glenelg_plot1, data_glenelg_plot2, data_glenelg_plot3, data_glenelg_plot4)

# predict model results into dataframe
data_glenelg_plot <- cbind(data_glenelg_plot, predict.gam(gam_g_fox, newdata = data_glenelg_plot, se.fit = TRUE, type = "link", exclude = c("s(station)", "s(survey_duration)")))
data_glenelg_plot <- rename(data_glenelg_plot, fox_predicted = fit,  fox_predicted_se = se.fit) # rename

# fox plot
g_fox_plot <- ggplot(aes(x, y, fill = fox_predicted),
                     data = data_glenelg_plot) +
  geom_tile()+
  scale_fill_viridis("log(fox)", option = "viridis", limits = c(-2.19917294, 0.8425221)) +
  geom_point(data = records_glenelg, fill = NA, col = "white", size = 0.7, alpha = 0.7, shape = 3) +
  theme_bw(10) + 
  ggtitle("Glenelg region, 2018") +
  theme(axis.title = element_blank()) + 
  annotate("text", x = mean(data_glenelg_plot$x) + 6200, y = mean(data_glenelg_plot$y) + 13000, label = "Replicate 1") + 
  annotate("text", x = mean(data_glenelg_plot$x) - 10000, y = mean(data_glenelg_plot$y) - 500, label = "Replicate 1") +
  annotate("text", x = mean(data_glenelg_plot$x) - 7200, y = mean(data_glenelg_plot$y) + 26000, label = "Replicate 2") +
  annotate("text", x = mean(data_glenelg_plot$x) + 10000, y = mean(data_glenelg_plot$y) - 7500, label = "Replicate 2") + 
  annotate("text", x = mean(data_glenelg_plot$x) - 9000, y = mean(data_glenelg_plot$y) + 13000, label = "NI") + 
  annotate("text", x = mean(data_glenelg_plot$x) + 3450, y = mean(data_glenelg_plot$y) + 300, label = "NI") + 
  annotate("text", x = mean(data_glenelg_plot$x) - 12500, y = mean(data_glenelg_plot$y) - 14200, label = "I") + 
  annotate("text", x = mean(data_glenelg_plot$x) + 9000, y = mean(data_glenelg_plot$y) - 21500, label = "I") 
g_fox_plot



# OTWAY PLOTS -----------------------------------------------------------
## STEP 1) Make buffer zone around cameras to restrict plotting - optional 
#take just the cams
records_otways_cams <- distinct(records_otways, station, .keep_all = TRUE)
# change it to sf class
records_otways_cams <- st_as_sf(records_otways_cams, coords = c("x", "y"), crs = 32754) 
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
          length=50),
  y = seq(min(otways_buffer_df$y),
          max(otways_buffer_df$y),
          length=50),
  station = "T053", 
  station_year = "T053x2017",
  year = seq(2017,2019,by=1),
  survey_duration = 60)

# subset to just locations within buffer zone
data_otways_plot = data_otways_plot[with(data_otways_plot, inSide(otways_buffer_df, x, y)),]

data_otways_plot <- cbind(data_otways_plot, predict.gam(gam_o_fox, newdata = data_otways_plot, se.fit = TRUE, type = "link", exclude = c("s(station)", "s(survey_duration)")))
data_otways_plot <- rename(data_otways_plot, fox_predicted = fit,  fox_predicted_se = se.fit) # rename


## STEP 3) PLOTS
# fox plot
par(mar = c(5.1, 20, 4.1, 2.1))
o_fox_plot <- ggplot(aes(x, y, fill = fox_predicted),
                     data = data_otways_plot) +
  geom_tile()+
  scale_fill_viridis("log(fox)", option = "viridis", limits = c(-2.19917294, 0.8425221)) +
  facet_wrap(~year, nrow = 1) +
  geom_point(data = records_otways, fill = NA, col = "white", size = 0.7, alpha = 0.7, shape = 3) +
  theme_bw(10) + 
  ggtitle("Otway region") +
  theme(axis.title = element_blank()) + 
  annotate("text", x = mean(data_otways_plot$x) + 7500, y = mean(data_otways_plot$y) + 7500, label = "NI") + 
  annotate("text", x = mean(data_otways_plot$x) - 7500, y = mean(data_otways_plot$y) - 7500, label = "I") 
o_fox_plot



# COMBINE AND PLOT --------------------------------------------------------

png("C2-manuscript/figs/fig2A_600dpi.png", width = 7, height = 7, res = 600, units = "in")
g_fox_plot +
  plot_annotation(title = 'A') 
dev.off()

png("C2-manuscript/figs/fig2B_600dpi.png", width = 10, height = 4.75, res = 600, units = "in")
o_fox_plot +
  plot_annotation(title = 'B') 
dev.off()

# to merge, using imagemagick - type in the terminal:
# convert C2-manuscript/figs/fig2A_600dpi.png C2-manuscript/figs/fig2B_600dpi.png -gravity center -append C2-manuscript/figs/fig2_600dpi.png
