# REFORMAT SPECIES RECORDS DATA FOR GAMS

rm(list = ls())
options(scipen = 999)

library(readr)
library(dplyr)
library(reshape2)
library(sf)
library(mgcv)
library(DHARMa)

#load species records table
records <- read_csv("/Users/mrees2/Dropbox/personal/matt/github/compile-camera-records-matt/derived_data/records_matt_clean.csv")
head(records)


# GROUP PREY ------------------------------------------------------
# couldn't be identified as dasyurid or rodent in many cases, so lump them all together (along with eastern pygmy possums, as they are more rodent than brushtail like in the eye of a cat)
records$species <- if_else(records$species == "antechinus_agile" | records$species == "antechinus_dusky" | records$species == "antechinus_swamp" | records$species == "antechinus_unid" | records$species == "dasyurid_sp" | records$species == "dunnart_wf" | records$species == "mouse_house" | records$species == "rat_black" | records$species == "rat_bush" | records$species == "rat_swamp" | records$species == "rodent_unid" | records$species == "small_mammal" | records$species == "possum_ep",   
                           "small_mammal", records$species)

# group long-nosed and southern brown bandicoots in the Otways -- very few SBB detections here (Only SBBs in glenelg)
records$species <- if_else(records$species == "bandicoot_ln" | records$species == "bandicoot_sb", "bandicoot", records$species)

# add xy coordinates
records_sf <- st_as_sf(records, coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = 32754)
records$x <- st_coordinates(records_sf)[,1]
records$y <- st_coordinates(records_sf)[,2]
head(records)


## REFORMAT RECORDS FOR GAMS ----------------------------------------------------------------------------
records$occ <- 1
# melt the dataframe
gam_data <- melt(records, id.var = colnames(records), measure.var = "occ")
# cast to derive presence-absence per year per station_year, long format
gam_data = dcast(gam_data, station_year ~ species, fill = 0, fun = min)
# add in site covariates
sites_covs <- subset(records, select = c("station", "x", "y", "grid", "region", "year", "station_year", "survey_duration"))
# merge gam_data with site covs
gam_data <- merge(gam_data, sites_covs, by = "station_year", all.y = TRUE)
# remove duplicate rows
gam_data <- distinct(gam_data, .keep_all = TRUE)
# no NA's?
unique(is.na(gam_data))

# also remove species not of interest
gam_data <- subset(gam_data, select = -c(bird, deer_unid, echidna, emu, false_trigger, fault, glider_ft, glider_sugar, goat, kangaroo_eg, kangaroo_unid, koala, macropod_unknown, other, people, pig, unknown, wallaby_rn, wallaby_swamp))
head(gam_data)

# save
write.csv(gam_data, "derived_data/gam_data.csv")

# specify factors for models
gam_data$year <- factor(gam_data$year, ordered = FALSE)
gam_data$station <- factor(gam_data$station, ordered = FALSE)

# split by region
gam_data_glenelg <- filter(gam_data, region == "glenelg")
gam_data_otways <- filter(gam_data, region == "otways")



# GAMS - GLENELG -------------------------------------------------------
# fox
gam_g_fox <- bam(fox ~ s(x, y, bs = "ds", m = c(1, 0.5), k = 100) + 
                   offset(log(survey_duration)), 
                 data = gam_data_glenelg, family = binomial, discrete = TRUE)   

# small mammal
gam_g_sm <- bam(small_mammal ~ s(x, y, bs = "ds", m = c(1, 0.5), k = 100) + 
                  offset(log(survey_duration)), 
                data = gam_data_glenelg, family = binomial, discrete = TRUE)   

# southern brown bandicoot
gam_g_sbb <- bam(bandicoot ~ s(x, y, bs = "ds", m = c(1, 0.5), k = 100) + 
                   offset(log(survey_duration)), 
                 data = gam_data_glenelg, family = binomial, discrete = TRUE)   

# long-nosed potoroo
gam_g_pot <- bam(potoroo_ln ~ s(x, y, bs = "ds", m = c(1, 0.5), k = 100) + 
                   offset(log(survey_duration)), 
                 data = gam_data_glenelg, family = binomial, discrete = TRUE)   

# brushtail possum
gam_g_btp <- bam(possum_brushtail ~ s(x, y, bs = "ds", m = c(1, 0.5), k = 100) + 
                   offset(log(survey_duration)), 
                 data = gam_data_glenelg, family = binomial, discrete = TRUE)   

# ringtail possum
gam_g_rtp <- bam(possum_ringtail ~ s(x, y, bs = "ds", m = c(1, 0.5), k = 100) + 
                   offset(log(survey_duration)), 
                 data = gam_data_glenelg, family = binomial, discrete = TRUE)   



# GAMS - OTWAYS -------------------------------------------------------

# fox
gam_o_fox <- bam(fox ~ year + s(x, y, by = year,  bs = "ds", m = c(1, 0.5), k = 100) + 
                    s(station, bs = "re") + 
                    offset(log(survey_duration)), 
                  data = gam_data_otways, family = binomial, discrete = TRUE)   

# small mammal
gam_o_sm <- bam(small_mammal ~ year + s(x, y, by = year,  bs = "ds", m = c(1, 0.5), k = 100) + 
                   s(station, bs = "re") + 
                   offset(log(survey_duration)), 
                 data = gam_data_otways, family = binomial, discrete = TRUE)   

# long-nosed and (a couple) southern brown bandicoots
gam_o_band <- bam(bandicoot ~ year + s(x, y, by = year,  bs = "ds", m = c(1, 0.5), k = 100) + 
                     s(station, bs = "re") + 
                     offset(log(survey_duration)), 
                   data = gam_data_otways, family = binomial, discrete = TRUE)   

# long-nosed potoroo
gam_o_pot <- bam(potoroo_ln ~ year + s(x, y, by = year,  bs = "ds", m = c(1, 0.5), k = 100) + 
                    s(station, bs = "re") + 
                    offset(log(survey_duration)), 
                  data = gam_data_otways, family = binomial, discrete = TRUE)   

# brushtail possum
gam_o_btp <- bam(possum_brushtail ~ year + s(x, y, by = year,  bs = "ds", m = c(1, 0.5), k = 100) + 
                    s(station, bs = "re") + 
                    offset(log(survey_duration)), 
                  data = gam_data_otways, family = binomial, discrete = TRUE)   

# ringtail possum
gam_o_rtp <- bam(possum_ringtail ~ year + s(x, y, by = year,  bs = "ds", m = c(1, 0.5), k = 100) + 
                    s(station, bs = "re") + 
                    offset(log(survey_duration)), 
                  data = gam_data_otways, family = binomial, discrete = TRUE)   



# SAVE MODELS -------------------------------------------------------------
# glenelg
saveRDS(gam_g_fox, "models/gams/glenelg/gam_g_fox.RData")  
saveRDS(gam_g_sm , "models/gams/glenelg/gam_g_sm.RData") 
saveRDS(gam_g_sbb, "models/gams/glenelg/gam_g_sbb.RData")  
saveRDS(gam_g_pot, "models/gams/glenelg/gam_g_pot.RData")  
saveRDS(gam_g_btp, "models/gams/glenelg/gam_g_btp.RData")  
saveRDS(gam_g_rtp, "models/gams/glenelg/gam_g_rtp.RData")  

#otways
saveRDS(gam_o_fox, "models/gams/otways/gam_o_fox.RData")  
saveRDS(gam_o_sm , "models/gams/otways/gam_o_sm.RData") 
saveRDS(gam_o_band, "models/gams/otways/gam_o_band.RData")  
saveRDS(gam_o_pot, "models/gams/otways/gam_o_pot.RData")  
saveRDS(gam_o_btp, "models/gams/otways/gam_o_btp.RData")  
saveRDS(gam_o_rtp, "models/gams/otways/gam_o_rtp.RData")  

# END