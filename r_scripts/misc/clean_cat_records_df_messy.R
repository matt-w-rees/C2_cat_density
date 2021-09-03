
library(dplyr)

glenelg <- bind_rows(rt_id_a, rt_id_c, rt_id_h, rt_id_mc,     
               rt_unid_a, rt_unid_c, rt_unid_h, rt_unid_mc, 
               rt_unm_a, rt_unm_c, rt_unm_h, .id = NULL)

glenelg$region <- "glenelg"

otways <- bind_rows(rt_id_n_17, rt_id_n_18, rt_id_n_19, rt_id_s_17, rt_id_s_18, rt_id_s_19, 
                     rt_unid_n_17, rt_unid_n_18, rt_unid_n_19, rt_unid_s_17, rt_unid_s_18, rt_unid_s_19, 
                     rt_unm_n_17,  rt_unm_n_18,  rt_unm_n_19,  rt_unm_s_17,  rt_unm_s_18,  rt_unm_s_19, .id = NULL)

otways$region <- "otways"

cat_records <- rbind(glenelg, otways)
head(cat_records)


# clean
cat_records2 <- subset(cat_records, select = -c(Date, Time, delta.time.secs, delta.time.mins, delta.time.hours, delta.time.days, Directory, FileName))
cat_records2 <- subset(cat_records2, select = -c(HierarchicalSubject))



records <- cat_records2
# load camdata (site information )
camdata <- read.csv("/Users/mrees2/Dropbox/personal/matt/github/compile-camera-records-matt/raw_data/camdata_matt.csv")
# merge
records <- left_join(records, camdata)
# clean
records$date_time <- records$DateTimeOriginal
records$DateTimeOriginal <- NULL
records$station <- records$Station
records$Station <- NULL
# remove year from station col in the otways
records$station <- as.character(records$station)
records$station <- ifelse(records$region == "otways", substr(records$station, 1, nchar(records$station) - 5), records$station)
head(unique(records$station))
# add station_year col
library(lubridate)
records$station_year <- paste0(records$station, "_", year(ymd_hms(records$date_time)))
head(unique(records$station_year))

# FIX DATE MISTAKES -------------------------------------------------------
# move H80_2018 retrieval date forward 2 days
records$date_end <- if_else(records$station_year == "H80_2018", ymd(records$date_end) + days(2), ymd(records$date_end))
# move H80_2018 retrieval date forward 6 days
records$date_end <- if_else(records$station_year == "T085_2019", ymd(records$date_end) + days(6), ymd(records$date_end))
# fix U024_2018 small mammal records (long story)


# REORDER & SAVE -------------------------------------------------------
# reorder
records <- records %>%
  arrange(region, date_start, station_year, date_time)


head(records)
records$fox_control <- NULL

records <- records %>% rename(id_status = Species, individual = metadata_individual_ID, coat_type = metadata_cat_coat, latitude = lat, longitude = long)

records$date_start <- ymd(records$date_start)
head(records)
unique(records$individual)


write.csv(records, "derived_data/cat_id_detections.csv")
write.csv(records, "/Users/mrees2/Dropbox/personal/matt/github/diel-activity-invasive-predators/raw_data/cat_id_detections.csv")


