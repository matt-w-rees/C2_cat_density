# EXTRACT SECR DATA FROM FERAL CATS IMAGES (ID'S TAGGED VIA METADATA) - GLENELG & OTWAYS (1 SESSION PER GRID)

library(readr)
library(camtrapR)
library(dplyr)
library(sf)
library(secr)

# Load survey info and do some cleaning (dirtying more like it)
#camdata <- read_csv("/Users/mrees2/Dropbox/personal/matt/github/compile-camera-records-matt/derived_data/capthists/records_matt_clean.csv")
camdata <- read_csv("/Users/ree140/Dropbox/personal/matt/github/unlinked/compile-camera-records-matt/derived_data/records_matt_clean.csv") 
camdata <- distinct(camdata, station_year, .keep_all = TRUE)
camdata$Station <- ifelse(camdata$region == "otways", gsub("_", "x", camdata$station_year), camdata$station)
camdata$block <- camdata$grid 

#  Load CT model info and do some cleaning 
cam_models <- read_csv("/Users/ree140/Dropbox/personal/matt/github/unlinked/compile-camera-records-matt/derived_data/sp_records_model.csv") %>%
  select(c(Station, model = Model))  %>%
  distinct(Station, .keep_all = TRUE)

# fix na's where cam model is diplayed in picture
cam_models$model <- if_else(cam_models$Station == "T005x2019", "PC900 PROFESSIONAL", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T034x2019", "PC900 PROFESSIONAL", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T035x2019", "HC600 HYPERFIRE", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T041x2019", "HC600 HYPERFIRE", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T045x2019", "PC900 PROFESSIONAL", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T050x2019", "HC600 HYPERFIRE", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T062x2019", "PC900 PROFESSIONAL", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T099x2019", "PC900 PROFESSIONAL", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "T100x2019", "HC500 HYPERFIRE", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "TN01x2019", "PC900 PROFESSIONAL", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "TN10x2019", "HC600 HYPERFIRE", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "TN12x2019", "HC600 HYPERFIRE", cam_models$model)
cam_models$model <- if_else(cam_models$Station == "UND11x2019", "HC600 HYPERFIRE", cam_models$model)
# hows that leave us?
table(is.na(cam_models$model))
# 30 stations with unknown NA's - I know none of them are HYPERFIRE HYPERFIRE 2, can't leave them as NA - so assume HC600 HYPERFIRE
cam_models$model <- if_else(is.na(cam_models$model), "HC600 HYPERFIRE", cam_models$model)
## append camdata with models
camdata <- left_join(camdata, cam_models)
# lets make cam models per region into a nice dataframe for an appendix table
df <- as.data.frame(xtabs(~camdata$model + camdata$region)) %>% 
  rename("Camera Model" = camdata.model, "Region" = camdata.region, "Frequency" = Freq) %>%
  relocate("Region") %>%
  filter(!(Frequency == 0))
df$Region <- if_else(df$Region == "glenelg", "Glenelg", "Otway")
write.csv(df, "derived_data/cam_models.csv")

# also add xy coords
camdata_sf <- st_as_sf(camdata, coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = 32754)
camdata$xcoord <- st_coordinates(camdata_sf)[,1]
camdata$ycoord <- st_coordinates(camdata_sf)[,2]


# LOAD IMAGES ---------------------------------------------------------------
images <- file.path("raw_data/tagged_cat_images")
length(list.files(images, pattern = "JPG", recursive = TRUE))

images_mc   <- file.path("raw_data/tagged_cat_images/glenelg/mt_clay")
images_h    <- file.path("raw_data/tagged_cat_images/glenelg/hotspur")
images_a    <- file.path("raw_data/tagged_cat_images/glenelg/annya")
images_c    <- file.path("raw_data/tagged_cat_images/glenelg/cobbob")
images_n_17 <- file.path("raw_data/tagged_cat_images/otways/2017/north")
images_n_18 <- file.path("raw_data/tagged_cat_images/otways/2018/north")
images_n_19 <- file.path("raw_data/tagged_cat_images/otways/2019/north")
images_s_17 <- file.path("raw_data/tagged_cat_images/otways/2017/south")
images_s_18 <- file.path("raw_data/tagged_cat_images/otways/2018/south")
images_s_19 <- file.path("raw_data/tagged_cat_images/otways/2019/south")


# PREPARE CAMERA-TRAP SURVEY INFO ----------------------------------------
# split for each grid
camdata_mc   <- camdata[which(camdata$block == "mtclay"),]
camdata_h    <- camdata[which(camdata$block == "hotspur"),]
camdata_a    <- camdata[which(camdata$block == "annya"),]
camdata_c    <- camdata[which(camdata$block == "cobb"),]
camdata_n_17 <- camdata[which(camdata$block == "north" & camdata$year == "2017"),]
camdata_n_18 <- camdata[which(camdata$block == "north" & camdata$year == "2018"),]
camdata_n_19 <- camdata[which(camdata$block == "north" & camdata$year == "2019"),]
camdata_s_17 <- camdata[which(camdata$block == "south" & camdata$year == "2017"),]
camdata_s_18 <- camdata[which(camdata$block == "south" & camdata$year == "2018"),]
camdata_s_19 <- camdata[which(camdata$block == "south" & camdata$year == "2019"),]

camop_mc    <- cameraOperation(CTtable = camdata_mc,   stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_h     <- cameraOperation(CTtable = camdata_h,    stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_a     <- cameraOperation(CTtable = camdata_a,    stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_c     <- cameraOperation(CTtable = camdata_c,    stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_n_17  <- cameraOperation(CTtable = camdata_n_17, stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_n_18  <- cameraOperation(CTtable = camdata_n_18, stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_n_19  <- cameraOperation(CTtable = camdata_n_19, stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_s_17  <- cameraOperation(CTtable = camdata_s_17, stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_s_18  <- cameraOperation(CTtable = camdata_s_18, stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camop_s_19  <- cameraOperation(CTtable = camdata_s_19, stationCol = "Station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")


# MAKE RECORD TABLES FOR ID'D CATS FOR EACH SESSION -----------------------
rt_id_mc   <- recordTable(inDir = images_mc,   IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_h    <- recordTable(inDir = images_h,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_a    <- recordTable(inDir = images_a,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_c    <- recordTable(inDir = images_c,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_n_17 <- recordTable(inDir = images_n_17, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_n_18 <- recordTable(inDir = images_n_18, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_n_19 <- recordTable(inDir = images_n_19, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_s_17 <- recordTable(inDir = images_s_17, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_s_18 <- recordTable(inDir = images_s_18, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_id_s_19 <- recordTable(inDir = images_s_19, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "unidentifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)


# MAKE RECORD TABLES FOR UNID'D CATS FOR EACH SESSION ---------------------
rt_unid_mc   <- recordTable(inDir = images_mc,   IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_h    <- recordTable(inDir = images_h,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_a    <- recordTable(inDir = images_a,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_c    <- recordTable(inDir = images_c,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_n_17 <- recordTable(inDir = images_n_17, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_n_18 <- recordTable(inDir = images_n_18, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_n_19 <- recordTable(inDir = images_n_19, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_s_17 <- recordTable(inDir = images_s_17, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_s_18 <- recordTable(inDir = images_s_18, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unid_s_19 <- recordTable(inDir = images_s_19, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unmarked", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)


# MAKE RECORD TABLES FOR UNMARKED CATS FOR EACH SESSION ---------------------
# no mt clay - no black cats detected here 
rt_unm_h    <- recordTable(inDir = images_h,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_a    <- recordTable(inDir = images_a,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_c    <- recordTable(inDir = images_c,    IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_n_17 <- recordTable(inDir = images_n_17, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_n_18 <- recordTable(inDir = images_n_18, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_n_19 <- recordTable(inDir = images_n_19, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_s_17 <- recordTable(inDir = images_s_17, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_s_18 <- recordTable(inDir = images_s_18, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)
rt_unm_s_19 <- recordTable(inDir = images_s_19, IDfrom = "metadata", metadataSpeciesTag = "mark_status", exclude = c("unidentifiable", "identifiable", "kitten"), timeZone = "Australia/Brisbane", writecsv = FALSE)


# MAKE SECR CAPTHIST ------------------------------------------------------
ch_mc    <- spatialDetectionHistory(recordTableIndividual = rt_id_mc,   camOp = camop_mc,   CTtable = camdata_mc,   species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_h     <- spatialDetectionHistory(recordTableIndividual = rt_id_h,    camOp = camop_h,    CTtable = camdata_h,    species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_a     <- spatialDetectionHistory(recordTableIndividual = rt_id_a,    camOp = camop_a,    CTtable = camdata_a,    species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_c     <- spatialDetectionHistory(recordTableIndividual = rt_id_c,    camOp = camop_c,    CTtable = camdata_c,    species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_n_17  <- spatialDetectionHistory(recordTableIndividual = rt_id_n_17, camOp = camop_n_17, CTtable = camdata_n_17, species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_n_18  <- spatialDetectionHistory(recordTableIndividual = rt_id_n_18, camOp = camop_n_18, CTtable = camdata_n_18, species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_n_19  <- spatialDetectionHistory(recordTableIndividual = rt_id_n_19, camOp = camop_n_19, CTtable = camdata_n_19, species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_s_17  <- spatialDetectionHistory(recordTableIndividual = rt_id_s_17, camOp = camop_s_17, CTtable = camdata_s_17, species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_s_18  <- spatialDetectionHistory(recordTableIndividual = rt_id_s_18, camOp = camop_s_18, CTtable = camdata_s_18, species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")
ch_s_19  <- spatialDetectionHistory(recordTableIndividual = rt_id_s_19, camOp = camop_s_19, CTtable = camdata_s_19, species = "cat", output = "binary", stationCol = "Station", stationCovariateCols = "model", speciesCol = "metadata_species", Xcol = "xcoord", Ycol = "ycoord", individualCol = "metadata_individual_ID", recordDateTimeCol = "DateTimeOriginal", recordDateTimeFormat = "%Y-%m-%d %H:%M:%S", occasionLength = 1,  occasionStartTime = 12, day1 = "station", includeEffort = TRUE, timeZone = "Australia/Brisbane")

# add markocc attribute for mark-resight
markocc(traps(ch_mc)) <- rep(0, 58)
markocc(traps(ch_h )) <- rep(0, 67)
markocc(traps(ch_a )) <- rep(0, 78)
markocc(traps(ch_c )) <- rep(0, 75)
markocc(traps(ch_n_17)) <- rep(0, 59)
markocc(traps(ch_n_18)) <- rep(0, 74)
markocc(traps(ch_n_19)) <- rep(0, 90)
markocc(traps(ch_s_17)) <- rep(0, 68)
markocc(traps(ch_s_18)) <- rep(0, 75)
markocc(traps(ch_s_19)) <- rep(0, 94)

# looks right?
par(mfrow = c(2,2))
plot(ch_a,  tracks = TRUE)
plot(ch_c,  tracks = TRUE)
plot(ch_h,  tracks = TRUE)
plot(ch_mc, tracks = TRUE)
par(mfrow = c(3,2))
plot(ch_n_17, tracks = TRUE)
plot(ch_s_17, tracks = TRUE)
plot(ch_n_18, tracks = TRUE)
plot(ch_s_18, tracks = TRUE)
plot(ch_n_19, tracks = TRUE)
plot(ch_s_19, tracks = TRUE)

# save capthist as RData file
saveRDS(ch_mc,   "derived_data/capthists/ch_mc.RData")
saveRDS(ch_h,    "derived_data/capthists/ch_h.RData")
saveRDS(ch_a,    "derived_data/capthists/ch_a.RData")
saveRDS(ch_c,    "derived_data/capthists/ch_c.RData")
saveRDS(ch_n_17, "derived_data/capthists/ch_n_17.RData")
saveRDS(ch_n_18, "derived_data/capthists/ch_n_18.RData")
saveRDS(ch_n_19, "derived_data/capthists/ch_n_19.RData")
saveRDS(ch_s_17, "derived_data/capthists/ch_s_17.RData")
saveRDS(ch_s_18, "derived_data/capthists/ch_s_18.RData")
saveRDS(ch_s_19, "derived_data/capthists/ch_s_19.RData")


# MAKE UNID'D MATRICES ----------------------------------------------------
unid_mc    <- detectionHistory(recordTable = rt_unid_mc,    camOp = camop_mc,    stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_h     <- detectionHistory(recordTable = rt_unid_h,     camOp = camop_h,     stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_a     <- detectionHistory(recordTable = rt_unid_a,     camOp = camop_a,     stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_c     <- detectionHistory(recordTable = rt_unid_c,     camOp = camop_c,     stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_n_17  <- detectionHistory(recordTable = rt_unid_n_17,  camOp = camop_n_17,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_n_18  <- detectionHistory(recordTable = rt_unid_n_18,  camOp = camop_n_18,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_n_19  <- detectionHistory(recordTable = rt_unid_n_19,  camOp = camop_n_19,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_s_17  <- detectionHistory(recordTable = rt_unid_s_17,  camOp = camop_s_17,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_s_18  <- detectionHistory(recordTable = rt_unid_s_18,  camOp = camop_s_18,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unid_s_19  <- detectionHistory(recordTable = rt_unid_s_19,  camOp = camop_s_19,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)

# change to dataframe
unid_mc   <- data.frame(unid_mc)        
unid_h    <- data.frame(unid_h)        
unid_a    <- data.frame(unid_a)        
unid_c    <- data.frame(unid_c)        
unid_n_17 <- data.frame(unid_n_17)        
unid_n_18 <- data.frame(unid_n_18)        
unid_n_19 <- data.frame(unid_n_19)        
unid_s_17 <- data.frame(unid_s_17)        
unid_s_18 <- data.frame(unid_s_18)        
unid_s_19 <- data.frame(unid_s_19)        

# replace NA's with 0
unid_mc[(is.na(unid_mc))] = 0
unid_h[(is.na(unid_h))] = 0
unid_a[(is.na(unid_a))] = 0
unid_c[(is.na(unid_c))] = 0
unid_n_17[(is.na(unid_n_17))] = 0
unid_n_18[(is.na(unid_n_18))] = 0
unid_n_19[(is.na(unid_n_19))] = 0
unid_s_17[(is.na(unid_s_17))] = 0
unid_s_18[(is.na(unid_s_18))] = 0
unid_s_19[(is.na(unid_s_19))] = 0

# save as csv
write.csv(unid_mc, file = "derived_data/capthists/unid_mc.csv")
write.csv(unid_h, file =  "derived_data/capthists/unid_h.csv")
write.csv(unid_a, file = "derived_data/capthists/unid_a.csv")
write.csv(unid_c, file = "derived_data/capthists/unid_c.csv")
write.csv(unid_n_17, file = "derived_data/capthists/unid_n_17.csv")
write.csv(unid_n_18, file = "derived_data/capthists/unid_n_18.csv")
write.csv(unid_n_19, file = "derived_data/capthists/unid_n_19.csv")
write.csv(unid_s_17, file = "derived_data/capthists/unid_s_17.csv")
write.csv(unid_s_18, file = "derived_data/capthists/unid_s_18.csv")
write.csv(unid_s_19, file = "derived_data/capthists/unid_s_19.csv")


# MAKE UNM'D MATRICES ----------------------------------------------------
# again - skip mt clay - no black cats
unm_h     <- detectionHistory(recordTable = rt_unm_h,     camOp = camop_h,     stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_a     <- detectionHistory(recordTable = rt_unm_a,     camOp = camop_a,     stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_c     <- detectionHistory(recordTable = rt_unm_c,     camOp = camop_c,     stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_n_17  <- detectionHistory(recordTable = rt_unm_n_17,  camOp = camop_n_17,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_n_18  <- detectionHistory(recordTable = rt_unm_n_18,  camOp = camop_n_18,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_n_19  <- detectionHistory(recordTable = rt_unm_n_19,  camOp = camop_n_19,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_s_17  <- detectionHistory(recordTable = rt_unm_s_17,  camOp = camop_s_17,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_s_18  <- detectionHistory(recordTable = rt_unm_s_18,  camOp = camop_s_18,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)
unm_s_19  <- detectionHistory(recordTable = rt_unm_s_19,  camOp = camop_s_19,  stationCol = "Station", speciesCol = "metadata_species", recordDateTimeCol = "DateTimeOriginal", species = "cat", occasionLength = 1, occasionStartTime = 12, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE)

# change to dataframe
unm_h    <- data.frame(unm_h)        
unm_a    <- data.frame(unm_a)        
unm_c    <- data.frame(unm_c)        
unm_n_17 <- data.frame(unm_n_17)        
unm_n_18 <- data.frame(unm_n_18)        
unm_n_19 <- data.frame(unm_n_19)        
unm_s_17 <- data.frame(unm_s_17)        
unm_s_18 <- data.frame(unm_s_18)        
unm_s_19 <- data.frame(unm_s_19)        

# replace NA's with 0
unm_h[(is.na(unm_h))] = 0
unm_a[(is.na(unm_a))] = 0
unm_c[(is.na(unm_c))] = 0
unm_n_17[(is.na(unm_n_17))] = 0
unm_n_18[(is.na(unm_n_18))] = 0
unm_n_19[(is.na(unm_n_19))] = 0
unm_s_17[(is.na(unm_s_17))] = 0
unm_s_18[(is.na(unm_s_18))] = 0
unm_s_19[(is.na(unm_s_19))] = 0

# save as csv
write.csv(unm_h,    file =  "derived_data/capthists/unm_h.csv")
write.csv(unm_a,    file = "derived_data/capthists/unm_a.csv")
write.csv(unm_c,    file = "derived_data/capthists/unm_c.csv")
write.csv(unm_n_17, file = "derived_data/capthists/unm_n_17.csv")
write.csv(unm_n_18, file = "derived_data/capthists/unm_n_18.csv")
write.csv(unm_n_19, file = "derived_data/capthists/unm_n_19.csv")
write.csv(unm_s_17, file = "derived_data/capthists/unm_s_17.csv")
write.csv(unm_s_18, file = "derived_data/capthists/unm_s_18.csv")
write.csv(unm_s_19, file = "derived_data/capthists/unm_s_19.csv")


# END 