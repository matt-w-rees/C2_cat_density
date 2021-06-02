
rm(list = ls())
options(scipen = 999)

library(secr)
library(dplyr)
library(mgcv)
library(DHARMa)


# LOAD DATA ------------------------------------------------------------
# load pre-prepared fox record table
#records <- read_csv("/Users/mrees2/Dropbox/personal/matt/github/merge-records-add-variables-reformat/derived_data/records_presence_absence.csv")
## subset to my data
#records <- filter(records, data_source == "matt")
#records <- subset(records, select = c("region", "station_year", "station", "fox", "cat", "potoroo_ln", "bandicoot_sb", "date_start", "date_end", "survey_duration", "year", "latitude", "longitude", "XGROUPNAME", "foxbaits"))
#records$station_year <- gsub("_", "x", records$station_year)
#head(records)
## add xy coords
#records_sp <- st_as_sf(records, coords = c("longitude", "latitude"), crs = 4326) %>% 
#  st_transform(crs = 32754)
#records$x <- st_coordinates(records_sp)[,1]
#records$y <- st_coordinates(records_sp)[,2]
## save
#write.csv(records, "derived_data/spp_records_pa.csv")
#rm(records_sp)

# load records
records <- read.csv("derived_data/spp_records_pa.csv")

# capthists
mrch_a  <- readRDS("derived_data/capthists/mrch_a.RData")
mrch_c  <- readRDS("derived_data/capthists/mrch_c.RData")
mrch_h  <- readRDS("derived_data/capthists/mrch_h.RData")
mrch_mc <- readRDS( "derived_data/capthists/mrch_mc.RData")
mrch_s_17 <- readRDS("derived_data/capthists/mrch_s_17.RData")
mrch_n_17 <- readRDS("derived_data/capthists/mrch_n_17.RData")
mrch_s_18 <- readRDS("derived_data/capthists/mrch_s_18.RData")
mrch_n_18 <- readRDS("derived_data/capthists/mrch_n_18.RData")
mrch_s_19 <- readRDS("derived_data/capthists/mrch_s_19.RData")
mrch_n_19 <- readRDS("derived_data/capthists/mrch_n_19.RData")


# FIT GAMS ----------------------------------------------------------------
# specify records class first
records <- mutate(records, fox = as.integer(fox),
                           year = factor(year, ordered = FALSE),
                           region = factor(region, ordered = FALSE),
                           station = factor(station, ordered = FALSE),
                           survey_duration = as.integer(survey_duration))

# subset to data to each region
records_glenelg <- filter(records, region == "glenelg")
records_otways <- filter(records, region == "otways")

# fit GAM - glenelg 
gam_g_fox <- bam(fox ~ s(x, y, bs = "ds", m = c(1, 0.5), k = 200) + 
                       offset(log(survey_duration)), 
                 data = records_glenelg, family = binomial, discrete = TRUE) 
plot(gam_g_fox, pages = 1, scheme = 2)
summary(gam_g_fox)

# fit GAM - otways 
gam_o_fox <- bam(fox ~ s(x, y, year, bs = "fs", xt = list(bs = "ds", m = c(1, 0.5)), k = 100) + 
                       s(station, bs = "re") + 
                       offset(log(survey_duration)), 
                 data = records_otways, family = binomial)   
summary(gam_o_fox)



# GLENELG PREDICTIONS -----------------------------------------------------
# PREDICT GAM INTO HABITAT MASK -------------------------------------------
# separate masks per session
mask_a  <- make.mask(traps(mrch_a),  buffer = 4000, spacing = 474 * 0.6, type = 'trapbuffer')     
mask_c  <- make.mask(traps(mrch_c),  buffer = 4000, spacing = 474 * 0.6, type = 'trapbuffer')     
mask_h  <- make.mask(traps(mrch_h),  buffer = 4000, spacing = 474 * 0.6, type = 'trapbuffer')     
mask_mc <- make.mask(traps(mrch_mc),  buffer = 4000, spacing = 474 * 0.6, type = 'trapbuffer') 

# make multisession files
mrch_glenelg <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc)
glenelg_mask <- make.mask(traps(mrch_glenelg), buffer = 4000, spacing = 474 * 0.6, type = 'trapbuffer')       

# make a new dataframe with each habitat mask cell so we can add covariates from the GAMs
glenelg_mask_df <- do.call(rbind.data.frame, glenelg_mask)

# add survey duration dummy variable
glenelg_mask_df$survey_duration <- 60 

# predict GAM into the dataframe 
glenelg_mask_df <- cbind(glenelg_mask_df, predict.gam(gam_g_fox, newdata = glenelg_mask_df, se.fit = TRUE, type = "response"))

## convert df to list per each session
# but first give row values a column name
glenelg_mask_df <- setNames(cbind(rownames(glenelg_mask_df), glenelg_mask_df, row.names = NULL), c("sess", "x", "y", "station", "fox_predicted", "fox_predicted_se"))

# make sess column the same per session
glenelg_mask_df$sess <- as.character(glenelg_mask_df$sess)
glenelg_mask_df$sess <- gsub('[0-9]+', '', glenelg_mask_df$sess)
glenelg_mask_df$sess <- gsub('[.]+', '', glenelg_mask_df$sess)
#head(glenelg_mask_df)


# because we have a multi-session CH - we need to add in covariates separately for each session then merge back into multi-session mask...
# split dataframe by session 
glenelg_mask_df_a   <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_a"),]
glenelg_mask_df_c   <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_c"),]
glenelg_mask_df_h   <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_h"),]
glenelg_mask_df_mc  <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_mc"),]

# add covariates to each session mask separately 
covariates(mask_a) <- glenelg_mask_df_a[,5:6]
covariates(mask_c) <- glenelg_mask_df_c[,5:6]
covariates(mask_h) <- glenelg_mask_df_h[,5:6]
covariates(mask_mc) <- glenelg_mask_df_mc[,5:6]


# ADD FOX VALUES AS A TRAP COVARIATE --------------------------------------
# as we are going to test whether foxes influence cat detectability - we need detectability covariates to be input as trap covariates
# convert traps list into dataframe and do some cleaning
glenelg_traps_df <- do.call(rbind.data.frame, traps(mrch_glenelg))
glenelg_traps_df <- setNames(cbind(rownames(glenelg_traps_df), glenelg_traps_df, row.names = NULL), c("station", "x", "y"))
glenelg_traps_df$station <- as.character(glenelg_traps_df$station)
glenelg_traps_df$station <- gsub('mrch_mc.', "", glenelg_traps_df$station)
glenelg_traps_df$station <- gsub('mrch_h.', "", glenelg_traps_df$station)
glenelg_traps_df$station <- gsub('mrch_a.', "", glenelg_traps_df$station)
glenelg_traps_df$station <- gsub('mrch_c.', "", glenelg_traps_df$station)

# add survey duration dummy variable
glenelg_traps_df$survey_duration <- 60 

# predict GAM model values into the traps dataframe
# predict into the dataframe - exclude fox predictor in the cat prey model to avoid double counting this effect
glenelg_traps_df <- cbind(glenelg_traps_df, predict.gam(gam_g_fox, newdata = glenelg_traps_df, se.fit = TRUE, type = "response"))

# rename fox_predicted to fox_predicted_trapcov for model specification
colnames(glenelg_traps_df) <- c("station", "x", "y", "survey_duration", "fox_predicted_trapcov", "fox_predicted_trapcov_se")

# need to split dataframes and capthists into single sessions to add covariates 
glenelg_traps_df_a   <- glenelg_traps_df[which(grepl("^A", glenelg_traps_df$station)),]
glenelg_traps_df_c   <- glenelg_traps_df[which(grepl("^C", glenelg_traps_df$station)),]
glenelg_traps_df_h   <- glenelg_traps_df[which(grepl("^H", glenelg_traps_df$station)),]
glenelg_traps_df_mc  <- glenelg_traps_df[which(grepl("^M", glenelg_traps_df$station)),]

# add values into capthists - fox
covariates(traps(mrch_a)) <- glenelg_traps_df_a[,5:6]
covariates(traps(mrch_c)) <- glenelg_traps_df_c[,5:6]
covariates(traps(mrch_h)) <- glenelg_traps_df_h[,5:6]
covariates(traps(mrch_mc)) <- glenelg_traps_df_mc[,5:6]



# OTWAY PREDICTIONS -----------------------------------------------------
# merge to single capthist (in order of deployment) and make mask
mrch_otways <- MS.capthist(mrch_s_17, mrch_n_17, mrch_s_18, mrch_n_18, mrch_s_19, mrch_n_19)
otways_mask <- make.mask(traps(mrch_otways), buffer = 4000, spacing = 312 * 0.6, type = 'trapbuffer')       

# also separate masks per session
mask_n_17 <- make.mask(traps(mrch_n_17), buffer = 4000, spacing = 312 * 0.6, type = 'trapbuffer')     
mask_n_18 <- make.mask(traps(mrch_n_18), buffer = 4000, spacing = 312 * 0.6, type = 'trapbuffer')     
mask_n_19 <- make.mask(traps(mrch_n_19), buffer = 4000, spacing = 312 * 0.6, type = 'trapbuffer')   
mask_s_17 <- make.mask(traps(mrch_s_17), buffer = 4000, spacing = 312 * 0.6, type = 'trapbuffer')     
mask_s_18 <- make.mask(traps(mrch_s_18), buffer = 4000, spacing = 312 * 0.6, type = 'trapbuffer')     
mask_s_19 <- make.mask(traps(mrch_s_19), buffer = 4000, spacing = 312 * 0.6, type = 'trapbuffer')   


# make a new dataframe with each habitat mask cell so we can add covariates from the GAMs
otways_mask_df <- do.call(rbind.data.frame, otways_mask)
# make sess column the same per session
otways_mask_df <- tibble::rownames_to_column(otways_mask_df, "sess")
otways_mask_df$sess <- as.character(otways_mask_df$sess)
otways_mask_df$sess <- substr(otways_mask_df$sess, 1, 9)
unique(otways_mask_df$sess)

# add a year column 
otways_mask_df$year <- ifelse(otways_mask_df$sess == "mrch_n_17", "2017", NA)
otways_mask_df$year <- ifelse(otways_mask_df$sess == "mrch_s_17", "2017", otways_mask_df$year)
otways_mask_df$year <- ifelse(otways_mask_df$sess == "mrch_n_18", "2018", otways_mask_df$year)
otways_mask_df$year <- ifelse(otways_mask_df$sess == "mrch_s_18", "2018", otways_mask_df$year)
otways_mask_df$year <- ifelse(otways_mask_df$sess == "mrch_n_19", "2019", otways_mask_df$year)
otways_mask_df$year <- ifelse(otways_mask_df$sess == "mrch_s_19", "2019", otways_mask_df$year)
#head(otways_mask_df)

# add dummy variable for random effect and survey duration (note exclude in predict function below - but we still need to give it a column)
otways_mask_df$station <- "T052" 
otways_mask_df$survey_duration <- 60 

# predict into the dataframe 
otways_mask_df <- cbind(otways_mask_df, predict.gam(gam_o_fox, newdata = otways_mask_df, se.fit = TRUE, type = "response", exclude = c("s(station)")))
otways_mask_df <- rename(otways_mask_df, fox_predicted = fit,  fox_predicted_se = se.fit) # rename
#head(otways_mask_df)

# because we have a multi-session CH - we need to add in covariates seperately for each session then merge back into multi-session mask...
#split dataframe by session 
otways_mask_df_n_17 <- otways_mask_df[which(otways_mask_df$sess == "mrch_n_17"),]
otways_mask_df_s_17 <- otways_mask_df[which(otways_mask_df$sess == "mrch_s_17"),]
otways_mask_df_n_18 <- otways_mask_df[which(otways_mask_df$sess == "mrch_n_18"),]
otways_mask_df_s_18 <- otways_mask_df[which(otways_mask_df$sess == "mrch_s_18"),]
otways_mask_df_n_19 <- otways_mask_df[which(otways_mask_df$sess == "mrch_n_19"),]
otways_mask_df_s_19 <- otways_mask_df[which(otways_mask_df$sess == "mrch_s_19"),]

# add covariates to each session mask separately 
covariates(mask_n_17) <- otways_mask_df_n_17[,7:8]
covariates(mask_n_18) <- otways_mask_df_n_18[,7:8]
covariates(mask_n_19) <- otways_mask_df_n_19[,7:8]
covariates(mask_s_17) <- otways_mask_df_s_17[,7:8]
covariates(mask_s_18) <- otways_mask_df_s_18[,7:8]
covariates(mask_s_19) <- otways_mask_df_s_19[,7:8]



# ADD FOX VALUES AS A TRAP COVARIATE --------------------------------------
# as we are going to test whether foxes influence cat detectability - we need detectability covariates for input as trap covariates
# convert traps list into dataframe and do some cleaning
otways_traps_df <- do.call(rbind.data.frame, traps(mrch_otways))
otways_traps_df <- setNames(cbind(rownames(otways_traps_df), otways_traps_df, row.names = NULL), c("station_year", "x", "y"))
otways_traps_df$station2 <- as.character(otways_traps_df$station_year) # need this column for later!
otways_traps_df$station_year <- as.character(otways_traps_df$station_year)
otways_traps_df$year <- substr(otways_traps_df$station_year, nchar(otways_traps_df$station_year) - 3, nchar(otways_traps_df$station_year))
# more cleaning 
otways_traps_df$station_year <- gsub('mrch_n_17.', "", otways_traps_df$station_year)
otways_traps_df$station_year <- gsub('mrch_n_18.', "", otways_traps_df$station_year)
otways_traps_df$station_year <- gsub('mrch_n_19.', "", otways_traps_df$station_year)
otways_traps_df$station_year <- gsub('mrch_s_17.', "", otways_traps_df$station_year)
otways_traps_df$station_year <- gsub('mrch_s_18.', "", otways_traps_df$station_year)
otways_traps_df$station_year <- gsub('mrch_s_19.', "", otways_traps_df$station_year)
#head(otways_traps_df)
# add dummy variables for GAM prediction
otways_traps_df$station <- "T052"
otways_traps_df$survey_duration <- 60

# predict GAM model values into the traps dataframe
otways_traps_df <- cbind(otways_traps_df, predict.gam(gam_o_fox, newdata = otways_traps_df, se.fit = TRUE, type = "response", exclude = c("s(station)")))
otways_traps_df <- rename(otways_traps_df, fox_predicted_trapcov = fit,  fox_predicted_trapcov_se = se.fit) # rename
#head(otways_traps_df)

# need to split dataframes and capthists into single sessions to add covariates 
otways_traps_df_n_17 <- otways_traps_df[which(grepl("^mrch_n_17", otways_traps_df$station2)),]
otways_traps_df_n_18 <- otways_traps_df[which(grepl("^mrch_n_18", otways_traps_df$station2)),]
otways_traps_df_n_19 <- otways_traps_df[which(grepl("^mrch_n_19", otways_traps_df$station2)),]
otways_traps_df_s_17 <- otways_traps_df[which(grepl("^mrch_s_17", otways_traps_df$station2)),]
otways_traps_df_s_18 <- otways_traps_df[which(grepl("^mrch_s_18", otways_traps_df$station2)),]
otways_traps_df_s_19 <- otways_traps_df[which(grepl("^mrch_s_19", otways_traps_df$station2)),]

# add values into capthists 
covariates(traps(mrch_n_17)) <- otways_traps_df_n_17[,8:9]
covariates(traps(mrch_n_18)) <- otways_traps_df_n_18[,8:9]
covariates(traps(mrch_n_19)) <- otways_traps_df_n_19[,8:9]
covariates(traps(mrch_s_17)) <- otways_traps_df_s_17[,8:9]
covariates(traps(mrch_s_18)) <- otways_traps_df_s_18[,8:9]
covariates(traps(mrch_s_19)) <- otways_traps_df_s_19[,8:9]
#head(covariates(traps(mrch_s_17)))







# SAVE EVERYTHING -------------------------------------------------------------
# save GAMs
saveRDS(gam_g_fox, "models/gam_g_fox.RData")  
saveRDS(gam_o_fox, "models/gam_o_fox.RData")  

# save mask dataframes  (useful for plotting etc in this format)
saveRDS(glenelg_mask_df, "derived_data/glenelg_mask_df.RData")
saveRDS(otways_mask_df, "derived_data/otway_mask_df.RData")

# convert into multi-session habitat mask (same order as capthist!)
mask_glenelg <- list(mask_a, mask_c, mask_h, mask_mc)
mask_otways <- list(mask_s_17, mask_n_17, mask_s_18, mask_n_18, mask_s_19, mask_n_19)
mrch_glenelg <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc)
mrch_otways <- MS.capthist(mrch_s_17, mrch_n_17, mrch_s_18, mrch_n_18, mrch_s_19, mrch_n_19)
# check covs all worked
#head(covariates(mask_otways)) 
#head(covariates(traps(mrch_glenelg[1])))
#summary(mrch_glenelg, terse = TRUE)
# now save 
saveRDS(mask_glenelg, "derived_data/mask_glenelg.RData")
saveRDS(mask_otways, "derived_data/mask_otways.RData")
saveRDS(mrch_glenelg, "derived_data/mrch_glenelg.RData")
saveRDS(mrch_otways, "derived_data/mrch_otways.RData")

## END