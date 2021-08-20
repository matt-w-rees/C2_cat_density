## ADD SPATIAL COVARIATES TO SECR MASK AND TRAPFILES
# EVC GROUP (MASK ONLY)
# FOX SPATIAL OCCUPANCY (MASK AND TRAPS) VIA GAMS USING CAMERA-TRAP DETECTIONS 
# Matthew Rees


rm(list = ls())
options(scipen = 999)

library(secr)
library(dplyr)
library(mgcv)
library(DHARMa)
library(sp)
library(sf)
library(rgdal)


# ADD VEGETATION CLASS ----------------------------------------------------
## LOAD AND MANIPULATE SPATIAL VEGETATION TYPE LAYER
evc <- st_read("raw_data/vegetation_layer/EVC_32754_4000M.shp")
# subset to just EVC group name
evc <- subset(evc, select = c("XGROUPNAME"))
evc$XGROUPNAME <- as.character(evc$XGROUPNAME)
# merge cleared types - different name for the same thing
evc$XGROUPNAME <- if_else(evc$XGROUPNAME == "No native vegetation recorded", "CLEARED", evc$XGROUPNAME)
evc$XGROUPNAME <- gsub("CLEARED", "Cleared", evc$XGROUPNAME)

# merge rainforest / wet forest
evc$XGROUPNAME <- if_else(evc$XGROUPNAME == "Rainforests", "Wet or Damp Forests", evc$XGROUPNAME)
# merge heathlands / heathy woodlands
evc$XGROUPNAME <- if_else(evc$XGROUPNAME == "Heathlands", "Heathy Woodlands", evc$XGROUPNAME)
# merge dry forests / herb-rich woodlands and lowland forests
evc$XGROUPNAME <- if_else(evc$XGROUPNAME == "Lowland Forests", "Dry Forests", evc$XGROUPNAME)
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Herb-rich Woodlands", "Dry Forests", evc$XGROUPNAME)

# make snakey lil ones NA - we will interpolate these values from surrounds  
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Riparian Scrubs or Swampy Scrubs and Woodlands", NA, evc$XGROUPNAME)
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Coastal Scrubs Grasslands and Woodlands", NA, evc$XGROUPNAME)
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Wetlands", NA, evc$XGROUPNAME)
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Riverine Grassy Woodlands or Forests", NA, evc$XGROUPNAME)
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Plains Woodlands or Forests", NA, evc$XGROUPNAME)
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Rocky Outcrop or Escarpment Scrubs", NA, evc$XGROUPNAME)
evc$XGROUPNAME <- ifelse(evc$XGROUPNAME == "Salt-tolerant and/or succulent Shrublands", NA, evc$XGROUPNAME)
table(evc$XGROUPNAME)
# subset to just groupname
evc <- subset(evc, select = c(XGROUPNAME))
# do some renaming
evc <- rename(evc, vegetation = XGROUPNAME)
evc$vegetation <- gsub("Wet or Damp Forests", "Wet Forests", evc$vegetation)
unique(evc$vegetation)

# Otways only has a smattering of dry forest and heathy woodlands - make a new file for otways and condense these
evc_otways <- evc
evc_otways$vegetation <- ifelse(evc_otways$vegetation == "Heathy Woodlands", "Dry Forests", evc_otways$vegetation)
evc_otways$vegetation <- ifelse(evc_otways$vegetation == "Lowland Forests", "Dry Forests", evc_otways$vegetation)
unique(evc_otways$vegetation)

# make evc SpatialPolygonsDataFrame for secr
evc <- as_Spatial(evc)
evc_otways <- as_Spatial(evc_otways)

# add covariates to secr mask file
# load capthists
mrch_a    <- readRDS("derived_data/capthists/mrch_a.RData")
mrch_c    <- readRDS("derived_data/capthists/mrch_c.RData")
mrch_h    <- readRDS("derived_data/capthists/mrch_h.RData")
mrch_mc   <- readRDS( "derived_data/capthists/mrch_mc.RData")
mrch_s_17 <- readRDS("derived_data/capthists/mrch_s_17.RData")
mrch_n_17 <- readRDS("derived_data/capthists/mrch_n_17.RData")
mrch_s_18 <- readRDS("derived_data/capthists/mrch_s_18.RData")
mrch_n_18 <- readRDS("derived_data/capthists/mrch_n_18.RData")
mrch_s_19 <- readRDS("derived_data/capthists/mrch_s_19.RData")
mrch_n_19 <- readRDS("derived_data/capthists/mrch_n_19.RData")

# make capthists for lower glenelg blocks
ch_lgnpn <- read.capthist(captfile = "raw_data/LGNPN_caphist.txt", trapfile = "raw_data/LGNPN_traps.txt", detector = "proximity")
ch_lgnps <- read.capthist(captfile = "raw_data/LGNPS_caphist.txt", trapfile = "raw_data/LGNPS_traps.txt", detector = "proximity")

# combine
mrch <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc, ch_lgnpn, ch_lgnps, mrch_n_17, mrch_s_17, mrch_n_18, mrch_s_18, mrch_n_19, mrch_s_19)

# but need to do is seperately per session
mask_a  <- make.mask(traps(mrch_a),  buffer = 4000, spacing = 200, type = 'trapbuffer')     
mask_c  <- make.mask(traps(mrch_c),  buffer = 4000, spacing = 200, type = 'trapbuffer')     
mask_h  <- make.mask(traps(mrch_h),  buffer = 4000, spacing = 200, type = 'trapbuffer')     
mask_mc <- make.mask(traps(mrch_mc), buffer = 4000, spacing = 200, type = 'trapbuffer') 
mask_n_17 <- make.mask(traps(mrch_n_17), buffer = 4000, spacing = 200, type = 'trapbuffer')     
mask_n_18 <- make.mask(traps(mrch_n_18), buffer = 4000, spacing = 200, type = 'trapbuffer')     
mask_n_19 <- make.mask(traps(mrch_n_19), buffer = 4000, spacing = 200, type = 'trapbuffer')   
mask_s_17 <- make.mask(traps(mrch_s_17), buffer = 4000, spacing = 200, type = 'trapbuffer')     
mask_s_18 <- make.mask(traps(mrch_s_18), buffer = 4000, spacing = 200, type = 'trapbuffer')     
mask_s_19 <- make.mask(traps(mrch_s_19), buffer = 4000, spacing = 200, type = 'trapbuffer')  

# make fancy masks which account for glenelg river for lower glenelg sites- exclude areas south / north of glenelg river - cats can't cross
poly_north <- readOGR(dsn='raw_data/lower_glenelg_shp', layer='LGNPN_mask')
poly_south <- readOGR(dsn='raw_data/lower_glenelg_shp', layer='LGNPS_mask')
mask_lgnpn <- make.mask(traps(ch_lgnpn), spacing = 200, buffer = 4000, type='trapbuffer', poly = poly_north, poly.habitat	
                        = TRUE, keep.poly = FALSE) 
mask_lgnps <- make.mask(traps(ch_lgnps), spacing = 200, buffer = 4000, type ='trapbuffer', poly = poly_south, poly.habitat	
                        = TRUE, keep.poly = FALSE)

# add covariates to each mask
mask_a     <- addCovariates(mask_a, evc)
mask_c     <- addCovariates(mask_c, evc)
mask_h     <- addCovariates(mask_h, evc)
mask_mc    <- addCovariates(mask_mc, evc)
mask_lgnpn <- addCovariates(mask_lgnpn, evc)
mask_lgnps <- addCovariates(mask_lgnps, evc)
mask_n_17  <- addCovariates(mask_n_17, evc_otways)
mask_n_18  <- addCovariates(mask_n_18, evc_otways)
mask_n_19  <- addCovariates(mask_n_19, evc_otways)
mask_s_17  <- addCovariates(mask_s_17, evc_otways)
mask_s_18  <- addCovariates(mask_s_18, evc_otways)
mask_s_19  <- addCovariates(mask_s_19, evc_otways)



# FIX NA VALUES (extrapolate category from nearest cell)
# function to do so - from Efford secr tutorial on habitat masks
copynearest <- function(mask, covariate){
  NAcov <- is.na(covariates(mask)[,covariate])
  OK <- subset(mask, !NAcov)
  i <- nearesttrap(mask, OK)
  covariates(mask)[,covariate][NAcov] <- covariates(OK)[i[NAcov],covariate] 
  mask
}

mask_a     <- copynearest(mask_a,     'vegetation')
mask_c     <- copynearest(mask_c,     'vegetation')
mask_h     <- copynearest(mask_h,     'vegetation')
mask_mc    <- copynearest(mask_mc,    'vegetation')
mask_lgnpn <- copynearest(mask_lgnpn, 'vegetation')
mask_lgnps <- copynearest(mask_lgnps, 'vegetation')
mask_n_17  <- copynearest(mask_n_17,  'vegetation')
mask_n_18  <- copynearest(mask_n_18,  'vegetation')
mask_n_19  <- copynearest(mask_n_19,  'vegetation')
mask_s_17  <- copynearest(mask_s_17,  'vegetation')
mask_s_18  <- copynearest(mask_s_18,  'vegetation')
mask_s_19  <- copynearest(mask_s_19,  'vegetation')


# FOX GAMS ----------------------------------------------------------------
# load fox presence absence records
records <- read.csv("raw_data/spp_records_pa.csv")
# lower glenelg:
records_lg <- read.csv("raw_data/lower_glenelg_fox_pa.csv")
records_lg$station_year <- records_lg$station
# merge together
records <- bind_rows(records, records_lg)

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
gam_o_fox <- bam(fox ~ year + s(x, y, by = year, bs = "ds", m = c(1, 0.5), k = 100) + 
                       s(station, bs = "re") + 
                       offset(log(survey_duration)), 
                 data = records_otways, family = binomial, discrete = TRUE)   
plot(gam_o_fox, pages = 1, scheme = 2)
summary(gam_o_fox)



# GLENELG PREDICTIONS -----------------------------------------------------
# PREDICT GAM INTO HABITAT MASK -------------------------------------------
# make multisession files
mrch_glenelg <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc, ch_lgnpn, ch_lgnps)
# glenelg masks
#glenelg_mask <- list(mask_a, mask_c, mask_h, mask_mc, mask_lgnpn, mask_lgnps)
glenelg_mask <- make.mask(traps(mrch_glenelg), buffer = 4000, spacing = 200, type = 'trapbuffer')       

# make a new dataframe with each habitat mask cell so we can add covariates from the GAMs
glenelg_mask_df <- do.call(rbind.data.frame, glenelg_mask)
#  first give row values a column name
glenelg_mask_df <- setNames(cbind(rownames(glenelg_mask_df), glenelg_mask_df, row.names = NULL), c("sess", "x", "y"))
# make sess column the same per session
glenelg_mask_df$sess <- as.character(glenelg_mask_df$sess)
glenelg_mask_df$sess <- gsub('[0-9]+', '', glenelg_mask_df$sess)
glenelg_mask_df$sess <- gsub('[.]+', '', glenelg_mask_df$sess)
head(glenelg_mask_df)

# now do the same with lower glenelg
lgnpn_mask_df <- as.data.frame(mask_lgnpn)
lgnps_mask_df <- as.data.frame(mask_lgnps)
lgnpn_mask_df$sess <- "north"
lgnps_mask_df$sess <- "south"
lgnp_mask_df <- rbind(lgnpn_mask_df, lgnps_mask_df)
lgnp_mask_df <- relocate(lgnp_mask_df, sess, x, y)
head(lgnp_mask_df)

# add to glenelg_mask_df
glenelg_mask_df <- rbind(glenelg_mask_df, lgnp_mask_df)

# add survey duration dummy variable
glenelg_mask_df$survey_duration <- 60 

# predict GAM into the dataframe 
glenelg_mask_df <- cbind(glenelg_mask_df, predict.gam(gam_g_fox, newdata = glenelg_mask_df, se.fit = TRUE, type = "link"))
# rename fox_predicted to fox_predicted_trapcov for model specification
colnames(glenelg_mask_df) <- c("sess", "x", "y", "survey_duration", "fox_predicted", "fox_predicted_se")

# because we have a multi-session CH - we need to add in covariates separately for each session then merge back into multi-session mask...
# split dataframe by session 
glenelg_mask_df_a   <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_a"),]
glenelg_mask_df_c   <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_c"),]
glenelg_mask_df_h   <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_h"),]
glenelg_mask_df_mc  <- glenelg_mask_df[which(glenelg_mask_df$sess == "mrch_mc"),]
glenelg_mask_df_lgnpn  <- glenelg_mask_df[which(glenelg_mask_df$sess == "north"),]
glenelg_mask_df_lgnps  <- glenelg_mask_df[which(glenelg_mask_df$sess == "south"),]

# add covariates to each session mask separately 
covariates(mask_a)  <- cbind(glenelg_mask_df_a[,5:6],  covariates(mask_a))
covariates(mask_c)  <- cbind(glenelg_mask_df_c[,5:6],  covariates(mask_c))
covariates(mask_h)  <- cbind(glenelg_mask_df_h[,5:6],  covariates(mask_h))
covariates(mask_mc) <- cbind(glenelg_mask_df_mc[,5:6], covariates(mask_mc))
covariates(mask_lgnpn) <- cbind(glenelg_mask_df_lgnpn[,5:6], covariates(mask_lgnpn))
covariates(mask_lgnps) <- cbind(glenelg_mask_df_lgnps[,5:6], covariates(mask_lgnps))


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
glenelg_traps_df$station <- gsub('ch_lgnpn.', "", glenelg_traps_df$station)
glenelg_traps_df$station <- gsub('ch_lgnps.', "", glenelg_traps_df$station)

# add survey duration dummy variable
glenelg_traps_df$survey_duration <- 60 

# predict GAM model values into the traps dataframe
# predict into the dataframe - exclude fox predictor in the cat prey model to avoid double counting this effect
glenelg_traps_df <- cbind(glenelg_traps_df, predict.gam(gam_g_fox, newdata = glenelg_traps_df, se.fit = TRUE, type = "link"))

# rename fox_predicted to fox_predicted_trapcov for model specification
colnames(glenelg_traps_df) <- c("station", "x", "y", "survey_duration", "fox_predicted_trapcov", "fox_predicted_trapcov_se")

# need to split dataframes and capthists into single sessions to add covariates 
glenelg_traps_df_a   <- glenelg_traps_df[which(grepl("^A", glenelg_traps_df$station)),]
glenelg_traps_df_c   <- glenelg_traps_df[which(grepl("^C", glenelg_traps_df$station)),]
glenelg_traps_df_h   <- glenelg_traps_df[which(grepl("^H", glenelg_traps_df$station)),]
glenelg_traps_df_mc  <- glenelg_traps_df[which(grepl("^M", glenelg_traps_df$station)),]
# lower glenelg:
glenelg_traps_df_lgnp  <- glenelg_traps_df[which(grepl("^L", glenelg_traps_df$station)),]
glenelg_traps_df_lgnp$st_no <- glenelg_traps_df_lgnp$station
glenelg_traps_df_lgnp$st_no <- as.integer(gsub('LGNP', "", glenelg_traps_df_lgnp$st_no))
glenelg_traps_df_lgnpn  <- filter(glenelg_traps_df_lgnp, st_no >= 65)
glenelg_traps_df_lgnps  <- filter(glenelg_traps_df_lgnp, st_no < 65)

# add values into capthists - fox
covariates(traps(mrch_a)) <- glenelg_traps_df_a[,5:6]
covariates(traps(mrch_c)) <- glenelg_traps_df_c[,5:6]
covariates(traps(mrch_h)) <- glenelg_traps_df_h[,5:6]
covariates(traps(mrch_mc)) <- glenelg_traps_df_mc[,5:6]
covariates(traps(ch_lgnpn)) <- glenelg_traps_df_lgnpn[,5:6]
covariates(traps(ch_lgnps)) <- glenelg_traps_df_lgnps[,5:6]


# OTWAY PREDICTIONS -----------------------------------------------------
# merge to single capthist (in order of deployment) and make mask
mrch_otways <- MS.capthist(mrch_n_17, mrch_s_17, mrch_n_18, mrch_s_18, mrch_n_19, mrch_s_19)
otways_mask <- make.mask(traps(mrch_otways), buffer = 4000, 200, type = 'trapbuffer')       

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
otways_mask_df <- cbind(otways_mask_df, predict.gam(gam_o_fox, newdata = otways_mask_df, se.fit = TRUE, type = "link", exclude = c("s(station)")))
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
covariates(mask_n_17) <- cbind(otways_mask_df_n_17[,7:8], covariates(mask_n_17))
covariates(mask_n_18) <- cbind(otways_mask_df_n_18[,7:8], covariates(mask_n_18))
covariates(mask_n_19) <- cbind(otways_mask_df_n_19[,7:8], covariates(mask_n_19))
covariates(mask_s_17) <- cbind(otways_mask_df_s_17[,7:8], covariates(mask_s_17))
covariates(mask_s_18) <- cbind(otways_mask_df_s_18[,7:8], covariates(mask_s_18))
covariates(mask_s_19) <- cbind(otways_mask_df_s_19[,7:8], covariates(mask_s_19))



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
otways_traps_df <- cbind(otways_traps_df, predict.gam(gam_o_fox, newdata = otways_traps_df, se.fit = TRUE, type = "link", exclude = c("s(station)")))
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



# COMBINE, CHECKS AND SAVE -------------------------------------------------------------
# merge masks and capthists back together: same order!!!
mrch <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc, ch_lgnpn, ch_lgnps, mrch_n_17, mrch_s_17, mrch_n_18, mrch_s_18, mrch_n_19, mrch_s_19)
mrch_glenelg <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc, ch_lgnpn, ch_lgnps)
mrch_otways <- MS.capthist(mrch_n_17, mrch_s_17, mrch_n_18, mrch_s_18, mrch_n_19, mrch_s_19)
masks <- list(mask_a, mask_c, mask_h, mask_mc, mask_lgnpn, mask_lgnps, mask_n_17, mask_s_17, mask_n_18, mask_s_18, mask_n_19, mask_s_19)
masks_glenelg <- list(mask_a, mask_c, mask_h, mask_mc, mask_lgnpn, mask_lgnps)
masks_otways <- list(mask_n_17, mask_s_17, mask_n_18, mask_s_18, mask_n_19, mask_s_19)



# plot Vegetation class - glenelg
forestcol <- terrain.colors(6, rev = TRUE)[c(1,3,5)]
names(masks) <- c("Glenelg: Annya", "Glenelg: Cobboboonee", "Glenelg: Hotspur", "Glenelg: Mt Clay", "Glenelg: LGNP north", "Glenelg: LGNP south", "Otways: North", "Otways: South", "Otways: North '18", "Otways: South '18", "Otways: North '19", "Otways: South '19")
par (mfrow = c(1, 7), mar = c(0,0,5,0))
for (sess in 1:6) {
  plot(masks[[sess]], covariate = 'vegetation', legend = FALSE, dots = FALSE, border = 100, col = forestcol) 
  plot(traps(mrch)[[sess]], add = TRUE, detpar = list(pch = 16, cex = 0.25, col = "black"))
  mtext(side=3, paste(session(masks)[[sess]])) }
# null plot for seperate legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend (x = "left", fill = forestcol, legend = levels(factor(evc$vegetation))[1:3], bty='n', title = "Vegetation type", pt.cex=3, cex=1.5)


# plot Vegetation class - Otways
forestcol <- terrain.colors(6, rev = TRUE)[c(1,3,6)]
names(masks) <- c("Glenelg: Annya", "Glenelg: Cobboboonee", "Glenelg: Hotspur", "Glenelg: Mt Clay", "Glenelg: LGNP north", "Glenelg: LGNP south", "Otways: North", "Otways: South", "Otways: North '18", "Otways: South '18", "Otways: North '19", "Otways: South '19")
par (mfrow = c(1,3), mar = c(0,0,5,0))
for (sess in 7:8) {
  plot(masks[[sess]], covariate = 'vegetation', legend = FALSE, dots = FALSE, border = 100, col = forestcol) 
  plot(traps(mrch)[[sess]], add = TRUE, detpar = list(pch = 16, cex = 0.5, col = "black"))
  mtext(side=3, paste(session(masks)[[sess]])) }
# null plot for seperate legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend (x = "left", fill = forestcol, legend = levels(factor(evc_otways$vegetation))[1:3], bty='n', title = "Vegetation type", pt.cex=3, cex=1.5)

# plot foxes - glenelg
par (mfrow = c(3,2), mar = c(2,2,4,5))
names(masks) <- c("Annya", "Cobboboonee", "Hotspur", "Mt Clay", "LGNP north", "LGNP south", "North '17", "South '17", "North '18", "South '18", "North '19", "South '19")
for (sess in 1:6) {
  plot(masks[[sess]], covariate = 'fox_predicted', legend = TRUE, dots = FALSE, border = 100)
  plot(traps(mrch)[[sess]], add = TRUE, detpar = list(pch = 16, cex = 0.25, col = "black"))
  mtext(side=3, paste(session(masks)[[sess]])) }

# plot foxes - otways
par (mfrow = c(3,2), mar = c(1,1,3,1))
names(masks) <- c("Annya", "Cobboboonee", "Hotspur", "Mt Clay", "LGNP north", "LGNP south", "North '17", "South '17", "North '18", "South '18", "North '19", "South '19")
for (sess in 7:12) {
  plot(masks[[sess]], covariate = 'fox_predicted', legend = TRUE, dots = FALSE, border = 100)
  plot(traps(mrch)[[sess]], add = TRUE, detpar = list(pch = 16, cex = 0.25, col = "black"))
  mtext(side=3, paste(session(masks)[[sess]])) }


## SAVE EVERYTHING
# save GAMs
saveRDS(gam_g_fox, "models/gam_g_fox.RData")  
saveRDS(gam_o_fox, "models/gam_o_fox.RData")  
# save mask dataframes  (useful for plotting etc in this format)
saveRDS(glenelg_mask_df, "derived_data/glenelg_mask_df.RData")
saveRDS(otways_mask_df, "derived_data/otway_mask_df.RData")
# save capthists
saveRDS(mrch, "derived_data/mrch.RData")
saveRDS(mrch_glenelg, "derived_data/mrch_glenelg.RData")
saveRDS(mrch_otways, "derived_data/mrch_otways.RData")
# save masks
saveRDS(masks, "derived_data/masks.RData")
saveRDS(masks_glenelg, "derived_data/masks_glenelg.RData")
saveRDS(masks_otways, "derived_data/masks_otways.RData")

## END