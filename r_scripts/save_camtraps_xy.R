
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

x <- as.data.frame(traps(mrch))
mask_df <- do.call(rbind.data.frame, traps(mrch))
mask_df <- setNames(cbind(rownames(mask_df), mask_df, row.names = NULL), c("station", "x", "y"))

library(tidyverse)
mask_df <- mask_df %>% separate(station, into = c("block", "station"), sep = '[.]')
head(mask_df)


write.csv(mask_df, "cat_density_cams_xy.csv")

