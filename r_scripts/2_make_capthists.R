# make standard secr capthists - seperately per session
# then add unmarked cats to capture-recapture capture history
# instead of ignoring mark status unknown (i.e. unidentified tabby cats), include them with the black cats as unmarked category (mt clay - no black cats)

library(secr)

# CAPTHISTS PER SESSION ---------------------------------------------------
## GLENELG
ch_mc <- readRDS("derived_data/capthists/ch_mc.RData")
unid_mc <- as.matrix(read.csv("derived_data/capthists/unid_mc.csv")[, -1]) 
mrch_mc <- addSightings(ch_mc, unmarked = unid_mc)
summary(mrch_mc, terse = TRUE)
# hotspur 
ch_h <- readRDS("derived_data/capthists/ch_h.RData")
unid_h <- as.matrix(read.csv("derived_data/capthists/unid_h.csv")[, -1]) 
unm_h <- as.matrix(read.csv("derived_data/capthists/unm_h.csv")[, -1]) 
all_unid_h = unid_h + unm_h   
mrch_h <- addSightings(ch_h, unmarked = all_unid_h)
summary(mrch_h, terse = TRUE)
# annya 
ch_a <- readRDS("derived_data/capthists/ch_a.RData")
unid_a <- as.matrix(read.csv("derived_data/capthists/unid_a.csv")[, -1]) 
unm_a <- as.matrix(read.csv("derived_data/capthists/unm_a.csv")[, -1]) 
all_unid_a = unid_a + unm_a   
mrch_a <- addSightings(ch_a, unmarked = all_unid_a)
summary(mrch_a, terse = TRUE)
# cobbob 
ch_c <- readRDS("derived_data/capthists/ch_c.RData")
unid_c <- as.matrix(read.csv("derived_data/capthists/unid_c.csv")[, -1]) 
unm_c <- as.matrix(read.csv("derived_data/capthists/unm_c.csv")[, -1]) 
all_unid_c = unid_c + unm_c   
mrch_c <- addSightings(ch_c, unmarked = all_unid_c)
summary(mrch_c, terse = TRUE)

## OTWAYS
# northern grid 2017 
ch_n_17 <- readRDS("derived_data/capthists/ch_n_17.RData")
unid_n_17 <- as.matrix(read.csv("derived_data/capthists/unid_n_17.csv")[, -1]) 
unm_n_17 <- as.matrix(read.csv("derived_data/capthists/unm_n_17.csv")[, -1]) 
all_unid_n_17 = unid_n_17 + unm_n_17
mrch_n_17 <- addSightings(ch_n_17, unmarked = all_unid_n_17)
summary(mrch_n_17, terse = TRUE)
# northern grid 2018 
ch_n_18 <- readRDS("derived_data/capthists/ch_n_18.RData")
unid_n_18 <- as.matrix(read.csv("derived_data/capthists/unid_n_18.csv")[, -1]) 
unm_n_18 <- as.matrix(read.csv("derived_data/capthists/unm_n_18.csv")[, -1]) 
all_unid_n_18 = unid_n_18 + unm_n_18
mrch_n_18 <- addSightings(ch_n_18, unmarked = all_unid_n_18)
summary(mrch_n_18, terse = TRUE)
# northern grid 2019 
ch_n_19 <- readRDS("derived_data/capthists/ch_n_19.RData")
unid_n_19 <- as.matrix(read.csv("derived_data/capthists/unid_n_19.csv")[, -1]) 
unm_n_19 <- as.matrix(read.csv("derived_data/capthists/unm_n_19.csv")[, -1]) 
all_unid_n_19 = unid_n_19 + unm_n_19
mrch_n_19 <- addSightings(ch_n_19, unmarked = all_unid_n_19)
summary(mrch_n_19, terse = TRUE)
# southern grid 2017 
ch_s_17 <- readRDS("derived_data/capthists/ch_s_17.RData")
unid_s_17 <- as.matrix(read.csv("derived_data/capthists/unid_s_17.csv")[, -1]) 
unm_s_17 <- as.matrix(read.csv("derived_data/capthists/unm_s_17.csv")[, -1]) 
all_unid_s_17 = unid_s_17 + unm_s_17
mrch_s_17 <- addSightings(ch_s_17, unmarked = all_unid_s_17)
summary(mrch_s_17, terse = TRUE)
# southern grid 2018 
ch_s_18 <- readRDS("derived_data/capthists/ch_s_18.RData")
unid_s_18 <- as.matrix(read.csv("derived_data/capthists/unid_s_18.csv")[, -1]) 
unm_s_18 <- as.matrix(read.csv("derived_data/capthists/unm_s_18.csv")  [, -1]) 
all_unid_s_18 = unid_s_18 + unm_s_18
mrch_s_18 <- addSightings(ch_s_18, unmarked = all_unid_s_18)
summary(mrch_s_18, terse = TRUE)
# southern grid 2019 
ch_s_19 <- readRDS("derived_data/capthists/ch_s_19.RData")
unid_s_19 <- as.matrix(read.csv("derived_data/capthists/unid_s_19.csv")[, -1]) 
unm_s_19 <- as.matrix(read.csv("derived_data/capthists/unm_s_19.csv")[, -1]) 
all_unid_s_19 = unid_s_19 + unm_s_19
mrch_s_19 <- addSightings(ch_s_19, unmarked = all_unid_s_19)
summary(mrch_s_19, terse = TRUE)



# MERGE CAPTHISTS ---------------------------------------------------------
# merge to single capthist
# in order of region (glenelg, otways), then deployment start dates (south, north 2018, and so on)
mrch <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc, mrch_s_17, mrch_n_17, mrch_s_18, mrch_n_18, mrch_s_19, mrch_n_19)
saveRDS(mrch, "derived_data/capthists/mrch.RData")
summary(mrch, terse = TRUE)

# END