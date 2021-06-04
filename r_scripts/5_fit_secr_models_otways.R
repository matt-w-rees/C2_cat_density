# RUN FERAL CAT SECR MODELS - OTWAY RANGES 

rm(list = ls())
options(scipen = 999)

library(secr)  # version 4.3.3
library(dplyr)

# load files
mrch_otways <- readRDS("derived_data/mrch_otways.RData")
masks_otways <- readRDS("derived_data/masks_otways.RData")


## specify session covariates
# get mean fox occ estimates
otways_mask_df <-  readRDS("derived_data/otway_mask_df.RData")
x <- otways_mask_df %>% 
  group_by(sess) %>% 
  mutate(fox_predicted_mean = mean(fox_predicted))
x <- round(unique(x$fox_predicted_mean), digits = 2)
x


## specify session covariates
sesscov <- data.frame(fox_predicted_mean = x,
                      foxbaited = factor(c("0unbaited", "0unbaited", "1baited", "0unbaited", "1baited", "0unbaited")), 
                      grid = factor(c("south", "north", "south", "north", "south", "north")), 
                      year = factor(c("2017", "2017", "2018", "2018", "2019", "2019")), 
                      before_after = factor(c("b", "b", "a", "a", "a", "a")))


# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 0, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1))
fit1_hr <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
fit1_ex <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_hr, fit1_ex)
saveRDS(df_fits, "models/otways_df_fits.RData")
df_fits <- readRDS("models/otways_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex --> Use this detector function for all subsequent model fitting
# also check buffer size is good
esa.plot(df_fits, max.buffer = 4500, ylim = c(0.001,0.003))
# also go back and double check maks spacing fits with this sigma estimate. 


# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 1, trace = TRUE, ncores = 6, details = list(knownmarks = FALSE, nsim = 10000), fixed = list(pID = 1), start = fit1_hr)
saveRDS(fit1_adj, "models/otways_fit1_adj.RData")
fit1_adj <- readRDS("models/otways_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#3.853155 1.000000 


# 3) fit models
# first - before we consider foxes - choose best "null" model - cat density constant? depdends on vegetation? or vary over the years?
fit1a <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ 1, sigma ~ 1)) 
fit1b <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ 1, sigma ~ 1)) 
fit1c <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ 1, sigma ~ 1)) 
AIC(fit1a, fit1b, fit1c, criterion = "AICc")[,-2]
# no evidence of year or vegetation impact - compare next sets of models to fit1a

# run traditional BACI models - with and without detectability impacts
fit_baci_null_D    <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after + grid, g0 ~ 1, sigma ~ 1)) 
fit_baci_null_Ddet <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after + grid, g0 ~ before_after + grid, sigma ~ before_after + grid)) 
fit_baci_D         <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after * grid, g0 ~ 1, sigma ~ 1)) 
fit_baci_Ddet      <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after * grid, g0 ~ before_after * grid, sigma ~ before_after * grid)) 
AIC(fit_baci_null_D, fit_baci_null_Ddet, fit_baci_D, fit_baci_Ddet, criterion = "AICc")[,-2]

# run fox models - with and without detectability impacts
fit_fox_baited_D    <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ 1, sigma ~ 1)) 
fit_fox_baited_Ddet <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ foxbaited, sigma ~ foxbaited)) 
fit_fox_avg_D       <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean)) 
fit_fox_avg_Ddet    <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean)) 
fit_fox_fs_D        <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1, sigma ~ 1)) 
fit_fox_fs_Ddet     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
AIC(fit1a, fit_fox_baited_D, fit_fox_baited_Ddet, fit_fox_avg_D, fit_fox_avg_Ddet, fit_fox_fs_D, fit_fox_fs_Ddet, criterion = "AICc")[,-2]

# put all models in a list, compare again, and save
otways_fits <- secrlist(fit1a, fit1b, fit1c, fit_baci_null_D, fit_baci_null_Ddet, fit_baci_D, fit_baci_Ddet, fit1a, fit_fox_baited_D, fit_fox_baited_Ddet, fit_fox_avg_D, fit_fox_avg_Ddet, fit_fox_fs_D, fit_fox_fs_Ddet)
saveRDS(otways_fits, "models/otways_fits.RData")
AIC(otways_fits, criterion = "AICc")[,-2] 


# END