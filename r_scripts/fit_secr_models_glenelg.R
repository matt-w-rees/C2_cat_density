# RUN FERAL CAT SECR MODELS - GLENELG REGION 2018 
# Matthew Rees

rm(list = ls())
options(scipen = 999)

library(secr)  
library(dplyr)

# load secr data
masks_glenelg <- readRDS("derived_data/masks_glenelg.RData")
mrch_glenelg <- readRDS("derived_data/mrch_glenelg.RData")

## specify session covariates
# get mean fox occ estimates
glenelg_mask_df <-  readRDS("derived_data/glenelg_mask_df.RData")
x <- glenelg_mask_df %>% 
  group_by(sess) %>% 
  mutate(fox_predicted_mean = mean(fox_predicted))
x <- round(unique(x$fox_predicted_mean), digits = 2)
x

# make session cov 
sesscov <- data.frame(fox_predicted_mean = x, foxbaited = factor(c("0unbaited", "1baited", "0unbaited", "1baited")))
str(sesscov)

# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_glenelg, mask = masks_glenelg, detectfn = 0, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE), fixed = list(pID = 1))
fit1_hr <- secr.fit(mrch_glenelg, mask = masks_glenelg, detectfn = 1, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
fit1_ex <- secr.fit(mrch_glenelg, mask = masks_glenelg, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_hr, fit1_ex)
saveRDS(df_fits, "models/glenelg_df_fits.RData")
df_fits <- readRDS("models/glenelg_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_hr wins --> Use this detector function for all subsequent model fitting
# also check buffer size is good
esa.plot(df_fits, max.buffer = 4000, ylim = c(0.001,0.003))


# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_glenelg, mask = masks_glenelg, detectfn = 1, trace = TRUE, ncores = 3, details = list(knownmarks = FALSE, nsim = 10000), fixed = list(pID = 1), method = "none", start = fit1_hr)
saveRDS(fit1_adj, "models/glenelg_fit1_adj.RData")
fit1_adj <- readRDS("models/glenelg_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#2.869429 1.000000 


# 3) fit actual models
# firstly, do we need the veg covariate?
fit1     <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, method = "none", model = list(D ~ 1, g0 ~ 1, sigma ~ 1))
fit1_veg <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, method = "none", model = list(D ~ vegetation, g0 ~ 1, sigma ~ 1))
AIC(fit1, fit1_veg, criterion = "AICc")[,-2]
# nah 

# run fox models - with and without detectability impacts
fit_fox_baited_D    <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ 1, sigma ~ 1)) 
fit_fox_baited_Ddet <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ foxbaited, sigma ~ foxbaited)) 
fit_fox_avg_D       <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean)) 
fit_fox_avg_Ddet    <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean)) 
fit_fox_fs_D        <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1, sigma ~ 1)) 
fit_fox_fs_Ddet     <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), method = "none", fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
AIC(fit1, fit_fox_baited_D, fit_fox_baited_Ddet, fit_fox_avg_D, fit_fox_avg_Ddet, fit_fox_fs_D, fit_fox_fs_Ddet, criterion = "AICc")[,-2]

# put all models in a list, compare again, and save
glenelg_fits <- secrlist(fit1, fit1_veg, fit_fox_baited_D, fit_fox_baited_Ddet, fit_fox_avg_D, fit_fox_avg_Ddet, fit_fox_fs_D, fit_fox_fs_Ddet)
saveRDS(glenelg_fits, "models/glenelg_fits.RData")
AIC(glenelg_fits, criterion = "AICc")[,-2] 

# END