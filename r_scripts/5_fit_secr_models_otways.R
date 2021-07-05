# RUN FERAL CAT SECR MODELS - GLENELG REGION 2017-19 
# Matthew Rees

# Set-up 
rm(list = ls())
options(scipen = 999)
library(secr)  
library(dplyr)

# load files
mrch_otways  <- readRDS("derived_data/mrch_otways.RData")
masks_otways <- readRDS("derived_data/masks_otways.RData")


## specify session covariates
sesscov <- data.frame(grid = factor(c("north", "south", "north", "south", "north", "south")), 
                      year = factor(c("2017", "2017", "2018", "2018", "2019", "2019")), 
                      baiting1 = factor(c("0ub", "0ub", "0ub", "1b", "0ub", "1b")),
                      baiting2 = factor(c("0ub", "0ub", "0ub", "1b", "0ub", "2b")))
sesscov




# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 0, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1), method = "Nelder-Mead")
fit1_ex <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_ex)
saveRDS(df_fits, "models/otways_df_fits.RData")
df_fits <- readRDS("models/otways_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex wins --> Use this detector function for all subsequent model fitting


# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, nsim = 10000), method = "Nelder-Mead", fixed = list(pID = 1), start = fit1_ex)
saveRDS(fit1_adj,   "models/otways_fit1_adj.RData")
fit1_adj <- readRDS("models/otways_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#5.713383 1.000000 


# 3) Fit actual models
## Null models
fit_null <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ 1, sigma ~ 1))
fit_null_Dveg <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + vegetation, g0 ~ 1, sigma ~ 1))
fit_null_g0T <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ T, sigma ~ 1))
AIC(fit_null, fit_null_Dveg, fit_null_g0T, criterion = "AICc")[,-2] 

## 4) Baiting models
fit_sess1 <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ 1, sigma ~ 1)) 
fit_sess2 <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ baiting2, sigma ~ baiting2)) 
fit_sess3 <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ year, sigma ~ year)) 
fit_sess4 <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ year + baiting2, sigma ~ year + baiting2)) 
AIC(fit_sess1, fit_sess2, fit_sess3, fit_sess4, criterion = "AICc")[,-2] 


## 5) Correlative models with fox occupancy
fit_fox_det <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
fit_fox_D <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + fox_predicted, g0 ~ 1, sigma ~ 1)) 
fit_fox_Ddet <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
fit_nl_fox_det <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ s(fox_predicted_trapcov, k = 3), sigma ~ s(fox_predicted_trapcov, k = 3))) 
fit_nl_fox_D <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + s(fox_predicted, k = 3), g0 ~ 1, sigma ~ 1)) 
fit_nl_fox_Ddet <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + s(fox_predicted, k = 3), g0 ~ s(fox_predicted_trapcov, k = 3), sigma ~ s(fox_predicted_trapcov, k = 3))) 
AIC(fit_null, fit_fox_det, fit_fox_D, fit_fox_Ddet, fit_nl_fox_det, fit_nl_fox_D, fit_nl_fox_Ddet, criterion = "AICc")[,-2] 


## 6) Combine models, save and compare AICc
otways_fits <- secrlist(fit_null, fit_null_Dveg, fit_null_g0T, fit_sess1, fit_sess2, fit_sess3, fit_sess4, fit_fox_det, fit_fox_D, fit_fox_Ddet, fit_nl_fox_det, fit_nl_fox_D, fit_nl_fox_Ddet)
saveRDS(otways_fits,   "models/otways_fits.RData")
otways_fits <- readRDS("models/otways_fits.RData")
AIC(otways_fits, criterion = "AICc")[,-2] 


# END