# RUN FERAL CAT SECR MODELS - GLENELG REGION 2018 
# Matthew Rees

# Set-up 
rm(list = ls())
options(scipen = 999)
library(secr)  
library(dplyr)

# Load secr data
masks_glenelg <- readRDS("derived_data/masks_glenelg.RData")
mrch_glenelg <- readRDS("derived_data/mrch_glenelg.RData")

# specify session covariates
sesscov <- data.frame(pair = factor(c("p1", "p1", "p2", "p2")),
                      foxbaiting_01  = factor(c("0unbaited", "1baited", "0unbaited", "1baited")),
                      foxbaiting_012 = factor(c("0unbaited", "1baited", "0unbaited", "2baited")))
sesscov


# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_glenelg, mask = masks_glenelg, detectfn = 0, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE), fixed = list(pID = 1), method = "Nelder-Mead")
fit1_ex <- secr.fit(mrch_glenelg, mask = masks_glenelg, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_ex)
saveRDS(df_fits, "models/glenelg_df_fits.RData")
df_fits <- readRDS("models/glenelg_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex wins --> Use this detector function for all subsequent model fitting


# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_glenelg, mask = masks_glenelg, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, nsim = 10000), method = "Nelder-Mead", fixed = list(pID = 1), start = df_fits$fit1_ex)
saveRDS(fit1_adj,   "models/glenelg_fit1_adj.RData")
fit1_adj <- readRDS("models/glenelg_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#3.843615 1.000000 


# 3) fit actual models
## Null models
# We chose the 4 grids based on similar landscape context and management. The main differences were the:
# a) relative proportion of vegetation type (cleared, lowland forest and heathy woodland), 
# b) and that one pair which was surveyed at a time (but 2nd pair was survey immediately after 1st) - this slight change in the time of year may have impacted cat detectability (but would not expect any numerical pop. changes in the two months).
# Was there an impact of vegetation on density (B)? or pair (slightly different time of the year) on detectability (C)? Compare these to the null (A). 
fit_null         <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1,          g0 ~ 1,    sigma ~ 1))
fit_null_Dveg    <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ 1,    sigma ~ 1))
fit_null_g0T     <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1,          g0 ~ T,    sigma ~ 1))
AIC(fit_null, fit_null_Dveg, fit_null_detpair, fit_null_g0T, criterion = "AICc")[,-2] 
# no model had an AICc any better than the null, so no need to carry any of these forward - and all good to compare baited treatments to both unbaited sites (they are basically the same forest anyway)

## Is there an effect of foxbaiting on density and/or detectability? Was this consistent across two baited sites (fit 2) or different (fit 3) - compare to null model. 
fit_sess1 <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ 1, sigma ~ 1)) 
fit_sess2 <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ foxbaiting_012, sigma ~ foxbaiting_012)) 
fit_sess3 <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ pair, sigma ~ pair)) 
fit_sess4 <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat, contrasts = list(session = MASS::contr.sdif)), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ pair + foxbaiting_012, sigma ~ pair + foxbaiting_012)) 
AIC(fit_null, fit_sess1, fit_sess2, fit_sess3, fit_sess4, criterion = "AICc")[,-2] 
  
## 5) Correlative models with fox occupancy 
fit_fox_det  <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1,             g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov))
fit_fox_D    <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1,                     sigma ~ 1)) 
fit_fox_Ddet <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
AIC(fit_null, fit_fox_det, fit_fox_D,  fit_fox_Ddet, criterion = "AICc")[,-2] 

## 6) Combine models, save and compare AICc altogether
glenelg_fits <- secrlist(fit_null, fit_null_Dveg, fit_null_g0T, fit_sess1, fit_sess2, fit_sess3, fit_sess4, fit_fox_det, fit_fox_D, fit_fox_Ddet)
saveRDS(glenelg_fits,   "models/glenelg_fits.RData")
glenelg_fits <- readRDS("models/glenelg_fits.RData")
AIC(glenelg_fits, criterion = "AICc")[,-2] 


# END 