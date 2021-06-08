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
                      foxbaiting = factor(c("0unbaited", "1baited", "0unbaited", "1baited")))



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
fit1      <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ 1, sigma ~ 1))
fit1_Dveg <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ 1, sigma ~ 1))
fit1_g0T  <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ T, sigma ~ 1))
AIC(fit1, fit1_Dveg, fit1_g0T, criterion = "AICc")[,-2] 


## 4) Experimental models
# a) same density in the unbaited sites, same effect of baiting in the baited sites
# b) different density in the unbaited sites, same effect of baiting in the baited sites:
# c) different density in the unbaited sites, different effect of baiting in the baited sites (i.e. D ~ session)
fit_CIa_D  <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaiting, g0 ~ 1, sigma ~ 1)) 
fit_CIa_Dd <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaiting, g0 ~ foxbaiting, sigma ~ foxbaiting)) 
fit_CIb_D  <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair + foxbaiting, g0 ~ 1, sigma ~ 1)) 
fit_CIb_Dd <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair + foxbaiting, g0 ~ pair + foxbaiting, sigma ~ pair + foxbaiting)) 
fit_CIc_D  <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair * foxbaiting, g0 ~ 1, sigma ~ 1)) 
fit_CIc_Dd <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair * foxbaiting, g0 ~ pair * foxbaiting, sigma ~ pair * foxbaiting)) 
AIC(fit1, fit_CIa_D , fit_CIa_Dd, fit_CIb_D, fit_CIb_Dd, fit_CIc_D, fit_CIc_Dd, criterion = "AICc")[,-2] 


## 5) Correlative models with fox occupancy
fit_fox_D <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1, sigma ~ 1)) 
fit_fox_Dd <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
AIC(fit1, fit_fox_D, fit_fox_Dd, criterion = "AICc")[,-2] 


## 6) Combine models, save and compare AICc
glenelg_fits <- secrlist(fit1, fit1_Dveg, fit1_g0T, fit_CIa_D , fit_CIa_Dd, fit_CIb_D, fit_CIb_Dd, fit_CIc_D, fit_CIc_Dd, fit_fox_D, fit_fox_Dd)
saveRDS(glenelg_fits,   "models/glenelg_fits.RData")
#glenelg_fits <- readRDS("models/glenelg_fits.RData")
AIC(glenelg_fits, criterion = "AICc")[,-2] 


# END