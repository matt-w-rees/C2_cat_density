# RUN FERAL CAT SECR MODELS - GLENELG REGION 2018 
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
sesscov <- data.frame(grid = factor(c("south", "north", "south", "north", "south", "north")), 
                      year = factor(c("2017", "2017", "2018", "2018", "2019", "2019")), 
                      before_after = factor(c("b", "b", "a", "a", "a", "a")))




# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 0, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE), fixed = list(pID = 1), method = "Nelder-Mead")
fit1_ex <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_ex)
saveRDS(df_fits, "models/otways_df_fits.RData")
#df_fits <- readRDS("models/otways_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex wins --> Use this detector function for all subsequent model fitting


# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, nsim = 10000), method = "Nelder-Mead", fixed = list(pID = 1), start = fit1_ex)
saveRDS(fit1_adj,   "models/otways_fit1_adj.RData")
#fit1_adj <- readRDS("models/otways_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#3.843615 1.000000 


# 3) Fit actual models
## Null models
fit1      <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ 1, sigma ~ 1))
fit1_Dveg <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ 1, sigma ~ 1))
fit1_g0T  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ T, sigma ~ 1))
AIC(fit1, fit1_Dveg, fit1_g0T, criterion = "AICc")[,-2] 


## 4) Experimental models
# a) same density in the unbaited sites, same effect of baiting in the baited sites
# b) different density in the unbaited sites, same effect of baiting in the baited sites:
# c) different density in the unbaited sites, different effect of baiting in the baited sites (i.e. D ~ session)
fit_BACIa_D     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid + before_after, g0 ~ 1, sigma ~ 1)) 
fit_BACIa_Ddet  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid + before_after, g0 ~ grid + before_after, sigma ~ grid + before_after)) 
fit_BACIb_D     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid * before_after, g0 ~ 1, sigma ~ 1)) 
fit_BACIb_Ddet  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid * before_after, g0 ~ grid * before_after, sigma ~ grid * before_after)) 
fit_BACIc_D     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid + year, g0 ~ , sigma ~ 1)) 
fit_BACIc_Ddet  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid + year, g0 ~ grid + year, sigma ~ grid + year)) 
fit_BACId_D     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid * year, g0 ~ 1, sigma ~ 1)) 
fit_BACId_Ddet  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid * year, g0 ~ grid * year, sigma ~ grid * year)) 
AIC(fit1, fit_BACIa_D, fit_BACIa_Ddet, fit_BACIb_D, fit_BACIb_Ddet, fit_BACIc_D, fit_BACIc_Ddet, fit_BACId_D, fit_BACId_Ddet, criterion = "AICc")[,-2] 


## 5) Correlative models with fox occupancy
fit_fox_D <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1, sigma ~ 1)) 
fit_fox_Dd <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
AIC(fit1, fit_fox_D, fit_fox_Dd, criterion = "AICc")[,-2] 


## 6) Combine models, save and compare AICc
otways_fits <- secrlist(fit1, fit1_Dveg, fit1_g0T, fit_BACIa_D, fit_BACIa_Ddet, fit_BACIb_D, fit_BACIb_Ddet, fit_BACIc_D, fit_BACIc_Ddet, fit_BACId_D, fit_BACId_Ddet, fit_fox_D, fit_fox_Dd)
saveRDS(otways_fits,   "models/otways_fits.RData")
#otways_fits <- readRDS("models/otways_fits.RData")
AIC(otways_fits, criterion = "AICc")[,-2] 


# END




# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 0, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1))
fit1_hr <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
fit1_ex <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_hr, fit1_ex)
saveRDS(df_fits,   "models/otways_df_fits.RData")
df_fits <- readRDS("models/otways_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex --> Use this detector function for all subsequent model fitting
# also check buffer size is good
esa.plot(df_fits, max.buffer = 4500, ylim = c(0.001,0.003))
# also go back and double check masK spacing fits with this sigma estimate. 


# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_otways, mask = masks_otways, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, nsim = 10000), fixed = list(pID = 1), start = fit1_hr)
saveRDS(fit1_adj,  "models/otways_fit1_adj.RData")
fit1_adj <- readRDS("models/otways_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#3.853155 1.000000 


# 3) fit models
# first - before we consider foxes - choose best "null" model - cat density constant? depdends on vegetation? or vary over the years?
fit1a <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ 1, sigma ~ 1)) 
fit1b <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ 1, sigma ~ 1)) 
fit1c <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ 1, sigma ~ 1)) 
AIC(fit1a, fit1b, fit1c, criterion = "AICc")[,-2]
# no evidence of year or vegetation impact - compare next sets of models to fit1a

# run traditional BACI models on the fox control (averaged across the two after years) - with and without detectability impacts
fit_baci_null_D    <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after + grid, g0 ~ 1, sigma ~ 1)) 
fit_baci_null_Ddet <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after + grid, g0 ~ before_after + grid, sigma ~ before_after + grid)) 
fit_baci_D         <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after * grid, g0 ~ 1, sigma ~ 1)) 
fit_baci_Ddet      <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ before_after * grid, g0 ~ before_after * grid, sigma ~ before_after * grid)) 
AIC(fit_baci_null_D, fit_baci_null_Ddet, fit_baci_D, fit_baci_Ddet, criterion = "AICc")[,-2]

# run fox models - with and without detectability impacts
fit_foxbaited_D    <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ 1, sigma ~ 1)) 
fit_foxbaited_Ddet <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ foxbaited, sigma ~ foxbaited)) 
fit_fox_D          <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1, sigma ~ 1)) 
fit_fox_Ddet       <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
AIC(fit1a, fit_foxbaited_D, fit_foxbaited_Ddet, fit_fox_D, fit_fox_Ddet, criterion = "AICc")[,-2]

# put all models in a list, compare again, and save
otways_fits <- secrlist(fit1a, fit1b, fit1c, fit_baci_null_D, fit_baci_null_Ddet, fit_baci_D, fit_baci_Ddet,  fit_foxbaited_D, fit_foxbaited_Ddet, fit_fox_D, fit_fox_Ddet)
saveRDS(otways_fits, "models/otways_fits.RData")
AIC(otways_fits, criterion = "AICc")[,-2] 


# END