# RUN FERAL CAT SECR MODELS - otways REGION 2018 
# Matthew Rees

rm(list = ls())
options(scipen = 999)

# load packages
library(secr)  # version 4.3.3
library(mgcv)  
library(ggplot2)
library(patchwork)
library(dplyr)

# load secr data
mask_otways <- readRDS("derived_data/mask_otways.RData")
mrch_otways <- readRDS("derived_data/mrch_otways.RData")

# load models 
#df_fits <- readRDS("models/otways_df_fits.RData")
#fit1_adj <- readRDS("models/otways_fit1_adj.RData")
#otways_fits <- readRDS("models/otways_fits.RData")

## specify session covariates
# get mean fox occ estimates
otways_mask_df <-  readRDS("derived_data/otway_mask_df.RData")
x <- otways_mask_df %>% 
  group_by(sess) %>% 
  mutate(fox_predicted_mean = mean(fox_predicted))
x <- round(unique(x$fox_predicted_mean), digits = 2)

# make session cov 
sesscov <- data.frame(fox_predicted_mean = x, foxbaited = factor(c("0unbaited", "0unbaited", "1baited", "0unbaited", "1baited", "0unbaited")), treament = factor(c("t", "nt", "t", "nt", "t", "nt")), year = factor(c("2017", "2017", "2018", "2018", "2019", "2019")), before_after = factor(c("b", "b", "a", "a", "a", "a")))


# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_otways, mask = mask_otways, detectfn = 0, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1))
fit1_ex <- secr.fit(mrch_otways, mask = mask_otways, detectfn = 2, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_ex)
saveRDS(df_fits, "models/otways_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex --> Use this detector function for all subsequent model fitting
# also check buffer size is good
esa.plot(df_fits, max.buffer = 4500, ylim = c(0.001,0.003))
# also go back and double check maks spacing fits with this sigma estimate. 

# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_otways, mask = mask_otways, detectfn = 1, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, nsim = 10000), fixed = list(pID = 1))
saveRDS(fit1_adj, "models/otways_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#3.853155 1.000000 


# 3) fit actual models
fit1a <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ 1, sigma ~ 1)) 
#fit1b <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year, g0 ~ year, sigma ~ year)) 
fit2a <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + before_after + treatment, g0 ~ 1, sigma ~ 1)) 
#fit2b <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + before_after + treatment, g0 ~ year + before_after + treatment, sigma ~ year + before_after + treatment)) 
fit3a <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + before_after * treatment, g0 ~ 1, sigma ~ 1)) 
#fit3b <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + before_after * treatment, g0 ~ year + before_after * treatment, sigma ~ year + before_after * treatment)) 
fit4a <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + foxbaited, g0 ~ 1, sigma ~ 1)) 
#fit4b <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + foxbaited, g0 ~ year + foxbaited, sigma ~ year + foxbaited)) 
fit5a <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean)) 
#fit5b <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean)) 
fit6a <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + fox_predicted, g0 ~ 1, sigma ~ 1)) 
#fit6b <- secr.fit(mrch_otways, mask = mask_otways, sessioncov = sesscov, detectfn = 1, trace = FALSE, ncores = 6, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ year + fox_predicted, g0 ~ year + fox_predicted_trapcov, sigma ~ year + fox_predicted_trapcov)) 
otways_fits <- secrlist(fit1a , fit2a, fit2b, fit3a, fit3b, fit4a, fit4b, fit5a, fit5b, fit6a, fit6b)
saveRDS(otways_fits, "models/otways_fits.RData")
AIC(otways_fits, criterion = "AICc")[,-2] # model selection


# END