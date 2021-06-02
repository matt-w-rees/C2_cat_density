# RUN FERAL CAT SECR MODELS - GLENELG REGION 2018 
# Matthew Rees

rm(list = ls())
options(scipen = 999)

library(secr)  # version 4.3.3
library(mgcv)  
library(ggplot2)
library(patchwork)
library(dplyr)

# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

# load secr data
mask_glenelg <- readRDS("derived_data/mask_glenelg.RData")
mrch_glenelg <- readRDS("derived_data/mrch_glenelg.RData")

# load models 
df_fits <- readRDS("models/glenelg_df_fits.RData")
fit1_adj <- readRDS("models/glenelg_fit1_adj.RData")
glenelg_fits <- readRDS("models/glenelg_fits.RData")

## specify session covariates
# get mean fox occ estimates
glenelg_mask_df <-  readRDS("derived_data/glenelg_mask_df.RData")
x <- glenelg_mask_df %>% 
  group_by(sess) %>% 
  mutate(fox_predicted_mean = mean(fox_predicted))
x <- round(unique(x$fox_predicted_mean), digits = 2)

# make session cov 
sesscov <- data.frame(fox_predicted_mean = x, foxbaited = factor(c("0unbaited", "1baited", "0unbaited", "1baited")), replicate_pair = factor(c("pair1", "pair1", "pair2", "pair2")))
str(sesscov)

# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch_glenelg, mask = mask_glenelg, detectfn = 0, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE), fixed = list(pID = 1))
fit1_ex <- secr.fit(mrch_glenelg, mask = mask_glenelg, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE), fixed = list(pID = 1), start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_ex)
saveRDS(df_fits, "models/glenelg_df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex --> Use this detector function for all subsequent model fitting
# also check buffer size is good
esa.plot(df_fits, max.buffer = 4500, ylim = c(0.001,0.003))
# also go back and double check maks spacing fits with this sigma estimate. 

# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch_glenelg, mask = mask_glenelg, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, nsim = 10000), fixed = list(pID = 1))
saveRDS(fit1_adj, "models/glenelg_fit1_adj.RData")
fit1_adj$details$chat[1:2]
#Tu       Tm 
#2.869429 1.000000 


# 3) fit actual models
fit1  <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ 1, sigma ~ 1))
fit2a <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ 1, sigma ~ 1))
fit2b <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ session, g0 ~ session, sigma ~ session))
fit3a <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ 1, sigma ~ 1))
fit3b <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaited, g0 ~ foxbaited, sigma ~ foxbaited))
fit4a <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ 1, sigma ~ 1))
fit4b <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean))
fit5a <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1, sigma ~ 1))
fit5b <- secr.fit(mrch_glenelg, mask = mask_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 4, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted, sigma ~ fox_predicted))
glenelg_fits <- secrlist(fit1 , fit2a, fit2b, fit3a, fit3b, fit4a, fit4b, fit5a, fit5b)
saveRDS(glenelg_fits, "models/glenelg_fits.RData")
AIC(glenelg_fits, criterion = "AICc")[,-2] # model selection



