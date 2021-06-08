# RUN FERAL CAT SECR MODELS - OTWAY RANGES AND GLENELG REGION TOGETHER

rm(list = ls())
options(scipen = 999)

library(secr)  # version 4.3.3
library(dplyr)

# load files
mrch <- readRDS("derived_data/mrch.RData")
masks <- readRDS("derived_data/masks.RData")


# SESSION COVARIATES ------------------------------------------------------
# get an average fox occupancy value for each grid
# glenelg
glenelg_mask_df <- readRDS("derived_data/glenelg_mask_df.RData")
x_g <- glenelg_mask_df %>% 
  group_by(sess) %>% 
  mutate(fox_predicted_mean = mean(fox_predicted))
x_g <- round(unique(x_g$fox_predicted_mean), digits = 2)
# otways
otways_mask_df <-  readRDS("derived_data/otway_mask_df.RData")
x_o <- otways_mask_df %>% 
  group_by(sess) %>% 
  mutate(fox_predicted_mean = mean(fox_predicted))
x_o <- round(unique(x_o$fox_predicted_mean), digits = 2)
# merge together
fox_predicted_mean <- c(x_g, x_o)

# make session cov 
# order = mrch_a, mrch_c, mrch_h, mrch_mc, mrch_s_17, mrch_n_17, mrch_s_18, mrch_n_18, mrch_s_19, mrch_n_19
sesscov <- data.frame(region = factor(c("glenelg", "glenelg", "glenelg", "glenelg", "otways", "otways", "otways", "otways", "otways", "otways")),
                      year = factor(c("2018", "2018", "2018", "2018", "2017", "2017", "2018", "2018", "2019", "2019")),
                      foxbaited = factor(c("unbaited", "baited", "unbaited", "baited", "unbaited", "unbaited", "baited", "unbaited", "baited", "unbaited")),
                      fox_predicted_mean = as.numeric(fox_predicted_mean))
sesscov
str(sesscov)


# FIT MARK-RESIGHT MODELS -------------------------------------------------
# 1) Choose best detector function
fit1_hn <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 0, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, autoini = 9), fixed = list(pID = 1), model = list(D ~ region, g0 ~ region, sigma ~ region))
fit1_hr <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE), fixed = list(pID = 1), model = list(D ~ region, g0 ~ region, sigma ~ region), start = fit1_hn)
fit1_ex <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 2, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE), fixed = list(pID = 1), model = list(D ~ region, g0 ~ region, sigma ~ region), start = fit1_hn)
df_fits <- secrlist(fit1_hn, fit1_ex)
saveRDS(df_fits, "models/df_fits.RData")
df_fits <- readRDS("models/df_fits.RData")
AIC(df_fits, criterion = "AICc")[,-2]
# fit1_ex --> Use this detector function for all subsequent model fitting
# also check buffer size is good
esa.plot(df_fits, max.buffer = 4000, ylim = c(0.001,0.003))
# also go back and double check maks spacing fits with this sigma estimate. 


# 2) Adjust for overdispersion in the unmarked sightings
fit1_adj <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, nsim = 10000), start = fit1_hr, fixed = list(pID = 1), model = list(D ~ region, g0 ~ region, sigma ~ region))
saveRDS(fit1_adj, "models/fit1_adj.RData")
fit1_adj <- readRDS("models/fit1_adj.RData")

fit1_adj$details$chat[1:2]
#Tu       Tm 
#7.003713 1.000000


# 3) fit actual models
# firstly, is there an effect of vegetation type and year?
fit1a  <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ region, g0 ~ region))
fit1b  <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ region + year, g0 ~ region, sigma ~ region))
fit1c  <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ region + vegetation, g0 ~ region, sigma ~ region))
AIC(fit1a, fit1b, fit1c, criterion = "AICc")[,-2] 
## --> no need for vegetation 

# effect of respective fox measure on density 
fit2a <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ foxbaited + region, g0 ~ region, sigma ~ region))
fit3a <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ fox_predicted_mean + region, g0 ~ region, sigma ~ region))
fit4a <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 1, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ fox_predicted + region  g0 ~ region, sigma ~ region))
AIC(fit1b , fit2a, fit3a, fit4a, criterion = "AICc")[,-2]

# effect of respective fox measure on density and detectability 
fit2b <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 2, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ foxbaited + region, g0 ~ foxbaited + region, sigma ~ foxbaited + region))
fit3b <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 2, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ fox_predicted_mean + region, g0 ~ fox_predicted_mean + region, sigma ~ fox_predicted_mean + region))
fit4b <- secr.fit(mrch, mask = masks, sessioncov = sesscov, detectfn = 2, trace = TRUE, ncores = 10, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), method = "Nelder-Mead", start = fit1_adj, model = list(D ~ fox_predicted + region, g0 ~ fox_predicted_trapcov + region, sigma ~ fox_predicted_trapcov + region))
fits <- secrlist(fit1 , fit2a, fit2b, fit3a, fit3b, fit4a, fit4b)
saveRDS(fits, "models/fits.RData")
AIC(fits, criterion = "AICc")[,-2] # model selection


# END