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

# get mean fox_prediceted value per grid to add a session covariate
glenelg_mask_df <- readRDS("derived_data/glenelg_mask_df.RData")
x <- glenelg_mask_df %>% 
  group_by(sess) %>% 
  mutate(fox_predicted_mean = mean(fox_predicted)) %>%
  select(sess, fox_predicted_mean)
x <- round(unique(x)$fox_predicted_mean, digits = 2)
x

# specify session covariates
sesscov <- data.frame(fox_predicted_mean = x,
                      pair = factor(c("p1", "p1", "p2", "p2")),
                      foxbaiting = factor(c("0unbaited", "1baited", "0unbaited", "1baited")))
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
fits1_a <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ 1, sigma ~ 1))
fits1_b <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ 1, sigma ~ 1))
fits1_c <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ T, sigma ~ 1))
fits1_d <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ T, sigma ~ 1))
AIC(fits1_a, fits1_b, fits1_c, fits1_d, criterion = "AICc")[,-2] 


## 4) Experimental models
# a) same density in the unbaited sites, same effect of baiting in the baited sites
# b) different density in the unbaited sites, same effect of baiting in the baited sites:
# c) different density in the unbaited sites, different effect of baiting in the baited sites (i.e. D ~ session)
fits2_a <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ foxbaiting, sigma ~ foxbaiting)) 
fits2_b <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaiting, g0 ~ 1, sigma ~ 1)) 
fits2_c <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ foxbaiting, g0 ~ foxbaiting, sigma ~ foxbaiting)) 
fits2_d <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair + foxbaiting, g0 ~ 1, sigma ~ 1)) 
fits2_e <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair * foxbaiting, g0 ~ 1, sigma ~ 1)) 
fits2_f <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair + foxbaiting, g0 ~ foxbaiting, sigma ~ foxbaiting)) 
fits2_g <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ pair * foxbaiting, g0 ~ foxbaiting, sigma ~ foxbaiting)) 
AIC(fits1_a, fits2_a, fits2_b, fits2_c, fits2_d, fits2_e, fits2_f, fits2_g, criterion = "AICc")[,-2] 


## 5) Correlative models with fox occupancy
fits3_a <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean))
fits3_b <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ 1, sigma ~ 1)) 
fits3_c <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted_mean, g0 ~ fox_predicted_mean, sigma ~ fox_predicted_mean))
fits3_d <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov))
fits3_e <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ 1, sigma ~ 1)) 
fits3_f <- secr.fit(mrch_glenelg, mask = masks_glenelg, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ fox_predicted, g0 ~ fox_predicted_trapcov, sigma ~ fox_predicted_trapcov)) 
AIC(fits1_a, fits3_a, fits3_b, fits3_c, fits3_d, fits3_e, fits3_f, criterion = "AICc")[,-2] 


## 6) Combine models, save and compare AICc
glenelg_fits <- secrlist(fits1_a, fits1_b, fits1_c, fits1_d, fits2_a, fits2_b, fits2_c, fits2_d, fits2_e, fits2_f, fits2_g, fits3_a, fits3_b, fits3_c, fits3_d, fits3_e, fits3_f)
saveRDS(glenelg_fits,   "models/glenelg_fits.RData")
#glenelg_fits <- readRDS("models/glenelg_fits.RData")
AIC(glenelg_fits, criterion = "AICc")[,-2] 



# PLOTS -------------------------------------------------------------------
# load models
#glenelg_fits <- readRDS("models/secr/glenelg/glenelg_fits.RData")
#fit_ci2 <- readRDS("models/secr/glenelg/glenelg_fit_ci2.RData")
library(ggplot2)
# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

# 1) EXPERIMENTAL 
# B) DENSITY ESTIMATES
## Extract estimates
# pair 1: 
estimate <- as.data.frame(unlist(sapply(predict(fit_CIc_D), "[", "D","estimate")))
names(estimate)[1] <- "estimate"
estimate <- tibble::rownames_to_column(estimate, "grid")
lower_bound <- as.data.frame(unlist(sapply(predict(fit_CIc_D), "[", "D","lcl")))
names(lower_bound)[1] <- "lcl"
lower_bound <- tibble::rownames_to_column(lower_bound, "grid")
upper_bound <- as.data.frame(unlist(sapply(predict(fit_CIc_D), "[", "D","ucl")))
names(upper_bound)[1] <- "ucl"
upper_bound <- tibble::rownames_to_column(upper_bound, "grid")
temp <- merge(lower_bound, upper_bound, by = "grid")
x <- merge(estimate, temp, by = "grid")

# rename grid
x$grid <- c("Annya","Cobboboonee", "Hotspur", "Mt Clay")
# add treatment variable
x$Treatment <- c("Non-treatment", "Treatment", "Non-treatment", "Treatment")
# add pair variable
x$Replicate <- c("Replicate 1", "Replicate 1", "Replicate 2", "Replicate 2")
# covert from hectares to km2
x$estimate <- x$estimate*100
x$lcl <- x$lcl*100
x$ucl <- x$ucl*100
# plot
plot_g_cat <- ggplot(x, aes(x = Replicate, y = estimate, color = Treatment)) + 
  geom_point(size = 4, position = position_dodge(width = 0.25)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), position = position_dodge(width = 0.25)) + 
  ylim(0,1) + labs(title = "", x = "", y = bquote("Cats per km"^2)) +
  scale_color_manual(values=c('blue','red'),  labels = c("Non-treatment", "Treatment")) + 
  theme(plot.background = element_rect(fill = 'white'),
        legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
plot_g_cat



# 2) CORRELATION 
## D ~ FOX
# extract values
all_predicted <- predict(fit_fox_D, newdata = data.frame(fox_predicted = seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), by = 0.01)))
predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))*100
lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))*100
upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))*100
pr_occ <- seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), by = 0.01)
newdf <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
# plot
plot_cat_fox <- ggplot(newdf, aes(x = pr_occ, y = predicted_values)) + 
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "grey85") +
  geom_line(size = 1)+
  ylim(0,1) + 
  labs(title = "", x = "Fox Pr(occupancy)", y = bquote("Cats per km"^2)) +
  scale_color_manual(values=c('blue','red')) + 
  theme(plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
plot_cat_fox

## Assemble plots and save
png("figs/fig4_glenelg_600dpi.png", width = 13, height = 4, res = 600, units = "in")
(plot_g_cat + plot_effect + plot_cat_fox) + plot_annotation(tag_levels = 'A', title = "Glenelg region") 
dev.off() 

# END 

# END