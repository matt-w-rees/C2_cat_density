# RUN FERAL CAT SECR MODELS - OTWAY REGION 2018 
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
#5.713383 1.000000 


# 3) Fit actual models
## Null models
fit1      <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ 1, sigma ~ 1))
fit1_Dveg <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ 1, sigma ~ 1))
fit1_g0T  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ 1, g0 ~ T, sigma ~ 1))
fit1_Dveg_g0T  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ vegetation, g0 ~ T, sigma ~ 1))
AIC(fit1, fit1_Dveg, fit1_g0T, fit1_Dveg_g0T, criterion = "AICc")[,-2] 


## 4) Experimental models
# a) same density in the unbaited sites, same effect of baiting in the baited sites
# b) different density in the unbaited sites, same effect of baiting in the baited sites:
# c) different density in the unbaited sites, different effect of baiting in the baited sites (i.e. D ~ session)
fit_BACIa_D     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid + before_after, g0 ~ 1, sigma ~ 1)) 
fit_BACIa_Ddet  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid + before_after, g0 ~ grid + before_after, sigma ~ grid + before_after)) 
fit_BACIb_D     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid * before_after, g0 ~ 1, sigma ~ 1)) 
fit_BACIb_Ddet  <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid * before_after, g0 ~ grid * before_after, sigma ~ grid * before_after)) 
fit_BACIc_D     <- secr.fit(mrch_otways, mask = masks_otways, sessioncov = sesscov, detectfn = 2, trace = FALSE, ncores = 3, details = list(knownmarks = FALSE, chat = fit1_adj$details$chat), fixed = list(pID = 1), start = fit1_adj, model = list(D ~ grid + year, g0 ~ 1, sigma ~ 1)) 
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



# PLOTS -------------------------------------------------------------------
# load models
#otways_fits <- readRDS("models/secr/otways/otways_fits.RData")
#fit_baci <- readRDS("models/secr/otways/otways_fit_baci.RData")
otways_mask_df <- readRDS("derived_data/otway_mask_df.RData")

library(ggplot2)
# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

# 1) EXPERIMENTAL 
# A) EFFECT SIZE

# B) DENSITY ESTIMATES
## Extract estimates
estimate <- as.data.frame(unlist(sapply(predict(fit_BACId_Ddet), "[", "D","estimate")))
names(estimate)[1] <- "estimate"
estimate <- tibble::rownames_to_column(estimate, "grid")
lower_bound <- as.data.frame(unlist(sapply(predict(fit_BACId_Ddet), "[", "D","lcl")))
names(lower_bound)[1] <- "lcl"
lower_bound <- tibble::rownames_to_column(lower_bound, "grid")
upper_bound <- as.data.frame(unlist(sapply(predict(fit_BACId_Ddet), "[", "D","ucl")))
names(upper_bound)[1] <- "ucl"
upper_bound <- tibble::rownames_to_column(upper_bound, "grid")
temp <- merge(lower_bound, upper_bound, by = "grid")
fit_baci_vals <- merge(estimate, temp, by = "grid")
# rename grid
fit_baci_vals$grid <- gsub("session = ", "", fit_baci_vals$grid)
new_grid_names <- c(mrch_n_17 = "north 2017", mrch_n_18 = "north 2018", mrch_n_19 = "north 2019", mrch_s_17 = "south 2017", mrch_s_18 = "south 2018", mrch_s_19 = "south 2019")
fit_baci_vals$grid <- as.character(new_grid_names[fit_baci_vals$grid])
# add treatment variable
fit_baci_vals$Treatment <- c("Non-treatment", "Non-treatment", "Non-treatment", "Treatment", "Treatment", "Treatment")
fit_baci_vals$Year <- c("2017", "2018", "2019", "2017", "2018", "2019")
# covert from hectares to km2
fit_baci_vals$estimate <- fit_baci_vals$estimate*100
fit_baci_vals$lcl <- fit_baci_vals$lcl*100
fit_baci_vals$ucl <- fit_baci_vals$ucl*100
## Plot
plot_o_cat <- ggplot(fit_baci_vals, aes(x = Year, y = estimate, color = Treatment, group = Treatment)) + 
  geom_point(size = 4, position = position_dodge(width = 0.25)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), position = position_dodge(width = 0.25)) + 
  ylim(0,2.1) + 
  geom_line(linetype = "dashed", size = 0.5, position = position_dodge(width = 0.25)) + 
  labs(title = "", x = "Year", y = bquote("Cats per km"^2)) +
  scale_color_manual(values=c('blue','red')) + 
  theme(plot.background = element_rect(fill = 'white'),
        legend.position = "none",
        #legend.title = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
plot_o_cat


# 2) CORRELATION 
## D ~ FOX
# extract values
all_predicted <- predict(otways_fits$fit_fox_Dd, 
                         newdata = data.frame(fox_predicted = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01),
                                              fox_predicted_trapcov = mean(otways_mask_df$fox_predicted)))
predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))*100
lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))*100
upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))*100
pr_occ <- seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01)
newdf <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
# plot
plot_cat_fox <- ggplot(newdf, aes(x = pr_occ, y = predicted_values)) + 
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "grey85") +
  geom_line(size = 1)+
  ylim(0,1.5) + 
  labs(title = "", x = "Fox Pr(occupancy)", y = bquote("Cats per km"^2)) +
  scale_color_manual(values=c('blue','red')) + 
  theme(plot.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
plot_cat_fox



## Assemble plots and save
png("figs/fig5_otway_600dpi.png", width = 9, height = 8, res = 600, units = "in")
(plot_o_cat + plot_effect) / (plot_cat_fox + plot_cat_prey) + plot_annotation(tag_levels = 'A', title = "Otway region")
dev.off() 


# END 