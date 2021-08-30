## Plot estimates from spatial mark-resight models - glenelg and otways

library(secr)  
library(ggplot2)
library(patchwork)
library(dplyr)

# load models
glenelg_fits <- readRDS("models/glenelg_fits.RData")
glenelg_mask_df <- readRDS("derived_data/glenelg_mask_df.RData")
otways_fits <- readRDS("models/otways_fits.RData")
otways_mask_df <- readRDS("derived_data/otway_mask_df.RData")


# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())



#  EFFECT SIZES -----------------------------------------------------

# extract dataframe and change from log to response scale
rr_g <- coef(glenelg_fits$fit_sess2)[c(2,4,6),] %>%
  mutate(beta.resp = exp(beta),
         lcl.resp = exp(lcl),
         ucl.resp = exp(ucl),
         pair = factor(c("Replicate 1", "Replicate 2", "Replicate 3")),
         region = "Glenelg")

rr_o <- coef(otways_fits$fit_sess2)[c(2,4,6),] %>%
  mutate(beta.resp = exp(beta),
         lcl.resp = exp(lcl),
         ucl.resp = exp(ucl),
         year = factor(c("2017", "2018", "2019")),
         region = "Otway")

# mean response ratio
rr <- bind_rows(rr_g, rr_o)
rr
rr2 <- filter(rr, !(region == "Otway" & year == "2017"))
min(rr2$beta.resp)
mean(rr2$beta.resp)
round(max(rr2$beta.resp), digits = 1)


# plot - glenelg
plot_g_difference <- ggplot(rr_g, aes(y = beta.resp, x = pair)) + 
  geom_pointrange(aes(ymin = lcl.resp, ymax = ucl.resp), size = 1, col = "black") +  
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, by = 1)) +
  geom_hline(yintercept = 1, colour = "darkgrey", linetype = "dashed") + 
  labs(title = "", y = "Response ratio", x = "") +
  theme(plot.title = element_text(size=11),           
        axis.title = element_text(size = 10))
plot_g_difference

# plot - otways
plot_o_difference <- ggplot(rr_o, aes(y = beta.resp, x = year)) + 
  geom_pointrange(aes(ymin = lcl.resp, ymax = ucl.resp), size = 1, col = "black") +  
  geom_hline(yintercept = 1, colour = "darkgrey", linetype = "dashed") + 
  labs(title = "", y = "Response Ratio", x = "") +
  theme(plot.title = element_text(size=11),           axis.title = element_text(size = 10))
plot_o_difference



# DENSITY ESTIMATES -------------------------------------------------------

## GLENELG
estimate <- as.data.frame(unlist(sapply(predict(glenelg_fits$fit_sess2), "[", "D","estimate")))
names(estimate)[1] <- "estimate"
estimate <- tibble::rownames_to_column(estimate, "grid")
lower_bound <- as.data.frame(unlist(sapply(predict(glenelg_fits$fit_sess2), "[", "D","lcl")))
names(lower_bound)[1] <- "lcl"
lower_bound <- tibble::rownames_to_column(lower_bound, "grid")
upper_bound <- as.data.frame(unlist(sapply(predict(glenelg_fits$fit_sess2), "[", "D","ucl")))
names(upper_bound)[1] <- "ucl"
upper_bound <- tibble::rownames_to_column(upper_bound, "grid")
temp <- merge(lower_bound, upper_bound, by = "grid")
x <- merge(estimate, temp, by = "grid")
# rename grid
x$grid <- c("LGNP north", "LGNP south", "Annya","Cobboboonee", "Hotspur", "Mt Clay")
# add treatment variable
x$Treatment <- c("Non-impact", "Impact", "Non-impact", "Impact", "Non-impact", "Impact")
# add pair variable
x$Replicate <- c("Replicate 3", "Replicate 3", "Replicate 1", "Replicate 1", "Replicate 2", "Replicate 2")
# covert from hectares to km2
x$estimate <- x$estimate*100
x$lcl <- x$lcl*100
x$ucl <- x$ucl*100
# plot
plot_g_response <- ggplot(x, aes(x = Replicate, y = estimate, color = Treatment)) + 
  geom_point(size = 4, position = position_dodge(width = 0.25)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), position = position_dodge(width = 0.25)) + 
  scale_y_continuous(limits = c(0,0.9), breaks = seq(0, 1, by = 0.1)) +
  labs(title = "", x = "", y = bquote("Cats km"^-2)) +
  scale_color_manual(values=c('red','blue'),  labels = c("Impact", "Non-impact")) +
  theme(plot.title = element_text(size=11),           
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom")
plot_g_response


## OTWAYS
estimate <- as.data.frame(unlist(sapply(predict(otways_fits$fit_sess2), "[", "D","estimate")))
names(estimate)[1] <- "estimate"
estimate <- tibble::rownames_to_column(estimate, "grid")
lower_bound <- as.data.frame(unlist(sapply(predict(otways_fits$fit_sess2), "[", "D","lcl")))
names(lower_bound)[1] <- "lcl"
lower_bound <- tibble::rownames_to_column(lower_bound, "grid")
upper_bound <- as.data.frame(unlist(sapply(predict(otways_fits$fit_sess2), "[", "D","ucl")))
names(upper_bound)[1] <- "ucl"
upper_bound <- tibble::rownames_to_column(upper_bound, "grid")
temp <- merge(lower_bound, upper_bound, by = "grid")
fit_baci_vals <- merge(estimate, temp, by = "grid")
# rename grid
fit_baci_vals$grid <- gsub("session = ", "", fit_baci_vals$grid)
fit_baci_vals$grid <- substr(fit_baci_vals$grid, 0, 9)
new_grid_names <- c(mrch_n_17 = "north 2017", mrch_n_18 = "north 2018", mrch_n_19 = "north 2019", mrch_s_17 = "south 2017", mrch_s_18 = "south 2018", mrch_s_19 = "south 2019")
fit_baci_vals$grid <- as.character(new_grid_names[fit_baci_vals$grid])
# add treatment variable
fit_baci_vals$landscape <- c("Non-impact", "Non-impact", "Non-impact", "Impact", "Impact", "Impact")
fit_baci_vals$year <- c("2017", "2018", "2019", "2017", "2018", "2019")
# covert from hectares to km2
fit_baci_vals$estimate <- fit_baci_vals$estimate*100
fit_baci_vals$lcl <- fit_baci_vals$lcl*100
fit_baci_vals$ucl <- fit_baci_vals$ucl*100
# plot
plot_o_response <- ggplot(fit_baci_vals, aes(x = year, y = estimate, color = landscape, group = landscape)) + 
  geom_point(size = 4, position = position_dodge(width = 0.25)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), position = position_dodge(width = 0.25)) + 
  scale_y_continuous(limits = c(0,1.4), breaks = seq(0, 1.5, by = 0.2)) +
  geom_line(linetype = "dashed", size = 0.5, position = position_dodge(width = 0.25)) + 
  labs(title = "", x = "", y = bquote("Cats km"^-2)) +
  scale_color_manual(values=c('red','blue')) + 
  theme(plot.title = element_text(size=11),           
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom")
plot_o_response


# DENSITY - FOX-CAT CORRELATION -------------------------------------------------------------------
# get linear estimates:
# glenelg:
all_predicted <- predict(glenelg_fits$fit_fox_D, newdata = data.frame(fox_predicted = seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), by = 0.01)))
predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))*100
lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))*100
upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))*100
pr_occ <- round(seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), by = 0.01), digits = 2)
df_glenelg <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
# Otways:
all_predicted <- predict(otways_fits$fit_fox_Ddet, newdata = data.frame(fox_predicted =         seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                        fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                        year = factor(c("2017", "2018", "2019"))))
predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))*100
lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))*100
upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))*100
pr_occ <- round(seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), digits = 2)
df_otways <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
df_otways$year <- rep(c("2017", "2018", "2019"), nrow(df_otways)/3)
df_otways <- filter(df_otways, year == "2019")

# get non-linear estimates:
# glenelg:
all_predicted <- predict(glenelg_fits$fit_nl_fox_D, newdata = data.frame(fox_predicted = seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), by = 0.01)))
predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))*100
lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))*100
upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))*100
pr_occ <- round(seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), by = 0.01), digits = 2)
df_nl_glenelg <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
# Otways:
all_predicted <- predict(otways_fits$fit_nl_fox_Ddet, newdata = data.frame(fox_predicted =         seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                           fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                           year = factor(c("2017", "2018", "2019"))))
predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))*100
lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))*100
upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))*100
pr_occ <- round(seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), digits = 2)
df_nl_otways <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
df_nl_otways$year <- rep(c("2017", "2018", "2019"), nrow(df_nl_otways)/3)
df_nl_otways <- filter(df_nl_otways, year == "2019")


## PLOT
plot_cor <- ggplot(NULL, aes(x = pr_occ, y = predicted_values)) + 
  geom_ribbon(data = df_glenelg, aes(ymin = lower_bound, ymax = upper_bound, fill = "Glenelg"), alpha = 0.2) +
  geom_ribbon(data = df_otways, aes(ymin = lower_bound, ymax = upper_bound, fill = "Otway"), alpha = 0.2) + 
  geom_ribbon(data = df_nl_glenelg, aes(ymin = lower_bound, ymax = upper_bound, fill = "Glenelg"), alpha = 0.2) +
  geom_ribbon(data = df_nl_otways, aes(ymin = lower_bound, ymax = upper_bound, fill = "Otway"), alpha = 0.2) + 
  geom_line(data = df_glenelg, aes(color = "Glenelg"), linetype = 1, size = 1) +
  geom_line(data = df_otways,  aes(color = "Otway"), linetype = 1, size = 1) + 
  geom_line(data = df_nl_glenelg, aes(color = "Glenelg"), linetype = 2, size = 1) +
  geom_line(data = df_nl_otways,  aes(color = "Otway"), linetype = 2, size = 1) + 
  ylim(0,1.55) + 
  scale_fill_manual(name = "Region",
                     breaks = c("Otway", "Glenelg"),
                     values = c("Glenelg" = "#482E1B", "Otway" = "#384566"),
                     guide = "legend") +
  scale_color_manual(name = "Region",
                     breaks = c("Otway", "Glenelg"),
                     values = c("Glenelg" = "#482E1B", "Otway" = "#384566"),
                     guide = "legend") +
  labs(title = "", x = "log(fox occurrence)", y = bquote("Cats km"^-2)) +
  theme(axis.title = element_text(size = 12),       
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))
plot_cor 



# DETECTABILITY - FOX-CAT CORRELATION -------------------------------------------------------------------
# (otways)
## g0 
# linear 
all_predicted <- predict(otways_fits$fit_fox_Ddet, newdata = data.frame(fox_predicted =         seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                        fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                        year = factor(c("2017", "2018", "2019"))))
predicted_values <- unlist(sapply(all_predicted, "[", "g0","estimate"))
lower_bound <- unlist(sapply(all_predicted, "[", "g0","lcl"))
upper_bound <- unlist(sapply(all_predicted, "[", "g0","ucl"))
pr_occ <- seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01)
newdf <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
newdf$year <- rep(c("2017", "2018", "2019"), nrow(newdf)/3)
newdf <- filter(newdf, year == "2019")
# non-linear
all_predicted <- predict(otways_fits$fit_nl_fox_Ddet, newdata = data.frame(fox_predicted =         seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                           fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                           year = factor(c("2017", "2018", "2019"))))
predicted_values <- unlist(sapply(all_predicted, "[", "g0","estimate"))
lower_bound <- unlist(sapply(all_predicted, "[", "g0","lcl"))
upper_bound <- unlist(sapply(all_predicted, "[", "g0","ucl"))
pr_occ <- seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01)
newdf_nl <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
newdf_nl$year <- rep(c("2017", "2018", "2019"), nrow(newdf_nl)/3)
newdf_nl <- filter(newdf_nl, year == "2019")
# plot
plot_g0_fox <- ggplot(NULL, aes(x = pr_occ, y = predicted_values)) + 
  geom_ribbon(data = newdf, aes(ymin = lower_bound, ymax = upper_bound), fill = "#384566", alpha = 0.2) +
  geom_line(data = newdf, size = 1.2, col = "#384566") +
  geom_ribbon(data = newdf_nl, aes(ymin = lower_bound, ymax = upper_bound), fill = "#384566", alpha = 0.2) +
  geom_line(data = newdf_nl, size = 1.2, col = "#384566", linetype = 2)+
  ylim(0, 0.32) + 
  labs(title = "", x = "log(fox occurrence)", y = expression(paste(italic("g")["0"]))) +
  theme(plot.title = element_text(size=11),           axis.title = element_text(size = 10))

## sigma 
# linear
all_predicted <- predict(otways_fits$fit_fox_Ddet, newdata = data.frame(fox_predicted =         seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                        fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                        year = factor(c("2017", "2018", "2019"))))
predicted_values <- unlist(sapply(all_predicted, "[", "sigma","estimate"))
lower_bound <- unlist(sapply(all_predicted, "[", "sigma","lcl"))
upper_bound <- unlist(sapply(all_predicted, "[", "sigma","ucl"))
pr_occ <- seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01)
newdf <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
newdf$year <- rep(c("2017", "2018", "2019"), nrow(newdf)/3)
newdf <- filter(newdf, year == "2019")
# non-linear
all_predicted <- predict(otways_fits$fit_nl_fox_Ddet, newdata = data.frame(fox_predicted =         seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                           fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01), 
                                                                           year = factor(c("2017", "2018", "2019"))))
predicted_values <- unlist(sapply(all_predicted, "[", "sigma","estimate"))
lower_bound <- unlist(sapply(all_predicted, "[", "sigma","lcl"))
upper_bound <- unlist(sapply(all_predicted, "[", "sigma","ucl"))
pr_occ <- seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), by = 0.01)
newdf_nl <- cbind.data.frame(predicted_values, lower_bound, upper_bound, pr_occ)
newdf_nl$year <- rep(c("2017", "2018", "2019"), nrow(newdf_nl)/3)
newdf_nl <- filter(newdf_nl, year == "2019")

# plot
plot_sigma_fox <- ggplot(NULL, aes(x = pr_occ, y = predicted_values)) + 
  geom_ribbon(data = newdf, aes(ymin = lower_bound, ymax = upper_bound), fill = "#384566", alpha = 0.2) +
  geom_line(data = newdf, size = 1.2, col = "#384566") +
  geom_ribbon(data = newdf_nl, aes(ymin = lower_bound, ymax = upper_bound), fill = "#384566", alpha = 0.2) +
  geom_line(data = newdf_nl, size = 1.2, col = "#384566", linetype = 2)+
  ylim(100, 580) + 
  labs(title = "", x = "log(fox occurrence)", y = "sigma") +
  theme(plot.title = element_text(size=11),           axis.title = element_text(size = 10))



# ASSEMBLE PLOTS ----------------------------------------------------------

# 1) Correlation
# 1a) Density (both)
png("C2-manuscript/figs/foxD_600dpi.png", width = 6.5, height = 4, res = 600, units = "in") 
plot_cor + 
  theme(plot.margin = margin(20, 30, 30, 30))
dev.off() 

# 1b (Detectability (otways only)
png("C2-manuscript/figs/foxDet_otways_600dpi.png", width = 8.3, height = 5, res = 600, units = "in")
plot_g0_fox + plot_sigma_fox + plot_annotation(tag_levels = "a") + 
  theme(plot.margin = margin(40, 0, 40, 0))
dev.off() 

# 2) glenelg experimental 
png("C2-manuscript/figs/glenelg_estimates_600dpi.png", width = 8.3, height = 4.5, res = 600, units = "in")
plot_g_response + plot_g_difference + plot_annotation(tag_levels = "a") + 
  theme(plot.margin = margin(0, 0, 0, 0))
dev.off() 

# 3) otways experimental
png("C2-manuscript/figs/otways_estimates_600dpi.png", width = 8.3, height = 4.5, res = 600, units = "in")
plot_o_response + plot_o_difference + plot_annotation(tag_levels = "a") + 
  theme(plot.margin = margin(0, 0, 0, 0))
dev.off() 


# END