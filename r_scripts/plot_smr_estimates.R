library(secr)  
library(ggplot2)
library(patchwork)
library(dplyr)


# load models
glenelg_fits <- readRDS("models/glenelg_fits.RData")
glenelg_mask_df <- readRDS("derived_data/glenelg_mask_df.RData")
otways_fits <- readRDS("models/otways_fits.RData")
otways_mask_df <- readRDS("derived_data/otway_mask_df.RData")


# Customise ggplot theme -----------------------------------------------------
theme_matty <- function () { 
  theme_bw(base_size=14) %+replace% 
    theme(
           # background / panels
           panel.background = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_rect(colour = "lightgrey", fill=NA, size=1),
           # font sizes
           axis.text = element_text(size = 14),
           axis.title = element_text(size = 14),
           title = element_text(size = 12),
           # margins
           axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
           axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
           # legend
           legend.title = element_blank(),
           legend.text = element_text(size = 12),
           legend.position = "bottom"
  )
}



#  EFFECT SIZE -----------------------------------------------------
## GLENELG
# difference between NI and I landscapes: 
x <- coef(glenelg_fits$fit_sess1)[c(2,4,6),c(1,3:4)]
x$pair <- factor(c("Replicate 1", "Replicate 2", "Replicate 3"))
# plot
plot_g_difference <- ggplot(x, aes(x = beta, y = reorder(pair, desc(pair)))) + 
  geom_pointrange(aes(xmin = lcl, xmax = ucl), size = 1, col = "black") +  
  xlim(-2,2.2) + 
  geom_vline(xintercept = 0, colour = "darkgrey", linetype = "dashed") + 
  labs(title = "", x = "Difference between impact and non-impact landscapes", y = "") +
  theme_matty()
plot_g_difference

# density estimates 
# pair 1: 
estimate <- as.data.frame(unlist(sapply(predict(glenelg_fits$fit_sess1), "[", "D","estimate")))
names(estimate)[1] <- "estimate"
estimate <- tibble::rownames_to_column(estimate, "grid")
lower_bound <- as.data.frame(unlist(sapply(predict(glenelg_fits$fit_sess1), "[", "D","lcl")))
names(lower_bound)[1] <- "lcl"
lower_bound <- tibble::rownames_to_column(lower_bound, "grid")
upper_bound <- as.data.frame(unlist(sapply(predict(glenelg_fits$fit_sess1), "[", "D","ucl")))
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
  ylim(0,0.8) + 
  labs(title = "", x = "", y = bquote("Cats per km"^2)) +
  scale_color_manual(values=c('red','blue'),  labels = c("Impact", "Non-impact")) +
  theme_matty()
plot_g_response




## OTWAY
# extract coefficients 
x <- coef(otways_fits$fit_sess4)[c(2, 4, 6),c(1,3:4)]
x$year <- factor(c("2017", "2018", "2019"))
plot_o_difference <- ggplot(x, aes(x = beta, y = reorder(year, desc(year)))) + 
  geom_pointrange(aes(xmin = lcl, xmax = ucl), size = 1, col = "black") +  
  xlim(-1.2,1.2) + 
  labs(title = "", x = "Difference between impact and non-impact landscapes", y = "") +
  theme_matty()
plot_o_difference

## Extract estimates
estimate <- as.data.frame(unlist(sapply(predict(otways_fits$fit_sess4), "[", "D","estimate")))
names(estimate)[1] <- "estimate"
estimate <- tibble::rownames_to_column(estimate, "grid")
lower_bound <- as.data.frame(unlist(sapply(predict(otways_fits$fit_sess4), "[", "D","lcl")))
names(lower_bound)[1] <- "lcl"
lower_bound <- tibble::rownames_to_column(lower_bound, "grid")
upper_bound <- as.data.frame(unlist(sapply(predict(otways_fits$fit_sess4), "[", "D","ucl")))
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
## Plot
plot_o_response <- ggplot(fit_baci_vals, aes(x = year, y = estimate, color = landscape, group = landscape)) + 
  geom_point(size = 4, position = position_dodge(width = 0.25)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl), position = position_dodge(width = 0.25)) + 
  ylim(0,1.8) + 
  geom_line(linetype = "dashed", size = 0.5, position = position_dodge(width = 0.25)) + 
  labs(title = "", x = "Year", y = bquote("Cats per km"^2)) +
  scale_color_manual(values=c('blue','red')) + 
  theme_matty()
plot_o_response




# FOX CORRELATION PLOT -------------------------------------------------------------------
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
  geom_line(data = df_glenelg, aes(color = "Glenelg"), linetype = 1, size = 1.2) +
  geom_line(data = df_otways,  aes(color = "Otway"), linetype = 1, size = 1.2) + 
  geom_line(data = df_nl_glenelg, aes(color = "Glenelg"), linetype = 2, size = 1.2) +
  geom_line(data = df_nl_otways,  aes(color = "Otway"), linetype = 2, size = 1.2) + 
  ylim(0,1.55) + 
  scale_fill_manual(name = "Region",
                     breaks = c("Otway", "Glenelg"),
                     values = c("Glenelg" = "#482E1B", "Otway" = "#384566"),
                     guide = "legend") +
  scale_color_manual(name = "Region",
                     breaks = c("Otway", "Glenelg"),
                     values = c("Glenelg" = "#482E1B", "Otway" = "#384566"),
                     guide = "legend") +
  labs(title = "", x = "log(fox occurrence)", y = bquote("Cats per km"^2)) +
  theme_matty()
plot_cor 



# COR - DETECTABILITY -----------------------------------------------------
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
  theme_matty()

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
  theme_matty()



# ASSEMBLE PLOTS ----------------------------------------------------------

# 1) Correlation
# 1a) Density (both)
png("C2-manuscript/figs/foxD_600dpi.png", width = 7, height = 6, res = 600, units = "in")
plot_cor
dev.off() 

# 1b (Detectability (otways only)
png("C2-manuscript/figs/foxDet_otways_600dpi.png", width = 12, height = 5, res = 600, units = "in")
plot_g0_fox + plot_sigma_fox + plot_annotation(tag_levels = "a")
dev.off() 


# 2) glenelg experimental 
png("C2-manuscript/figs/glenelg_estimates_600dpi.png", width = 13, height = 6, res = 600, units = "in")
plot_g_difference + plot_g_response + plot_annotation(tag_levels = "a")
dev.off() 

# 3) otways experimental
png("C2-manuscript/figs/otways_estimates_600dpi.png", width = 13, height = 6, res = 600, units = "in")
plot_o_difference + plot_o_response +  plot_annotation(tag_levels = "a")
dev.off() 


# END