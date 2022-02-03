library(dplyr)
library(mgcv)
library(secr)
library(ggplot2)


# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())


# FOX RECORDS ----------------------------------------------------------------
# load fox presence absence records
records <- read.csv("raw_data/spp_records_pa.csv")
# lower glenelg:
records_lg <- read.csv("raw_data/lower_glenelg_fox_pa.csv")
records_lg$station_year <- records_lg$station
# merge together
records <- bind_rows(records, records_lg)


# add forest block (+ year for otways covariates)
records$block <- if_else(grepl("^A", records$station), "annya", "test")
records$block <- if_else(grepl("^C", records$station), "cobb", records$block)
records$block <- if_else(grepl("^H", records$station), "hotspur", records$block)
records$block <- if_else(grepl("^M", records$station), "mt_clay", records$block)
records$block <- if_else(grepl("^U", records$station), "o_ub", records$block)
records$block <- if_else(grepl("^T", records$station), "o_t", records$block)
records$block <- if_else(records$region == "otways", paste0(records$block, "_", records$year), records$block)
records$block <- if_else(grepl("^L", records$station) & as.integer(substr(records$station, 5, 8))  >= 65, "lg_n", records$block)
records$block <- if_else(grepl("^L", records$station) & as.integer(substr(records$station, 5, 8))  < 65, "lg_s", records$block)
unique(records$block)
unique(is.na(records$block))


## estimate occupancy per landscape

# specify records class first
records <- mutate(records, fox = as.integer(fox),
                  year = factor(year, ordered = FALSE),
                  region = factor(region, ordered = FALSE),
                  station = factor(station, ordered = FALSE),
                  block = factor(block, ordered = FALSE),
                  survey_duration = as.integer(survey_duration))


# fit GAM - glenelg 
gam_fox <- bam(fox ~ block + 
                     offset(log(survey_duration)), 
                 data = records, family = binomial, discrete = TRUE) 
plot(gam_fox, pages = 1, scheme = 2)
summary(gam_fox)

# predict estimates to new dataset
df <-  expand.grid(block = unique(records$block), survey_duration = 60)
estimates <- cbind(df, predict.gam(gam_fox, newdata = df, type = "link",  se.fit = TRUE))
estimates <- estimates %>% 
  mutate(fox_occ_lcl = fit - (1.96 * se.fit),
         fox_occ_ucl = fit + (1.96 * se.fit)) %>% 
  rename(fox_occ_pred = fit, fox_occ_se = se.fit)


# EXTRACT CAT DENSITY ESTIMATES -------------------------------------------

# get cat density estimates
glenelg_fits <- readRDS("models/glenelg_fits.RData")
otways_fits <- readRDS("models/otways_fits.RData")

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
# covert from hectares to km2
x$estimate <- x$estimate*100
x$lcl <- x$lcl*100
x$ucl <- x$ucl*100



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
# covert from hectares to km2
fit_baci_vals$estimate <- fit_baci_vals$estimate*100
fit_baci_vals$lcl <- fit_baci_vals$lcl*100
fit_baci_vals$ucl <- fit_baci_vals$ucl*100


# merge glenelg and otway estimates
cat <- bind_rows(x, fit_baci_vals)


# MERGE FOX CAT DATAFRAMES -----------------------------------------------------------
# match up grid to block
estimates$grid <- c("Annya", "Cobboboonee", "Hotspur", "Mt Clay", "north 2017", "north 2018", 
                    "north 2019", "south 2017", "south 2018", 
                    "south 2019", "LGNP south", "LGNP north")
estimates$survey_duration <- NULL
estimates$block <- NULL
estimates_df <- left_join(cat, estimates, by = "grid")
estimates_df$Region <- c(rep("Glenelg", 6), rep("Otway", 6))

estimates_df$baiting <- c("U", "B", "U", "B", "U", "B", "U", "U", "U", "U", "B", "B")


# PLOT --------------------------------------------------------------------

plot_cat_fox <- ggplot(estimates_df, aes(y = estimate, x = fox_occ_pred,  colour = grid, group = Region)) + 
  geom_pointrange(aes(ymin = lcl, ymax = ucl)) + 
  geom_pointrange(aes(xmin = fox_occ_lcl, xmax = fox_occ_ucl)) + 
  geom_point(size = 4) +
  geom_smooth(method = "lm", color = "black") + 
  geom_text(aes(label=baiting),hjust=1.5, vjust=1.5) +
  labs(title = "", x = "log(fox occurrence)", y = bquote("Cats km"^-2)) +
  theme(plot.title = element_text(size=11),           
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom")
plot_cat_fox
