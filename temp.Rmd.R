
```{r, load data, include = FALSE}

# load camera-trap station info
records <- read_csv("/Users/mrees2/Dropbox/personal/matt/github/compile-camera-records-matt/derived_data/records_matt_clean.csv")
camdata <- distinct(records, station_year, .keep_all = TRUE)
camdata_otways <- camdata[which(camdata$region == "otways"),] 
camdata_glenelg <- camdata[which(camdata$region == "glenelg"),]

# load GAM data  
gam_data <- read_csv("derived_data/gam_data.csv") 
# add treatment variable
gam_data$treatment <- if_else(gam_data$grid == "cobb" | gam_data$grid == "mtclay" | gam_data$grid == "south", "treatment", "non-treatment")
# split for each region
gam_data_otways <- gam_data[which(gam_data$region == "otways"),]  
gam_data_glenelg <- gam_data[which(gam_data$region == "glenelg"),]  
# add grid_year col for otways
gam_data_otways$grid_year <- paste0(gam_data_otways$grid, "_", gam_data_otways$year)

# load masks as dataframes (for covariates)
otways_mask_df <- readRDS("derived_data/otway_mask_df.RData")
glenelg_mask_df <- readRDS("derived_data/glenelg_mask_df.RData")

# load SECR capture histories - each session (but seperate unID from unmarked)
# mt clay (no black cats)
ch_mc <- readRDS("derived_data/ch_mc.RData")
unid_mc <- as.matrix(read.csv("derived_data/unid_mc.csv")[, -1]) 
mrch_mc <- addSightings(ch_mc, uncertain = unid_mc)
# hotspur 
ch_h <- readRDS("derived_data/ch_h.RData")
unid_h <- as.matrix(read.csv("derived_data/unid_h.csv")[, -1]) 
unm_h <- as.matrix(read.csv("derived_data/unm_h.csv")[, -1]) 
mrch_h <- addSightings(ch_h, unmarked = unm_h, uncertain = unid_h)
# annya 
ch_a <- readRDS("derived_data/ch_a.RData")
unid_a <- as.matrix(read.csv("derived_data/unid_a.csv")[, -1]) 
unm_a <- as.matrix(read.csv("derived_data/unm_a.csv")[, -1]) 
mrch_a <- addSightings(ch_a, unmarked = unm_a, uncertain = unid_a)
# cobbob 
ch_c <- readRDS("derived_data/ch_c.RData")
unid_c <- as.matrix(read.csv("derived_data/unid_c.csv")[, -1]) 
unm_c <- as.matrix(read.csv("derived_data/unm_c.csv")[, -1]) 
mrch_c <- addSightings(ch_c, unmarked = unm_c, uncertain = unid_c)
# northern grid 2017 
ch_n_17 <- readRDS("derived_data/ch_n_17.RData")
unid_n_17 <- as.matrix(read.csv("derived_data/unid_n_17.csv")[, -1]) 
unm_n_17 <- as.matrix(read.csv("derived_data/unm_n_17.csv")[, -1]) 
mrch_n_17 <- addSightings(ch_n_17, unmarked = unm_n_17, uncertain = unid_n_17)
# northern grid 2018 
ch_n_18 <- readRDS("derived_data/ch_n_18.RData")
unid_n_18 <- as.matrix(read.csv("derived_data/unid_n_18.csv")[, -1]) 
unm_n_18 <- as.matrix(read.csv("derived_data/unm_n_18.csv")[, -1]) 
mrch_n_18 <- addSightings(ch_n_18, unmarked = unm_n_18, uncertain = unid_n_18)
# northern grid 2019 
ch_n_19 <- readRDS("derived_data/ch_n_19.RData")
unid_n_19 <- as.matrix(read.csv("derived_data/unid_n_19.csv")[, -1]) 
unm_n_19 <- as.matrix(read.csv("derived_data/unm_n_19.csv")[, -1]) 
mrch_n_19 <- addSightings(ch_n_19, unmarked = unm_n_19, uncertain = unid_n_19)
# southern grid 2017 
ch_s_17 <- readRDS("derived_data/ch_s_17.RData")
unid_s_17 <- as.matrix(read.csv("derived_data/unid_s_17.csv")[, -1]) 
unm_s_17 <- as.matrix(read.csv("derived_data/unm_s_17.csv")[, -1]) 
mrch_s_17 <- addSightings(ch_s_17, unmarked = unm_s_17, uncertain = unid_s_17)
# southern grid 2018 
ch_s_18 <- readRDS("derived_data/ch_s_18.RData")
unid_s_18 <- as.matrix(read.csv("derived_data/unid_s_18.csv")[, -1]) 
unm_s_18 <- as.matrix(read.csv("derived_data/unm_s_18.csv")  [, -1]) 
mrch_s_18 <- addSightings(ch_s_18, unmarked = unm_s_18, uncertain = unid_s_18)
# southern grid 2019 
ch_s_19 <- readRDS("derived_data/ch_s_19.RData")
unid_s_19 <- as.matrix(read.csv("derived_data/unid_s_19.csv")[, -1]) 
unm_s_19 <- as.matrix(read.csv("derived_data/unm_s_19.csv")[, -1]) 
mrch_s_19 <- addSightings(ch_s_19, unmarked = unm_s_19, uncertain = unid_s_19)

# combine capthists per region
mrch_glenelg <- MS.capthist(mrch_a, mrch_c, mrch_h, mrch_mc)
mrch_otways <- MS.capthist(mrch_s_17, mrch_n_17, mrch_s_18, mrch_n_18, mrch_s_19, mrch_n_19)

# load SECR models
df_fits_glenelg <- readRDS("models/secr/glenelg/glenelg_df_fits.RData")
glenelg_fits <- readRDS("models/secr/glenelg/glenelg_fits.RData")
df_fits_otways <- readRDS("models/secr/otways/otways_df_fits.RData")
otways_fits <- readRDS("models/secr/otways/otways_fits.RData")

# load GAMs
gam_g_fox <- readRDS("models/gams/glenelg/gam_g_fox.RData")
gam_g_btp <- readRDS("models/gams/glenelg/gam_g_btp.RData")
gam_g_pot <- readRDS("models/gams/glenelg/gam_g_pot.RData")
gam_g_rtp <- readRDS("models/gams/glenelg/gam_g_rtp.RData")
gam_g_sbb <- readRDS("models/gams/glenelg/gam_g_sbb.RData")
gam_g_sm  <- readRDS("models/gams/glenelg/gam_g_sm.RData")
gam_o_fox <- readRDS("models/gams/otways/gam_o_fox.RData")
gam_o_btp <- readRDS("models/gams/otways/gam_o_btp.RData")
gam_o_pot <- readRDS("models/gams/otways/gam_o_pot.RData")
gam_o_rtp <- readRDS("models/gams/otways/gam_o_rtp.RData")
gam_o_band<- readRDS("models/gams/otways/gam_o_band.RData")
gam_o_sm  <- readRDS("models/gams/otways/gam_o_sm.RData")

```


\newpage


# Survey methods

## Timeline

In the Glenelg region, we consecutively surveyed two pairs of camera-trap grids during the first half of 2018 (Figure  \@ref(fig:camop)A). First, Annya and Cobboboonee; second, Hotspur and Mt Clay. By this time, continuous fox control had been occurring for approximately 13 years.

In the Otway Ranges, fox control commenced in November 2017, but had an approximate five month pause between July and December 2018. We surveyed one pair of camera-trap grids three times annually in June – September, from 2017 - 2019 (Figure  \@ref(fig:camop)B).  

$~$
  
  $~$
  
  ```{r, camop, echo = FALSE, fig.height = 6, fig.width=4.5, fig.cap = "Camera-trap operation times in the Glenelg region (A) and Otway Ranges (B). Each blue, horizontal line represents one camera-trap. Grey shading indicates periods of fox control in the Otway Ranges."}

par(mfrow = c(2,1), mar = c(2.6, 2.1, 2.1, 2.1))

## GLENELG
camop_glenelg <- cameraOperation(CTtable = camdata_glenelg, stationCol = "station", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
# plot
camOp <- camop_glenelg
which.tmp <- grep(as.Date(colnames(camOp)), pattern = "01$")
label.tmp <- format(as.Date(colnames(camOp))[which.tmp], "%Y-%m")
at.tmp <- which.tmp / ncol(camOp)
values_tmp <- sort(na.omit(unique(c(camOp))))
image(t(as.matrix(camOp)), xaxt = "n", yaxt = "n", col = "#1b98e0")
axis(1, at = at.tmp, labels = label.tmp)
# axis(2, at = seq(from = 0, to = 1, length.out = nrow(camOp)), labels = NA, las = 1)
abline(v = at.tmp, col = rgb(0,0,0, 0.1))
mtext("A", side = 3, line = 1, at = 0)

## OTWAYS
camdata_otways <- arrange(camdata_otways, date_start)
camop_otways    <- cameraOperation(CTtable = camdata_otways,  stationCol = "station_year", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
camOp <- camop_otways
which.tmp <- grep(as.Date(colnames(camOp)), pattern = "01$")
label.tmp <- format(as.Date(colnames(camOp))[which.tmp], "%Y-%m")
at.tmp <- which.tmp / ncol(camOp)
values_tmp <- sort(na.omit(unique(c(camOp))))
image(t(as.matrix(camOp)), xaxt = "n", yaxt = "n", col = "#1b98e0")
axis(1, at = at.tmp, labels = label.tmp)
# axis(2, at = seq(from = 0, to = 1, length.out = nrow(camOp)), labels = NA, las = 1)
abline(v = at.tmp, col = rgb(0,0,0, 0.1))
# plot baiting schedule
polygon(x = c(0.18, 0.458, 0.458, 0.18), y = c(0, 0, 2, 2), col = "grey80")  
polygon(x = c(0.645, 1, 1, 0.645), y = c(0, 0, 2, 2), col = "grey80")  
#replot cams
image(t(as.matrix(camOp)), xaxt = "n", yaxt = "n", col = "#1b98e0", add = TRUE)
abline(v = at.tmp, col = rgb(0,0,0, 0.1))
mtext("B", side = 3, line = 1, at = 0)

```

\newpage



## Camera-trap deployment
In the Glenelg region, we deployed camera sites once, whereas in the Otways, we redeployed cameras three times annually. All 2017 camera-sites were resurveyed each year, except for four logistically challenging sites in the southern grid. In 2018, we added 16 additional sites in the southern grid, as well as 36 additional sites in the northern grid. These additional sites were resurveyed in 2019. Camera-trap grids ranged from 3442 - 7819 camera-trap nights across both regions (Table 1).

At each site, we deployed a singular remote trail camera with infrared flash and temperature-in-motion detector. The vast majority of cameras were were Reconyx Hyperfire HC600, but a small proportion were PC900’s HF2X’s, infrared camera-traps (Reconyx, Holmen, Wisconsin). We programmed camera’s to their highest sensitivity and to take five consecutive photographs when triggered (no quiet period). We attached each camera to a tree, approximately 30 cm above the ground, and facing toward a lure 2 - 2.5 metres away. The lure comprised an oil-absorbing cloth doused in tuna oil and placed inside a PVC pipe container with a mesh top. We secured each lure to the top of a 1 metre wooden stake and attached a handful of small white feathers to the outside of the PVC pipe container. We cleared vegetation in the camera’s line-of-sight to reduce false triggers and avoid obscuring cat coat markings in images. While these cameras were aimed to target predators, they were also effective at detecting mammalian prey species.

$~$
  
  $~$
  
  
  ```{r cam, echo = FALSE, out.width="100%", fig.cap = "Example of a camera-trap set-up"}

knitr::include_graphics("figs/camtrap1.jpg")

```



\newpage


# Individual cat identification

We first labelled every camera-trap image with a species metadata tag using DigiKam software [www.digikam.org](https://www.digikam.org){target="_blank"}. We also tagged cats based on their coat type: black, spotty tabby, swirly tabby, ginger and other: coats with multiple colour blends (Fig. 3). This allowed us to summarise species detections and extract cat images using the "camtrapR" R package (Jürgen, Sollman, Courtiol & Wilting 2016; *Methods in Ecology and Evolution*).

We considered all black cats to be of the 'unmarked' category in spatial mark-resight models - even the few with white splotches on their underside (these couldn't always be seen as cats move with their head down).

In the remaining coat categories, we identified individual cats based on their unique coat markings. The ability to identify individuals substantially increased as the image library for each cat increased. Therefore we made the easiest identifications first to build up these libraries, before making decisions on the less obvious detections. We examined and matched all coats markings seen between two particular defections. Markings on the front legs were most often seen and particularly useful for ID's as the patterns do not skew as much with different body positions. On the whole, unidentifiable detections were mainly due to only part of a cat appearing in the frame, or because photos were blurry (fast moving cats or a fogged up camera lens).
                                                                                                                                                     
                                                                                                                                                     However, we were left with a small number of instances (less than ten) where only left or right flanks could be seen. In this case, the side with the most repeat detections was labelled as an individual, whereas the side with the least number of detections was considered unidentifiable. 
                                                                                                                                                     
                                                                                                                                                     Additionally, an extremely small portion of cats in the Otways had ginger coats. When ginger coats are photographed with an infrared flash, they become overexposed and appear consistently white (bottom-right corner in Fig. 3). We only had one detection of a ginger cat without an infrared coat. Therefore, if there were multiple ginger cat detections in a single grid, we treated them in same way as one-sided flank detections.
                                                                                                                                                     
                                                                                                                                                     One observer identified the feral cats in the Glenelg region (MR). In the 2017 and 2018 Otway datasets, where there were substantially more cat detections and fewer distinct coat patterns, two independent observers identified individual cats and discrepancies between observers were reviewed together until consensus was reached (MR, MLP, BH). If no consensus could be reached, the cat was considered unidentifiable. In the 2019 Otway dataset, many of the identified cats were sighted in the previous surveys – these larger individual libraries meant that cats could be identified more easily so only one observer was necessary (MR). We also made use of additional cat images taken within the Otway grids (just before each of our surveys) by white flash camera-traps from a complementary study (Z. Banikos, unpublished data). This provided additional and higher quality images (due to the white flash) of individuals.
                                                                                                                                                     
                                                                                                                                                     We were therefore left with three groups of cats: unmarked (black cats), marked (cats which could be identified to the individual-level with complete certainty) and mark status unknown (cats which were not black, but couldn't be identified to the individual level with complete certainty).

We discarded detections of cats which were obviously young enough to be dependent on a parent, as these kittens do not have independent activity centres or movements and not yet recruited into the adult population. 

![Feral cat coat categories.](figs/cat_coats.jpg)

\newpage


# Feral cat detection plots 

## Glenelg 

### Pair 1 
```{r, echo = FALSE, fig.cap = "Feral cat detections in the first grid pair of the Glenelg region. Camera-traps are indicated by black crosses. Circles with coloured fills represent individual cats and lines denote thier observed movements. Black circles represent black cats, red circles unidentifiable detections (the circle radius scales positively with the number daily detections). Fox control does not occur in Annya (A) but has since 2005 in Cobboboonee (B)."}

par(mfrow = c(1,2), mar = c(2.6, 1, 2.1, 1))

plot(mrch_a, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "A) Annya State Forest", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_a, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_a, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

plot(mrch_c, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "B) Cobboboonee National Park", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_c, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_c, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

```


\newpage


### Pair 2
```{r, echo = FALSE, fig.cap = "Feral cat detections in the second grid pair of the Glenelg region. Camera-traps are indicated by black crosses. Circles with coloured fills represent individual cats and lines denote thier observed movements. Black circles represent black cats, red circles unidentifiable detections (the circle radius scales positively with the number daily detections). Fox control does not occur in Hotspur (A) but has since 2005 in Mt Clay (B)."}

par(mfrow = c(1,2), mar = c(2.6, 1, 2.1, 1))

plot(mrch_h, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "A) Hotspur State Forest", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_h, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_h, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

plot(mrch_mc, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "B) Mt Clay Reserve", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_mc, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

```


\newpage




## Otway Ranges
### 2017
```{r, echo = FALSE,   fig.cap = "Feral cat detections in the Otway Ranges 2017. Camera-traps are indicated by black crosses. Circles with coloured fills represent individual cats and lines denote thier observed movements. Black circles represent black cats, red circles unidentifiable detections (the circle radius scales positively with the number daily detections). Fox control did not occur in either of the grids during this time."}

par(mfrow = c(1,2), mar = c(2.6, 1, 2.1, 1))

plot(mrch_n_17, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "A) Otways north 2017", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_n_17, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_n_17, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

plot(mrch_s_17, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "B) Otways south 2017", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_s_17, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_s_17, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

```


\newpage


### 2018
```{r, echo = FALSE,   fig.cap = "Feral cat detections in the Otway Ranges 2018. Camera-traps are indicated by black crosses. Circles with coloured fills represent individual cats and lines denote thier observed movements. Black circles represent black cats, red circles unidentifiable detections (the circle radius scales positively with the number daily detections)."}

par(mfrow = c(1,2), mar = c(2.6, 1, 2.1, 1))

plot(mrch_n_18, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "A) Otways north 2018", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_n_18, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_n_18, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

plot(mrch_s_18, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "B) Otways south 2018", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_s_18, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_s_18, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

```


\newpage


### 2019
```{r, echo = FALSE, fig.cap = "Feral cat detections in the Otway Ranges 2019. Camera-traps are indicated by black crosses. Circles with coloured fills represent individual cats and lines denote thier observed movements. Black circles represent black cats, red circles unidentifiable detections (the circle radius scales positively with the number daily detections)."}

par(mfrow = c(1,2), mar = c(2.6, 1, 2.1, 1))

plot(mrch_n_19, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "A) Otways north 2019", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_n_19, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_n_19, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

plot(mrch_s_19, rad = 2, tracks = T, bg = "white", detpar = list(col = "black"), title = "B) Otways south 2019", icolours = brewer.pal(12, "Paired"), randcol = TRUE, gridlines = FALSE, cappar = list(cex = 2), border = 10)
sightingPlot(mrch_s_19, type = "Tu", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "black")
sightingPlot(mrch_s_19, type = "Tn", add = TRUE, scale = 200, mean = FALSE, px = 1.1, py = 1.9, title = "", col = "red")

```


\newpage

# Fox and prey spatial layers

## Generalised additive model summaries

### Glenelg 

#### Fox
```{r, echo = FALSE}
summary(gam_g_fox)
```

#### Common brushtail possum
```{r, echo = FALSE}
summary(gam_g_btp)
```

#### Common ringtail possum
```{r, echo = FALSE}
summary(gam_g_rtp)
```

#### Long-nosed potoroo
```{r, echo = FALSE}
summary(gam_g_pot)
```

#### Southern brown bandicoot
```{r, echo = FALSE}
summary(gam_g_sbb)
```

#### Small mammals
```{r, echo = FALSE}
summary(gam_g_sm)
```

### Otways 

#### Fox
```{r, echo = FALSE}
summary(gam_o_fox)
```


#### Common brushtail possum
```{r, echo = FALSE}
summary(gam_o_btp)
```


#### Common ringtail possum
```{r, echo = FALSE}
summary(gam_o_rtp)
```


#### Long-nosed potoroo
```{r, echo = FALSE}
summary(gam_o_pot)
```


#### Bandicoots
The vast majority of bandicoots here were long-nosed bandicoots - only few sites detected southern brown bandicoots (on the edge of their habitat). We therefore combined these detections into a singular bandicoot group. 
```{r, echo = FALSE}
summary(gam_o_band)
```


#### Small mammals
```{r, echo = FALSE}
summary(gam_o_sm)
```


\newpage




# Density estimation 

## Glenelg region 

```{r, secr glenelg detectfn AIC, echo = FALSE}

df_fits_glenelg_aic <- AIC(df_fits_glenelg, criterion = "AICc")[,-2]
df_fits_glenelg_aic <- subset(df_fits_glenelg_aic, select = -c(model))
df_fits_glenelg_aic$logLik <- round(df_fits_glenelg_aic$logLik, digits = 2)
df_fits_glenelg_aic$AIC    <- round(df_fits_glenelg_aic$AIC, digits = 2)
df_fits_glenelg_aic$AICc   <- round(df_fits_glenelg_aic$AICc, digits = 2)
df_fits_glenelg_aic$dAICc   <- round(df_fits_glenelg_aic$dAICc, digits = 2)
df_fits_glenelg_aic$AICcwt  <- round(df_fits_glenelg_aic$AICcwt, digits = 2)

df_fits_glenelg_aic <- tibble::rownames_to_column(df_fits_glenelg_aic, "fit")
df_fits_glenelg_aic$fit <- ifelse(df_fits_glenelg_aic$fit == "fit1_ex", "exponential", "half-normal")

pander(df_fits_glenelg_aic, split.cell = 40, split.table = Inf, col.names = c("Detector function", "Parameters", "logLik", "AIC", "AICc", "dAICc", "AICcwt"), round = 3,
       caption = "Akaike's Information Criterion values for detector functions in the Glenelg region.")

```



```{r, secr glenelg fits AIC, echo = FALSE,  message = FALSE, warning = FALSE}

glenelg_fits_aic <- AIC(glenelg_fits, criterion = "AICc")[,-2]

glenelg_fits_aic$logLik <- round(glenelg_fits_aic$logLik, digits = 2)
glenelg_fits_aic$AIC <- round(glenelg_fits_aic$AIC, digits = 2)
glenelg_fits_aic$AICc <- round(glenelg_fits_aic$AICc, digits = 2)
glenelg_fits_aic$dAICc <- round(glenelg_fits_aic$dAICc, digits = 2)
glenelg_fits_aic$AICcwt <- round(glenelg_fits_aic$AICcwt, digits = 2)

glenelg_fits_aic$model <- gsub("_trapcov", "", glenelg_fits_aic$model)
glenelg_fits_aic$model <- gsub(" g0", ", g0", glenelg_fits_aic$model)
glenelg_fits_aic$model <- gsub(" sigma", ", sigma", glenelg_fits_aic$model)
glenelg_fits_aic$model <- gsub("~", " ~ ", glenelg_fits_aic$model)
glenelg_fits_aic$model <- gsub("_predicted", "", glenelg_fits_aic$model)
  
glenelg_fits_aic <- tibble::rownames_to_column(glenelg_fits_aic, "fit")

pander(glenelg_fits_aic, split.cell = 40, split.table = Inf, col.names = c("Fit", 
    "Model", "Parameters", "logLik", "AIC", "AICc", "dAICc", "AICcwt"), round = 2, caption = "Akaike's Information Criterion values for feral cat density models in the Glenelg Region (ordered in decreasing AICc scores).")

```

*Note:* We excluded fit13 (D ~ fox + prey, g0 ~ 1, sigma ~ 1) and fit15 (D ~ fox + prey, g0 ~ prey, sigma ~ prey) because they did not converge. 


```{r, echo = FALSE, fig.width = 8, fig.height=10, fig.cap = "The most complex model of feral cat density in the Glenelg region, Australia, as a function of the probability of fox occupancy (Panel A) and average prey occupancy (Panel B), as well as feral cat detectability; g0 (probability of detection per 24-hour occasion in the location of an individual’s activity centre; Panel C, E) and sigma (which is related to home range size; Panel D, F) as a function of fox (Panel C, D) and average prey (Panel E, F) probability of occupancy (fit16)."}

par(mar=c(5, 7, 4, 3.5), mfrow = c(3, 2), mgp = c(3, 1, 0))

# for figure labelling
# use the coords of the plot region as fractions of the figure region
# to define the adj= argument of mtext()
pplt <- par("plt")
adjx <- (0 - pplt[1]) / (pplt[2] - pplt[1])
# use the number of lines of margin to define the line= argument of
liney <- par("mar")[3] - 1.5

# D ~ FOX ACTIVITY ---------------------------------------------------------

# extract values
all_predicted <- predict(glenelg_fits$fit16, 
                         newdata = data.frame(fox_predicted = seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted)/30),
                                              prey_predicted = mean(glenelg_mask_df$prey_predicted), 
                                              prey_predicted_trapcov = mean(glenelg_mask_df$prey_predicted), 
                                              fox_predicted_trapcov = mean(glenelg_mask_df$fox_predicted)))
predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))
lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))
upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))
time_sequence <- seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted)/30)

# plot - FOX 
plot(predicted_values*100 ~ time_sequence, type = "l", lwd = 4.5, col = "black", ylim = c(0, 0.8), las = 1,  xlab = NA, ylab = NA, cex.axis = 1.25)
# CI as a shaded region
polygon(c(time_sequence, rev(time_sequence)), c(lower_bound*100, rev(upper_bound*100)),
        border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                                                                                                                                                                                                               col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values*100 ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Fox Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("Cats per km"^2), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("A", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     # D ~ PREY ACTIVITY ---------------------------------------------------------
                                                                                                                                                     # extract values
                                                                                                                                                     all_predicted <- predict(glenelg_fits$fit16, 
                                                                                                                                                                              newdata = data.frame(prey_predicted= seq(min(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted)/30),
                                                                                                                                                                                                   prey_predicted_trapcov = mean(glenelg_mask_df$prey_predicted), 
                                                                                                                                                                                                   fox_predicted = mean(glenelg_mask_df$fox_predicted),
                                                                                                                                                                                                   fox_predicted_trapcov = mean(glenelg_mask_df$fox_predicted)))
                                                                                                                                                     
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "D","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "D","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "D","ucl"))
                                                                                                                                                     time_sequence <- seq(min(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot - prey 
                                                                                                                                                     plot(predicted_values*100 ~ time_sequence, type = "l", lwd = 4.5, col = "black", ylim = c(0, 0.8), las = 1,  xlab = NA, ylab = NA, cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound*100, rev(upper_bound*100)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values*100 ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Average prey Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("Cats per km"^2), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("B", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     # DET ~ FOX  --------------------------------------------------
                                                                                                                                                     
                                                                                                                                                     # extract values
                                                                                                                                                     all_predicted <- predict(glenelg_fits$fit16, 
                                                                                                                                                                              newdata = data.frame(fox_predicted_trapcov = seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted)/30),
                                                                                                                                                                                                   prey_predicted = mean(glenelg_mask_df$prey_predicted), 
                                                                                                                                                                                                   fox_predicted = mean(glenelg_mask_df$fox_predicted),
                                                                                                                                                                                                   prey_predicted_trapcov = mean(glenelg_mask_df$prey_predicted)))
                                                                                                                                                     
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "g0","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "g0","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "g0","ucl"))
                                                                                                                                                     time_sequence <- seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 1), las = 1,  xlab = NA, ylab = NA,  cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Fox Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("g0"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("C", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     ## SIGMA
                                                                                                                                                     # extract values
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "sigma","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "sigma","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "sigma","ucl"))
                                                                                                                                                     time_sequence <- seq(min(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted), max(glenelg_mask_df$fox_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 2000), las = 1,  xlab = NA, ylab = NA, cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Fox Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("sigma"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("D", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     # DET ~ PREY  --------------------------------------------------
                                                                                                                                                     
                                                                                                                                                     # extract values
                                                                                                                                                     all_predicted <- predict(glenelg_fits$fit16, 
                                                                                                                                                                              newdata = data.frame(prey_predicted_trapcov = seq(min(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted)/30),
                                                                                                                                                                                                   prey_predicted = mean(glenelg_mask_df$prey_predicted), 
                                                                                                                                                                                                   fox_predicted = mean(glenelg_mask_df$fox_predicted),
                                                                                                                                                                                                   fox_predicted_trapcov = mean(glenelg_mask_df$fox_predicted)))
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "g0","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "g0","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "g0","ucl"))
                                                                                                                                                     time_sequence <- seq(min(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 1), las = 1,  xlab = NA, ylab = NA,  cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Average prey Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("g0"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("E", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     ## SIGMA
                                                                                                                                                     # extract values
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "sigma","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "sigma","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "sigma","ucl"))
                                                                                                                                                     time_sequence <- seq(min(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted), max(glenelg_mask_df$prey_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 2000), las = 1,  xlab = NA, ylab = NA, cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Average prey Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("sigma"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("F", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     ```
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     \newpage
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     ## Otway Ranges 
                                                                                                                                                     
                                                                                                                                                     ```{r, secr otways detectfn AIC, echo = FALSE}
                                                                                                                                                     
                                                                                                                                                     df_fits_otways_aic <- AIC(df_fits_otways, criterion = "AICc")[,-2]
                                                                                                                                                     df_fits_otways_aic <- subset(df_fits_otways_aic, select = -c(model))
                                                                                                                                                     df_fits_otways_aic$logLik <- round(df_fits_otways_aic$logLik, digits = 2)
                                                                                                                                                     df_fits_otways_aic$AIC    <- round(df_fits_otways_aic$AIC, digits = 2)
                                                                                                                                                     df_fits_otways_aic$AICc   <- round(df_fits_otways_aic$AICc, digits = 2)
                                                                                                                                                     df_fits_otways_aic$dAICc   <- round(df_fits_otways_aic$dAICc, digits = 2)
                                                                                                                                                     df_fits_otways_aic$AICcwt  <- round(df_fits_otways_aic$AICcwt, digits = 2)
                                                                                                                                                     
                                                                                                                                                     df_fits_otways_aic <- tibble::rownames_to_column(df_fits_otways_aic, "fit")
                                                                                                                                                     df_fits_otways_aic$fit <- ifelse(df_fits_otways_aic$fit == "fit1_ex", "exponential", "half-normal")
                                                                                                                                                     
                                                                                                                                                     pander(df_fits_otways_aic, split.cell = 40, split.table = Inf, col.names = c("Detector function", "Parameters", "logLik", "AIC", "AICc", "dAICc", "AICcwt"), round = 3, caption = "Akaike's Information Criterion values for detector functions in the Otway Ranges.")
                                                                                                                                                     
                                                                                                                                                     ```
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     \newpage
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     ```{r, secr otways fits AIC, echo = FALSE,  message = FALSE, warning = FALSE}
                                                                                                                                                     
                                                                                                                                                     otways_fits_aic <- AIC(otways_fits, criterion = "AICc")[,-2]
                                                                                                                                                     
                                                                                                                                                     otways_fits_aic$logLik <- round(otways_fits_aic$logLik, digits = 2)
                                                                                                                                                     otways_fits_aic$AIC <- round(otways_fits_aic$AIC, digits = 2)
                                                                                                                                                     otways_fits_aic$AICc <- round(otways_fits_aic$AICc, digits = 2)
                                                                                                                                                     otways_fits_aic$dAICc <- round(otways_fits_aic$dAICc, digits = 2)
                                                                                                                                                     otways_fits_aic$AICcwt <- round(otways_fits_aic$AICcwt, digits = 2)
                                                                                                                                                     
                                                                                                                                                     otways_fits_aic$model <- gsub("_trapcov", "", otways_fits_aic$model)
                                                                                                                                                     otways_fits_aic$model <- gsub(" g0", ", g0", otways_fits_aic$model)
                                                                                                                                                     otways_fits_aic$model <- gsub(" sigma", ", sigma", otways_fits_aic$model)
                                                                                                                                                     otways_fits_aic$model <- gsub("~", " ~ ", otways_fits_aic$model)
                                                                                                                                                     otways_fits_aic$model <- gsub("_predicted", "", otways_fits_aic$model)
                                                                                                                                                     
                                                                                                                                                     otways_fits_aic <- tibble::rownames_to_column(otways_fits_aic, "fit")
                                                                                                                                                     otways_fits_aic$fit <- otways_fits_aic$fit %>% stringr::str_remove(pattern = ".*fit")
                                                                                                                                                     otways_fits_aic$fit <- paste0("fit", otways_fits_aic$fit)
                                                                                                                                                     
                                                                                                                                                     pander(otways_fits_aic, split.cell = 40, split.table = Inf, col.names = c("Fit", 
                                                                                                                                                                                                                               "Model", "Parameters", "logLik", "AIC", "AICc", "dAICc", "AICcwt"), round = 3, caption = "Akaike's Information Criterion values for feral cat density models in the Otways Ranges (ordered in decreasing AICc scores).")
                                                                                                                                                     
                                                                                                                                                     ```
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     ```{r, echo = FALSE, fig.width = 8, fig.height=10, fig.cap = "The best-ranked model for the Otway region, Australia, Australia, predicted that feral cat detectability; g0 (probability of detection per 24-hour occasion in the location of an individual’s activity centre; row 1) and sigma (which is related to home range size; row 2) where impacted by fox (column 1) and prey (column 2) occupancy probability. Shaded region indicates 95% confidence intervals."}
                                                                                                                                                     
                                                                                                                                                     par(mar=c(5, 7, 4, 3.5), mfrow = c(3, 2), mgp = c(3, 1, 0))
                                                                                                                                                     
                                                                                                                                                     # for figure labelling
                                                                                                                                                     # use the coords of the plot region as fractions of the figure region
                                                                                                                                                     # to define the adj= argument of mtext()
                                                                                                                                                     pplt <- par("plt")
                                                                                                                                                     adjx <- (0 - pplt[1]) / (pplt[2] - pplt[1])
                                                                                                                                                     # use the number of lines of margin to define the line= argument of
                                                                                                                                                     liney <- par("mar")[3] - 1.5
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     # g0  --------------------------------------------------
                                                                                                                                                     ## FOXES
                                                                                                                                                     # extract values
                                                                                                                                                     all_predicted <- predict(otways_fits$fit28, 
                                                                                                                                                                              newdata = data.frame(fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted)/30),
                                                                                                                                                                                                   prey_predicted = mean(otways_mask_df$prey_predicted), 
                                                                                                                                                                                                   fox_predicted = mean(otways_mask_df$fox_predicted),
                                                                                                                                                                                                   prey_predicted_trapcov = mean(otways_mask_df$prey_predicted)))
                                                                                                                                                     
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "g0","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "g0","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "g0","ucl"))
                                                                                                                                                     time_sequence <- seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 0.25), las = 1,  xlab = NA, ylab = NA,  cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Fox Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("g0"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("A", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     ## PREY
                                                                                                                                                     # extract values
                                                                                                                                                     all_predicted <- predict(otways_fits$fit28, 
                                                                                                                                                                              newdata = data.frame(prey_predicted_trapcov = seq(min(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted)/30),
                                                                                                                                                                                                   prey_predicted = mean(otways_mask_df$prey_predicted), 
                                                                                                                                                                                                   fox_predicted = mean(otways_mask_df$fox_predicted),
                                                                                                                                                                                                   fox_predicted_trapcov = mean(otways_mask_df$fox_predicted)))
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "g0","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "g0","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "g0","ucl"))
                                                                                                                                                     time_sequence <- seq(min(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 0.25), las = 1,  xlab = NA, ylab = NA,  cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Average prey Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("g0"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("B", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     # SIGMA  --------------------------------------------------
                                                                                                                                                     ## FOXES
                                                                                                                                                     # extract values
                                                                                                                                                     all_predicted <- predict(otways_fits$fit28, 
                                                                                                                                                                              newdata = data.frame(fox_predicted_trapcov = seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted)/30),
                                                                                                                                                                                                   prey_predicted = mean(otways_mask_df$prey_predicted), 
                                                                                                                                                                                                   fox_predicted = mean(otways_mask_df$fox_predicted),
                                                                                                                                                                                                   prey_predicted_trapcov = mean(otways_mask_df$prey_predicted)))
                                                                                                                                                     
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "sigma","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "sigma","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "sigma","ucl"))
                                                                                                                                                     time_sequence <- seq(min(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted), max(otways_mask_df$fox_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 600), las = 1,  xlab = NA, ylab = NA, cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Fox Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("sigma"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("C", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     # DET ~ PREY  --------------------------------------------------
                                                                                                                                                     
                                                                                                                                                     # FOXES
                                                                                                                                                     all_predicted <- predict(otways_fits$fit28, 
                                                                                                                                                                              newdata = data.frame(prey_predicted_trapcov = seq(min(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted)/30),
                                                                                                                                                                                                   prey_predicted = mean(otways_mask_df$prey_predicted), 
                                                                                                                                                                                                   fox_predicted = mean(otways_mask_df$fox_predicted),
                                                                                                                                                                                                   fox_predicted_trapcov = mean(otways_mask_df$fox_predicted)))
                                                                                                                                                     
                                                                                                                                                     # extract values
                                                                                                                                                     predicted_values <- unlist(sapply(all_predicted, "[", "sigma","estimate"))
                                                                                                                                                     lower_bound <- unlist(sapply(all_predicted, "[", "sigma","lcl"))
                                                                                                                                                     upper_bound <- unlist(sapply(all_predicted, "[", "sigma","ucl"))
                                                                                                                                                     time_sequence <- seq(min(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted), max(otways_mask_df$prey_predicted)/30)
                                                                                                                                                     
                                                                                                                                                     # plot
                                                                                                                                                     plot(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black",  ylim = c(0, 600), las = 1,  xlab = NA, ylab = NA, cex.axis = 1.25)
                                                                                                                                                     # CI as a shaded region
                                                                                                                                                     polygon(c(time_sequence, rev(time_sequence)), c(lower_bound, rev(upper_bound)),
                                                                                                                                                             border = NA,      # NA here means don't include a border around the edge
                                                                                                                                                             col = "gray85")   # gray50 is just a 50% black/white mix, can change to any colour
                                                                                                                                                     # redraw the mean line
                                                                                                                                                     lines(predicted_values ~ time_sequence, type = "l", lwd = 4.5, col = "black")
                                                                                                                                                     mtext(expression("Average prey Pr(occupancy)"), side = 1, las = 1, line = 4, cex = 1)
                                                                                                                                                     mtext(expression("sigma"), side = 2, las = 3, line = 4, cex = 1)
                                                                                                                                                     mtext("D", side=3, adj=adjx, line=liney, cex = 1)
                                                                                                                                                     
                                                                                                                                                     ```
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     