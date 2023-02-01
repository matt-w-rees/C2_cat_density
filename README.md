This repository contains data and code used for the analysis of this manuscript. 

# Mesopredator release among invasive predators: controlling red foxes can increase feral cat density and alter their behaviour
### Matthew W. Rees, Jack H. Pascoe, Mark Le Pla, Alan Robley, Emma K. Birnbaum, Brendan A. Wintle, Bronwyn A. Hradsky

All analyses were conducted in R, numbered scripts (reflecting the order to run scripts in) are contained in the 'r_scripts' folder. These scripts format the raw data and fit models (there are also two scripts to generate the plots made for the manuscript).

The 'raw-data' folder contains:

'tagged_cat_images': raw images of feral cats with mark status and individuals identifications tagged to the image metadata
Lower Glenelg cat capture-histories (this data was extracted separately to the other landscapes). 
Spatial files of Vegetation types ('vegetation_layer/') and the Lower Glenelg National Park ('lower_glenelg_shp/') used for habitat mask input data and plots. 
Camera-trap presence-absence records (to fit fox occurrence GAMs): 'spp_records_pa.csv' and 'lower_glenelg_fox_pa.csv'. 
The 'derived-data/' folder contains secr habitat masks and capture histories produced in R scripts 1:3. 

The 'model/' folder contains saved secr models produced in R scripts 4:5. 