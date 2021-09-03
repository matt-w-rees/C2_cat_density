# Quantifying mesopredator release: lethal control of an invasive apex predator alters feral cat density and detectability

### Matthew W. Rees, Jack H. Pascoe, Mark Le Pla, Alan Robley, Emma K. Birnbaum, Brendan A. Wintle, Bronwyn A. Hradsky


This repository contains data and code used for the analysis of this submitted manuscript. 

## Abstract
The mesopredator release hypothesis predicts that subordinate predator density will increase as apex predators decline. We replicated field experiments across two regions with a simple predator guild comprising the introduced red fox *Vulpes vulpes* and feral cat *Felis catus*. We identified 160 individual cats from 68,504 camera-trap nights and estimated cat density using spatial mark-resight models. Targeted lethal fox control was associated with a negligible to 3.7-fold increase in feral cat density, mirroring variation in the duration and intensity of fox suppression. Correlative models confirmed this--feral cat density was negatively associated with fox occurrence at a fine spatial scale. We also observed changes in feral cat detectability across the (artificial) apex predator activity gradient. Our results suggest integrated predator management would help protect shared native prey and highlights that mesopredator release can manifest as changes in both behaviour and density, distorting inference if these processes are not distinguished.

## Running the code
All analyses were conducted in R, numbered scripts (reflecting the order to run scripts in) are contained in the 'r_scripts' folder. 
The 'raw-data' folder contains camera-trap species records, shapefiles and Lower Glenelg cat capture-histories (this data was extracted separately). Individual cat images are also contained in this repository, but not pushed to github (due to size limits).
The 'derived-data' folder contains secr habitat masks and capture histories produced in R scripts 1:3. 
Glenelg and Otway secr models are fitted separately (script 4 and 5, respectively - but there is no need to run script 4 before 5). 
