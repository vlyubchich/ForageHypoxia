library(dplyr)
library(geosphere)

benthomass <- read.csv("./data_benthos/benthos_biomass.csv") #data about segments (random only)
taxa_Data <- read.csv("./data_benthos/Taxa_group_IDs_updated_rjw_v2.xlsx - Invert_IDs.csv") #taxa data

df <- filter(benthomass, SITE_TYPE == 'RANDOM') #df stands for the main dataframe: all RANDOM sites.
df
