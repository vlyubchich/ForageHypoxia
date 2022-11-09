# Match benthic stations with corresponding nearest RCA model cell
rm(list = ls())

library(dplyr)
library(ggplot2)
library(geosphere)

# Load model cells, only water (FSM == 1)
Dcells <- read.csv("./data_rca/rca_cells_2.csv") %>%
    filter(FSM == 1)

# Load benthic stations, only random locations, aggregate to unique stations
Dstations <- read.csv("./data_benthos/benthos_biomass.csv") %>%
    filter(SITE_TYPE == "RANDOM") %>%
    group_by(STATION, Stratum, SAMPLE_DATE, Year, LATITUDE, LONGITUDE) %>% #
    summarise(Depth = mean(TOTAL_DEPTH)) %>%
    ungroup()

# Number of stations sampled per year
table(Dstations$Year)
# 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011
# 136  250  250  250  250  250  250  250  250  250  250  250  250  250  250  250  250
# 2012 2013 2014 2015 2016 2017 2018 2019
# 250  250  250  250  250  250  250  250

# Remove stations outside the model domain
# Here could use a spatial boundary of the model to cut out
# stations outside the model S-shaped domain. Alternatively,
# use ad-hoc approach of removing stations after calculating
# the distances. It might keep some stations that are really
# close to the boundary.

# Find geodesic distances from each station to cells and select the closest cell
# https://www.r-bloggers.com/2019/10/geographic-distance/
Dstations$CellID <- Dstations$Distance_km <- Dstations$Depth_cell <- NA
for (i in 1:nrow(Dstations)) { # i = 1
    # distance from this station to each cell, km
    d <- geosphere::distm(Dstations[i, c("LONGITUDE", "LATITUDE")],
                          Dcells[, c("LON", "LAT")]) / 1000L
    # index of the cell that is the closest to the station
    di <- which.min(d)
    # corresponding distance (smallest from all)
    Dstations$Distance_km[i] <- d[di]
    # corresponding CellID
    Dstations$CellID[i] <- Dcells$CellID[di]
    # corresponding depth recorded in the model
    Dstations$Depth_cell[i] <- Dcells$H[di]
}

# Remove stations outside the model domain
distance_threshold <- 10L
Dstations$Remove <- Dstations$Distance_km >= distance_threshold

# Percent flagged for removal
mean(Dstations$Remove) * 100L
# [1] 9.110169

p1 <- ggplot(Dstations, aes(x = LONGITUDE, y = LATITUDE, color = Remove)) +
    geom_point() +
    ggtitle(paste0("Stations ", distance_threshold,
                   " km or more from any grid cell\nmarked for removal")) +
    theme_light()

p2 <- ggplot(Dstations, aes(x = Distance_km)) +
    geom_histogram(binwidth = 1L) +
    ggtitle("Distances between stations and nearest cells") +
    ylab("Count") +
    theme_light()

png("./images/Stations for removal.png",
    width = 4, height = 5,
    units = "in", res = 300)
print(p1)
dev.off()

png("./images/Stations to cells distance.png",
    width = 5, height = 5,
    units = "in", res = 300)
print(p2)
dev.off()

# Finally remove those stations
Dstations <- Dstations %>%
    filter(!Remove) %>%
    select(-Remove)

write.csv(Dstations,
          './data_benthos/stations_cells.csv',
          row.names = FALSE)
