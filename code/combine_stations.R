
Dcells <- read.csv("./data_rca/rca_cells_2022-06-29.csv")
Dben <- read.csv("./data_benthos/benthos_biomass.csv")

sapply(unique(Dben$Year), function(y) length(unique(Dben$STATION[Dben$Year == y])))


# base R operations
Dben$CellID <- NA #create a new empty variable
Dcells = Dcells[Dcells$FSM == 1,] #select only water cells

# dplyr operations
library(dplyr)
Dcells2 <- read.csv("./data_rca/rca_cells_2022-06-29.csv") %>%
    filter(FSM == 1)



# https://www.r-bloggers.com/2019/10/geographic-distance/

# find geodesic distances from each missing airport to non-missing
# and select the closest non-missing as a Replacement
for (i in 1:nrow(Airports_replace)) { # i = 1
    # distance from this airport to existing (km)
    d <- geosphere::distm(Airports_replace[i, c("lon", "lat")],
                          Airports[, c("lon", "lat")], fun = distGeo) / 1000L
    di <- which.min(d)
    Airports_replace$DistanceToNearestKm[i] <- d[di]
    Airports_replace$Replacement[i] <- Airports$Airport_Code[di]
}
