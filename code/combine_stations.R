
Dcells <- read.csv("./data_rca/rca_cells_2022-06-29.csv")
Dben <- read.csv("./data_benthos/Station_Data.csv")

# base R operations
Dben$CellID <- NA #create a new empty variable
Dcells = Dcells[Dcells$FSM == 1,] #select only water cells

# dplyr operations
library(dplyr)
library(geosphere)
Dcells2 <- read.csv("./data_rca/rca_cells_2022-06-29.csv") %>%
    filter(FSM == 1)


# https://www.r-bloggers.com/2019/10/geographic-distance/

# find geodesic distances from each cell to station
# and select the closest non-missing as a Replacement
for (i in 1:nrow(Dben)) { # i = 1
    # distance from this airport to existing (km)
    d <- geosphere::distm(Dben[i, c("LONGITUDE", "LATITUDE")],
                          Dcells[, c("LON", "LAT")])
    di <- which.min(d)
    Dben$DistanceToNearestKm[i] <- d[di]
    Dben$Replacement[i] <- Dcells$CellID[di]
}
'Data compilation complete!'
write.csv(Dben, './Cell_Proximity.csv', row.names = FALSE)
