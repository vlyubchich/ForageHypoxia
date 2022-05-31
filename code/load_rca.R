# 2022-05-27

#Load Ches Bay shape file
shoreline <- readShapePoly(fn = "./CBnet_shoreline/CBnet_Shoreline_simplify")

library(ncdf4)
library(maptools)
library(rgdal)
library(sp)

grdf <- nc_open("./data_rca/grid_new.nc")
print(grdf)
names(grdf$var) # use "lat_rho", "lon_rho" as Jeremy suggested
lats <- ncvar_get(grdf, "lat_rho")
lons <- ncvar_get(grdf, "lon_rho")

STATIONS <- data.frame(Longitude = as.vector(lons),
                       Latitude = as.vector(lats))
#Use info from the shape file to find projections for the stations
###http://www.prj2epsg.org/search
###26918 - NAD_1983_UTM_Zone_18N
###http://spatialreference.org/ref/epsg/26918/
tmp <- project(cbind(360 + STATIONS$Longitude, STATIONS$Latitude),
               proj = "+proj=utm +zone=18 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
STATIONS <- cbind(STATIONS, PointX = tmp[,1], PointY = tmp[,2])
coordinates(STATIONS) <- c("PointX", "PointY")


plot(x = lons, y = lats)
plot(STATIONS, axes = TRUE, pch = 1) #, bg = "white", axes = FALSE, cex = 3, pch = 21, col = "black"
plot(shoreline, add = TRUE, border = "blue", lwd = 2)


cd /local/users/jtesta/RCA_output_benthos/
R

library(R.matlab)

data <- readMat("RCAoutput.mat")


cd /local/users/cshen/rca_30year/1987_out1


library(ncdf4)

setwd("./data_rca/")

nf <- nc_open("Y1987_eutr_0146.nc")
print(nf)


data <- ncvar_get(nf)
print("here is the data in the file:")
print(data)
nc_close( ncold )
