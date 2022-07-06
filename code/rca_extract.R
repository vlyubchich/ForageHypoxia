
cd /local/users/cshen/rca_30year/
R

# Packages and functions ----
rm(list = ls())
library(ncdf4)
library(data.table)
library(zoo)

mysum <- function(x) {
    c(mean(x), min(x), max(x), sd(x))
}
rcadata <- function(X, sel, dt, idt) {
    # X = WTEMP; sel = selection; idt = idates
    # X = input array
    # sel = selections
    # idt = indices for date subset
    vname <- deparse(substitute(X))
    if (dim(sel)[2] == 2) { # if select only one layer at a time
        D <- apply(sel, 1, function(i) { # i = sel[1000,]
            data.table(zoo::rollapply(X[i[1], i[2], LAYERS, idt],
                                      width = 6, mysum, by = 6),
                       CellID = CID[i[1], i[2]],
                       Date = dt)
        })
        D <- data.table::rbindlist(D) # do.call(rbind, D)
        names(D)[1:4] <- paste(vname, c("avg", "min", "max", "sd"), sep = "_")
    }
    return(D)
}

# Settings ----
# Set what data will extract
YEARS <- 2012:2015 # 30 years 1986:2015
LAYERS <- 20 # layers to extract
files_exclude <- c("./1987_out1//Y1987_eutr_0146.nc", # no time dimensionality
                   "./2006_out1//Y2006_eutr_0840.nc",
                   "./2008_out1//Y2008_eutr_0913.nc",
                   "./2012_out1//Y2012_eutr_1059.nc"
)

for (year in YEARS) { # year = 1989
    D <- data.table()

    # List files YXXXX_eutr_XXXX.nc, where X are the numbers
    files_year <- list.files(paste0("./", year, "_out1/"), full.names = TRUE)
    files_year <- sort(files_year[grepl("eutr\\_\\d+", files_year)])
    files_year <- setdiff(files_year, files_exclude)
    for (ify in 1:length(files_year)) { # fy = "./1986_out1//Y1986_eutr_0110.nc" # fy = files_year[1]
        fy <- files_year[ify]
        print(fy)
        nf <- nc_open(fy)
        # print(nf)

        # If this is the first file, save locations (cell) information
        if (year == YEARS[1] && ify == 1) {
            H <- ncvar_get(nf, varid = "H") # mean water depth (meter)
            DX <- ncvar_get(nf, varid = "DX") # length of each water cell (meter)
            DY <- ncvar_get(nf, varid = "DY") # width of each water cell (meter)
            LAT <- ncvar_get(nf, varid = "LAT") # latitude at cell center (degree)
            LON <- ncvar_get(nf, varid = "LON") # longitude at cell center (degree)
            FSM <- ncvar_get(nf, varid = "FSM") # land mask, 1=water, 0=land, -1=river BC, -2=ocean BC
            CID <- array(1:prod(dim(H)), dim = dim(H)) # cell IDs in the matrix format
            CELLS <- data.frame(CellID = as.vector(CID),
                                H = as.vector(H),
                                DX = as.vector(DX),
                                DY = as.vector(DY),
                                LAT = as.vector(LAT),
                                LON = as.vector(LON),
                                FSM = as.vector(FSM))
            write.csv(CELLS, row.names = FALSE,
                      file = paste0("/local/users/lyubchich/rca_cells_", Sys.Date(), ".csv"))
            # selection of cells
            selection <- which(FSM == 1, arr.ind = TRUE)
            # dim(selection) #  3386    2 # checked for 1986 and 1989
        }

        # Dates from the current file
        # To get a calendar date, add the TIME variable to 12 AM January 1 1983.
        Dates <- ncvar_get(nf, varid = "TIME") + as.Date("1983-01-01") + 0.84
        # dates index to use
        idates <- !(Dates < as.Date(paste0(year, "-01-01")))
        Dates <- names(table(Dates[idates]))

        # Extract the actual variables
        WTEMP <- ncvar_get(nf, varid = "HYDTEMP") # water temperature (degrees C)
        SAL <- ncvar_get(nf, varid = "SAL") # salinity
        TPOC <- ncvar_get(nf, varid = "TPOC") # particulate organic carbon (mg/L)
        CHLAVEG <- ncvar_get(nf, varid = "CHLAVEG") # chlorophyll-a (ug/L)
        DOAVEG <- ncvar_get(nf, varid = "DOAVEG") # dissolved oxygen (mg/L)

        # Process and combine data -- all columns for one variable,
        # remove index columns for other variables
        d <- cbind(rcadata(X = WTEMP, sel = selection, dt = Dates, idt = idates)[,1:4],
                   rcadata(X = SAL, sel = selection, dt = Dates, idt = idates)[,1:4],
                   rcadata(X = TPOC, sel = selection, dt = Dates, idt = idates)[,1:4],
                   rcadata(X = CHLAVEG, sel = selection, dt = Dates, idt = idates)[,1:4],
                   rcadata(X = DOAVEG, sel = selection, dt = Dates, idt = idates))
        D <- rbind(D, d)
        nc_close(nf)
    }
    write.csv(D, row.names = FALSE,
          file = paste0("/local/users/lyubchich/rca_data_", year, "_", Sys.Date(), ".csv"))
}
