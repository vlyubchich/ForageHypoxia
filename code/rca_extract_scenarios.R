# Extract data from ROMS-RCA model outputs saved on Jeremy's cluster.
# Specifically, save a table with cell information (lat, lon, dimensions, depth, type)
# and daily summaries of water quality variables for water cells
# (one file of the summaries per year).
# These resulting files should be placed in our folder ./data_rca/

# cd /local/users/cshen/rca_30year/
cd /local/users/adreiss/tst/
R

# Packages and functions ----
rm(list = ls())
library(abind)
library(ncdf4)
library(data.table)
library(zoo)

VER <- c(# "p2000_10pct_inc",
         # "p2000_10pct_redn",
         # "p2000_25pct_inc",
         # "p2000_25pct_redn",
         # "p2000_50pct_inc",
         # "p2000_50pct_redn",
         # "p2000_1D",
         # "p2000_2D",
         # "p2000_3D",
         # "p2000_4D",
         # "p2000_2D_50in",
         # "p2000_2D_75redn",
         # "p2000_4D_50in",
         # "p2000_4D_75redn",
    "p2000")

mysum <- function(x) {
    c(mean(x), sd(x)) # min(x), max(x),
}
rcadata <- function(X, sel, dt, idt) {
    # X = WTEMP; sel = selection; idt = idates
    # X = input array
    # sel = spatial selections
    # dt = vector of dates
    # idt = indices for date subset
    vname <- deparse(substitute(X))
    if (dim(sel)[2] == 2) { # if select only one layer at a time
        D <- apply(sel, 1, function(i) { # i = sel[1000,]
            data.table(zoo::rollapply(X[i[1], i[2], 1, idt],
                                      width = 6, mysum, by = 6),
                       CellID = CID[i[1], i[2]],
                       Date = dt)
        })
        D <- data.table::rbindlist(D) # do.call(rbind, D)
        names(D)[1:2] <- paste(vname, c("avg", "sd"), sep = "_") # "min", "max",
    }
    return(D)
}

is.leap <- function(year) {
    if ((year %% 4) == 0) {
        if ((year %% 100) == 0) {
            if ((year %% 400) == 0) {
                res = TRUE
            } else {
                res = FALSE
            }
        } else {
            res = TRUE
        }
    } else {
        res = FALSE
    }
    return(res)
}

# Settings ----
# Set what data will extract
YEARS <- 2000L
LAYERS <- 20 # layer to extract
files_exclude <- c("./1987_out1//Y1987_eutr_0146.nc", # no time dimensionality
                   "./2006_out1//Y2006_eutr_0840.nc",
                   "./2008_out1//Y2008_eutr_0913.nc",
                   "./2012_out1//Y2012_eutr_1059.nc"
)

for (ver in VER) {
    for (year in YEARS) { # year = 2000; ver = "p2000_10pct_inc"
        # List files YXXXX_eutr_XXXX.nc, where X are the numbers
        files_year <- list.files(paste0("./", ver, "/"),
                                 pattern = "eutr\\_\\d+",
                                 full.names = TRUE)
        files_year <- setdiff(files_year, files_exclude)

        for (ify in 1:length(files_year)) { # ify = 1; fy =  "./p2000_10pct_inc//Y2000_eutr_0621.nc"
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
                          file = paste0("/local/users/lyubchich/rca_cells_", ver, ".csv"))
                # spatial selection of cells
                selection <- which(FSM == 1, arr.ind = TRUE)
                # dim(selection) #  3386    2 # checked for 1986 and 1989
            }

            if (ify == 1) {
                # Create "starting" datasets to combine with data from each file, per year
                Dates <- ncvar_get(nf, varid = "TIME")
                WTEMP <- ncvar_get(nf, varid = "HYDTEMP")[,, LAYERS,, drop = FALSE]
                SAL <- ncvar_get(nf, varid = "SAL")[,, LAYERS,, drop = FALSE]
                TPOC <- ncvar_get(nf, varid = "TPOC")[,, LAYERS,, drop = FALSE]
                CHLAVEG <- ncvar_get(nf, varid = "CHLAVEG")[,, LAYERS,, drop = FALSE]
                DOAVEG <- ncvar_get(nf, varid = "DOAVEG")[,, LAYERS,, drop = FALSE]
            } else {
                # Extract dates
                Dates <- abind(Dates, ncvar_get(nf, varid = "TIME"))
                # Extract water quality variables
                # water temperature (degrees C)
                WTEMP <- abind(WTEMP, ncvar_get(nf, varid = "HYDTEMP")[,, LAYERS,, drop = FALSE], along = 4)
                # salinity
                SAL <- abind(SAL, ncvar_get(nf, varid = "SAL")[,, LAYERS,, drop = FALSE], along = 4)
                # particulate organic carbon (mg/L)
                TPOC <- abind(TPOC, ncvar_get(nf, varid = "TPOC")[,, LAYERS,, drop = FALSE], along = 4)
                # chlorophyll-a (ug/L)
                CHLAVEG <- abind(CHLAVEG, ncvar_get(nf, varid = "CHLAVEG")[,, LAYERS,, drop = FALSE], along = 4)
                # dissolved oxygen (mg/L)
                DOAVEG <- abind(DOAVEG, ncvar_get(nf, varid = "DOAVEG")[,, LAYERS,, drop = FALSE], along = 4)
            }
            nc_close(nf)
        }

        # To get calendar dates, add the TIME variable to 12 AM January 1 1983.
        Dates <- Dates + as.Date("1983-01-01")
        check_len_Dates <- length(Dates)
        # summary(Dates)
        # The RCA model has no leap years, remove Feb 29
        if (is.leap(year)) {
            f29 <- as.Date(paste0(year, "-02-29"))
            Dates[Dates >= f29] <- Dates[Dates >= f29] + 1
        }
        # Dates index to use
        idates <- !(Dates < as.Date(paste0(year, "-01-01")) |
                        Dates >= as.Date(paste0(year + 1, "-01-01")))
        # Check that each date has 6 records
        date_count <- table(Dates[idates])
        if (date_count[1] != 6) {
            date_count <- date_count[-1]
        }
        if (date_count[length(date_count)] != 6) {
            date_count <- date_count[-length(date_count)]
        }
        if (!all(date_count == 6)) {
            stop(paste0("Missing some obs in the middle of year ", year))
        }
        # Refresh the dates index
        idates <- is.element(as.character(Dates), names(date_count))
        Dates <- names(date_count)
        if (length(Dates) < 365) {
            print(paste0(year, " year is not full, some days are missing"))
        }
        if (dim(DOAVEG)[4] != check_len_Dates) {
            stop(paste0("Mismatch of the length of date indices (", check_len_Dates,
                        ") and dataset dimension (", dim(DOAVEG)[4], ")"))
        }

        # Process and combine data -- all columns for one variable,
        # remove index columns for other variables
        D <- cbind(rcadata(X = WTEMP, sel = selection, dt = Dates, idt = idates)[,1:2],
                   rcadata(X = SAL, sel = selection, dt = Dates, idt = idates)[,1:2],
                   rcadata(X = TPOC, sel = selection, dt = Dates, idt = idates)[,1:2],
                   rcadata(X = CHLAVEG, sel = selection, dt = Dates, idt = idates)[,1:2],
                   rcadata(X = DOAVEG, sel = selection, dt = Dates, idt = idates))

        write.csv(D, row.names = FALSE,
                  file = paste0("/local/users/lyubchich/rca_ts_", year, "_", ver, ".csv"))
    }
}
