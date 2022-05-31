# Extract benthic data, organize and save for easy reloading

# Packages ----
rm(list = ls()) #clean the environment
library(dplyr)
library(readxl)

# Functions ----
`%notin%` <- Negate(`%in%`)


# Load ----
# The data are from another project, on forage, that's why absolute paths.
# location of taxa data, use sheet 1 with data
pathtax <- "G:/My Drive/ResearchProjects/forage/dataraw/Taxa_group_IDs_updated_rjw_v2.xlsx"
# location of biomass data
pathfor <- "G:/My Drive/ResearchProjects/forage/dataraw/foraging_database.xlsx"
sheetsfor <- readxl::excel_sheets(pathfor)
# vector of options what an Excel cell with missing data can look like
nas <- c("", ".")


## Strata ----
STRATA <- data.frame(
    StratumFull = c("Mid Bay Mainstem", "Eastern Tributaries", "Western Tributaries",
                    "Upper Bay Mainstem", "Patuxent River", "Potomac River", "Mainstem",
                    "Rappahannock River", "York River", "James River"),
    State = c(rep("MD", 6), rep("VA", 4)),
    Area_km2 = c(2552, 534, 292, 785, 128, 1276, 4120, 372, 187, 684),
    Stratum = c("MMS", "MET", "MWT", "UPB", "PXR", "PMR",
                "VBY", "RAP", "YRK", "JAM")
)
STRATA <- STRATA %>%
    as_tibble() %>%
    mutate(AreaWeight = Area_km2 / sum(Area_km2))
# save for reuse
write.csv(STRATA,
          file = "./data_benthos/benthos_strata.csv",
          row.names = FALSE)


## Taxa ----
# read in taxa data
TAXA <- read_xlsx(pathtax, sheet = 1, na = nas)
# remove empty rows, if any
i <- apply(TAXA, 1, function(x) all(is.na(x))) #index of rows with all NAs
TAXA <- TAXA[!i,]
# convert species names to lowercase
TAXA$LBL <- base::tolower(TAXA$LBL)
# save for reuse
write.csv(TAXA,
          file = "./data_benthos/benthos_taxa.csv",
          row.names = FALSE)

# # select names of all polychaete species
# spPOL <- TAXA %>% filter(General_grouping == "polychaete") %>% pull(LBL) %>% unique()
# # select names of all nereididae species
# spNER <- TAXA %>% filter(General_grouping == "polychaete",
#                          Family == "Nereididae") %>% pull(LBL) %>% unique()
# # number of species per group
# length(spPOL)
# length(spNER)


## Biomass ----
### biomass ----

# select and stack MD bm
## select sheet names
sheets_MDBM <- grep("MD.+bm", sheetsfor, value = TRUE)
## load the sheets
MDBMs <- lapply(sheets_MDBM, function(x) read_xlsx(pathfor, sheet = x, na = nas))
## if names do not match across, rename mismatching cases
#lapply(MDBMs, names)
MDBMs[[25]] <- MDBMs[[25]] %>% rename(NET_MESH = NETMESH)
MDBMs[[21]]$SAMPLE_DATE <- gsub("2315", "2015", MDBMs[[21]]$SAMPLE_DATE)
## now stack the data
MDBM <- do.call(rbind, MDBMs)

# select and stack VA bm
## select sheet names
sheets_VABM <- grep("VA.+bm", sheetsfor, value = TRUE)
## load the sheets
VABMs <- lapply(sheets_VABM, function(x) read_xlsx(pathfor, sheet = x, na = nas))
## now stack the data
VABM <- do.call(rbind, VABMs)

# Merge MD and VA biomass data.
# names(MDBM)
# names(VABM)
# rename MD columns to match VA for merging
MDBM <- MDBM %>%
    rename(Stratum = STAEQ89)
VABM <- VABM %>%
    rename(Stratum = STRATUM)
# merge MD and VA using matching columns
matchingcols <- base::intersect(names(MDBM), names(VABM))
BM <- rbind(MDBM[,matchingcols], VABM[,matchingcols]) %>%
    mutate(SAMPLE_DATE = as.Date(SAMPLE_DATE))
# convert species names to lowercase
BM$LBL <- base::tolower(BM$LBL)
# add the year column
BM$Year <- as.numeric(format(BM$SAMPLE_DATE, "%Y"))

# check if all species in BM are in the species (TAXA) table.
all(is.element(BM$LBL, TAXA$LBL))

# stations with missing stratum info
StMissingStratum <- unique(with(BM, STATION[is.na(Stratum)]))
StMissingStratum
## these stations correspond to only once sampled stratum (see the next section on ev)
## so remove them from BM
BM <- BM %>% filter(!is.na(Stratum))

# stations not represented in EV, remove from BM
StMissingEV <- c("2171", "2174", "2179", "2185", "2048")
filter(BM, STATION %in% StMissingEV)
BM <- BM %>% filter(STATION %notin% StMissingEV)



### sampling events ----

# select and stack MD ev
## select sheet names
sheets_MDEV <- grep("MD.+ev", sheetsfor, value = TRUE)
## load the sheets
MDEVs <- lapply(sheets_MDEV, function(x) read_xlsx(pathfor, sheet = x, na = nas))
## if names do not match across, rename mismatching cases
#lapply(MDEVs, names)
MDEVs[[25]] <- MDEVs[[25]] %>% rename(SAMPTYPE = SAMP_TYPE)
## now stack the data
MDEV <- do.call(rbind, MDEVs)

# select and stack VA ev
## select sheet names
sheets_VAEV <- grep("VA.+ev", sheetsfor, value = TRUE)
## load the sheets
VAEVs <- lapply(sheets_VAEV, function(x) read_xlsx(pathfor, sheet = x, na = nas))
## now stack the data
VAEV <- do.call(rbind, VAEVs)

# merge MD and VA event data.
# names(MDEV)
# names(VAEV)
# rename MD columns to match VA for merging
# note that MDEV already has STRATUM (unlike bm sheets) that is the same as STAEQ89
all(MDEV$STRATUM == MDEV$STAEQ89)
MDEV <- MDEV %>% rename(SAMP_TYPE = SAMPTYPE)
# merge MD and VA using matching columns
# matchingcols <- base::intersect(names(MDEV), names(VAEV))
# actually, keep only columns that are needed
matchingcols <- c("STATION", "SAMPLE_DATE", "STRATUM",
                  "SITE_TYPE", #"SAMP_TYPE",
                  "LATITUDE", "LONGITUDE",
                  "TOTAL_DEPTH")
EV <- rbind(MDEV[,matchingcols], VAEV[,matchingcols]) %>%
    mutate(SAMPLE_DATE = as.Date(SAMPLE_DATE)) %>%
    rename(Stratum = STRATUM) %>%
    mutate(Year = as.numeric(format(SAMPLE_DATE, "%Y"))) #add the variable: year

# check the uniqueness of random stations and adjust EV, if needed
EVrand <- EV %>%
    filter(SITE_TYPE == "RANDOM") #select random sites in MD and VA
length(unique(EVrand$STATION)) == nrow(EVrand)
# find duplicated stations
i <- EVrand$STATION[duplicated(EVrand$STATION)]
# find it in ev
filter(EV, STATION == i)
# find it in bm
filter(BM, STATION == i)
# keep only one in ev, with the date as in bm
tmp <- filter(EV, STATION == i) #extract rows with the duplicated station
tmp[2, 4:ncol(tmp)] <- tmp[1, 4:ncol(tmp)] #copy part of the info from 1st to 2nd row
#replace a row in the data with complete info, remove the other (duplicated) rows
EV <- EV %>%
    filter(STATION != i) %>%
    rbind(tmp[2,])

# check stratum of StMissingStratum
with(EV, Stratum[STATION %in% StMissingStratum])
str2remove <- c("07T")
# remove data from str2remove
EV <- EV %>%
    filter(Stratum %notin% str2remove)


### combine ----
BMEV <- full_join(BM, EV, by = c("Year", "Stratum", "STATION", "SAMPLE_DATE"))
# summary(BMEV)

# there were some stations in MD 1995 that were not recorded in the event data,
# remove them
tmp <- BMEV %>%
    filter(is.na(LATITUDE))
table(tmp$Year)
BMEV <- BMEV %>%
    filter(!is.na(LATITUDE))
summary(BMEV)
write.csv(BMEV,
          file = "./data_benthos/benthos_biomass.csv",
          row.names = FALSE)
