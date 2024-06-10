# Hypoxia network for each 1 year -- time series analysis of network stats.

# Packages ----

rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)

library(ggplot2)
theme_set(theme_light())
library(GGally)
library(patchwork)

# Nets loop ----
# Load data to create networks and extract their statistics, for each year.

if (FALSE) {
    library(igraph)
    source("code/fun_create_adj.R")
    source("code/fun_netstats.R")

    rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
        filter(FSM == 1) %>%
        mutate(CellID = paste0("x", CellID)) %>%
        mutate(Atlantic = (LAT < 37.22 & LON > -75.96) |
                   (LAT < 37.5 & LON > -75.8) |
                   (LAT < 38.0 & LON > -75.6) |
                   (LAT < 36.96 & LON > -75.99)
        ) %>%
        filter(!Atlantic)

    YEARS = 1986L:2015L
    alpha = 0.05

    NS <- as.list(rep(NA, length(YEARS)))
    names(NS) <- YEARS
    i = 1
    for (year in YEARS) { # year = 2002
        print(year)

        # Load data obtained on cluster
        Mcoef <- readRDS(paste0("dataderived/CAUSpairs_", year, ".rds"))
        proc <- create_adj_pairs(Mcoef)

        # Thresholding / adjusting
        # proc$Acoef[proc$Acoef < 0.01] <- 0
        proc$Gp <- apply(proc$Gp, 2, p.adjust, method = "BY")
        diag(proc$Gp) <- 1

        # Unweighted directed adjacency matrix
        AM <- proc$Gp < alpha

        # Weighted directed adjacency matrix
        AMw <- proc$Acoef
        AMw[AM == 0] <- 0

        # Create graph, add attributes
        graphGranger <- igraph::graph_from_adjacency_matrix(AMw, weighted = TRUE
                                                            #AM, weighted = NULL
                                                            ,mode = "directed"
                                                            ,diag = FALSE)
        # igraph::is_weighted(graphGranger)

        # Stats
        NS[[i]] <- netstats(graphGranger)
        i <- i + 1
    }
    saveRDS(NS, "dataderived/NetStatsTS.rds")
}


# Analysis ----

YEARS = 1986L:2015L

## Network statistics ----
NS <- readRDS("dataderived/NetStatsTS.rds")

# Transform the list into a tibble with data
N <- tibble(Year = YEARS,
            indegree = sapply(NS, function(x) mean(x$indegree)),
            outdegree = sapply(NS, function(x) mean(x$outdegree)),
            closeness = sapply(NS, function(x) mean(x$closeness)),
            betweenness = sapply(NS, function(x) mean(x$betweenness)),
            # triads = ,
            # components = ,
            G_order = sapply(NS, function(x) x$G_order),
            G_size = sapply(NS, function(x) x$G_size),
            G_density = sapply(NS, function(x) x$G_density),
            G_reciprocity = sapply(NS, function(x) x$G_reciprocity),
            G_centralization = sapply(NS, function(x) x$G_centralization),
            G_pathLen = sapply(NS, function(x) x$G_pathLen),
            G_components = sapply(NS, function(x) x$G_components),
            G_isolated = sapply(NS, function(x) x$G_isolated)
) %>%
    select(-closeness, -betweenness, -G_order)


## Bay stats ----
# Baywide data from Allison Dreiss, email 2023-08-09
HV <- read_csv("data_rca/Hypoxic_volume_all_data_corr.csv") %>%
    rename(HVolume_km3 = `Volume_ km^3`)
HA <- read_csv("data_rca/Hypoxic_area_all_data_corr.csv") %>%
    rename(HArea_km2 = `Area in km^2`)

### whole year
hv_y <- HV %>%
    group_by(Year) %>%
    summarise(HV_year_avg = mean(HVolume_km3),
              HV_year_max = max(HVolume_km3),
              HV_year_med = median(HVolume_km3))

ha_y <- HA %>%
    group_by(Year) %>%
    summarise(HA_year_avg = mean(HArea_km2, na.rm = TRUE),
              HA_year_max = max(HArea_km2, na.rm = TRUE),
              HA_year_med = median(HArea_km2, na.rm = TRUE))

### summer months
season <- 5L:8L
hv_s <- HV %>%
    filter(Month %in% season) %>%
    group_by(Year) %>%
    summarise(HV_s_avg = mean(HVolume_km3),
              HV_s_max = max(HVolume_km3),
              HV_s_med = median(HVolume_km3))

ha_s <- HA %>%
    filter(Month %in% season) %>%
    group_by(Year) %>%
    summarise(HA_s_avg = mean(HArea_km2, na.rm = TRUE),
              HA_s_max = max(HArea_km2, na.rm = TRUE),
              HA_s_med = median(HArea_km2, na.rm = TRUE))

# Combine annual summaries
hlist <- list(hv_y, hv_s, ha_y, ha_s)
H <- hlist %>%
    purrr::reduce(full_join, by = "Year") %>%
    filter(Year %in% YEARS)

# Combine with network stats
D <- full_join(N, H, by = "Year")
Dlong <- D %>%
    pivot_longer(cols = 2:ncol(D), names_to = "Variable", values_to = "Value")

## Plots ----

# jpeg("images/hynet1_ts.jpeg", width = 11, height = 7, units = "in", res = 300)
png("images/hynet1_ts.png", width = 11, height = 7, units = "in", res = 300)
Dlong %>%
    filter(!grepl("HA_year_med", Variable)) %>%
    filter(!grepl("_max", Variable)) %>%
    ggplot(aes(x = Year, y = Value)) +
    geom_line() +
    facet_wrap(vars(Variable), ncol = 4, scales = "free_y")
dev.off()


png("images/hynet1_mat.png", width = 15, height = 15, units = "in", res = 300)
D %>%
    select(-Year) %>%
    select(-ends_with("_max")) %>%
    select(-G_size) %>%
    select(contains(c("G_", "HA_year_avg", "HV_year_avg"))) %>%
    GGally::ggpairs()
dev.off()

library(ggcorrplot)
corr <- D %>%
    select(-Year) %>%
    select(-ends_with("_max")) %>%
    select(-G_size) %>%
    select(contains(c("G_", "HA_year_avg", "HV_year_avg"))) %>%
    as.matrix() %>%
    cor() %>%
    round(digits = 2)
p.mat <- D %>%
    select(-Year) %>%
    select(-ends_with("_max")) %>%
    select(-G_size) %>%
    select(contains(c("G_", "HA_year_avg", "HV_year_avg"))) %>%
    as.matrix() %>%
    cor_pmat()
p1 <- ggcorrplot::ggcorrplot(
    corr,
    p.mat = p.mat,
    hc.order = FALSE,
    lab = TRUE,
    type = "lower",
    insig = "pch",
    colors = c("#6D9EC1", "white", "#E46726")
) +
    ggtitle("") + theme(plot.margin = unit(c(0,0,0,0),"mm"))

png("images/hynet1_cormat.png", width = 7, height = 5, units = "in", res = 300)
p1
dev.off()

apply(D, 2, mean)
# Year         indegree        outdegree           G_size        G_density    G_reciprocity G_centralization        G_pathLen     G_components
# 2.000500e+03     3.235483e+02     3.235483e+02     9.913521e+05     1.056312e-01     1.489208e-01     4.340684e-03     5.728271e-02     3.806667e+01
# G_isolated      HV_year_avg      HV_year_max      HV_year_med         HV_s_avg         HV_s_max         HV_s_med      HA_year_avg      HA_year_max
# 3.380000e+01     1.328267e+01     5.176651e+01     3.636431e+00     2.862888e+01     5.087745e+01     2.999562e+01     5.346682e+02     2.657092e+03
# HA_year_med         HA_s_avg         HA_s_max         HA_s_med
# 7.000000e-32     1.282893e+03     2.657092e+03     1.428856e+03
