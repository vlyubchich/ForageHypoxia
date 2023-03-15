## Primary working directory ----
pwd <- "/Users/ryan/Windows/Documents/Post UCB/Research/MD Sea Grant/Grant_2/ForageHypoxia"
# pwd <- "/home/ryan/ForageHypoxia"

## Primitives ----
# cat("\014")
# rm(list = ls())
# gc()
set.seed(42)

## Libraries ----
library("readr")
library("dplyr")
library("magrittr")
library("nonlinearTseries")
library("keras")
library("tensorflow")
library("ggplot2")

# ## Parallelization
# library("foreach")
# library("parallel")
# library("doParallel")

## Variables
num_lags <- 6 ## Will include one extra time step for the present
lag_gap <- 10
log_biomass <- TRUE
dither <- FALSE
train_percent <- 0.7

## UUID (1 second resolution) for each model run so we do not accidentally save over previous model fit output
tu_id <- gsub(x = Sys.time(), pattern = "-|:| ", replacement = "")

## Data ----
stations_cells <- readr::read_csv(paste0(pwd, "/data_benthos/stations_cells.csv"))
benthos_biomass <- readr::read_csv(paste0(pwd, "/data_benthos/benthos_biomass.csv"))
benthos_strata <- readr::read_csv(paste0(pwd, "/data_benthos/benthos_strata.csv"))
benthos_taxa <- readr::read_csv(paste0(pwd, "/data_benthos/Taxa_group_IDs_updated_rjw_v3_1.csv")) #benthos_taxa.csv"))

## Combine rca files from Slava (one per year) ----
years <- 1986:2015
setwd(paste0(pwd, "/data_rca/rca_ts_1986-2015_2"))
for (y in years) {
    print(y)

    rca_ts_addition <- readr::read_csv(paste0("rca_ts_", y, "_2.csv"))
    if (y == min(years)) {
        rca_ts <- rca_ts_addition
    } else {
        rca_ts %<>% dplyr::bind_rows(rca_ts_addition)
    }
}
setwd(pwd)

## Save combined rca data ----
## Commented out because only need to save once
setwd(paste0(pwd, "/data_rca"))
write.table(
    x = rca_ts,
    file = "rca_ts_1986-2015_2.csv",
    sep = ",",
    col.names = TRUE,
    row.names = FALSE
)
setwd(pwd)


## Model data prep ----
benthos_biomass_grp <- dplyr::left_join(benthos_biomass, benthos_taxa, by = "LBL") %>%
    dplyr::group_by(STATION, SAMPLE_DATE, Aggregate_grp) %>%
    dplyr::summarize(VALUE = sum(VALUE, na.rm = TRUE)) %>%
    ungroup()

benthos_biomass_grp %<>% dplyr::left_join(dplyr::select(stations_cells, STATION, CellID), by = "STATION")

## Remove NAs and . grouping names (check with Woodland about what the periods mean)
benthos_biomass_grp %<>% dplyr::filter(!is.na(Aggregate_grp), Aggregate_grp != ".")

## Cast data across station, cell ID, and date
benthos_biomass_grp_cast <- benthos_biomass_grp %>% 
    reshape2::dcast(STATION + CellID + SAMPLE_DATE ~ Aggregate_grp, value.var = "VALUE") %>%
    as_tibble()

## Limit benthos biomass data to rca environmental dates (stops in 2015)
benthos_biomass_grp_cast %<>% dplyr::filter(SAMPLE_DATE <= "2015-12-31")

## rca_ts is enormous, so help by removing rows with CellID not in the benthos data
rca_ts_inputs <- rca_ts %>% dplyr::filter(CellID %in% unique(benthos_biomass_grp_cast$CellID))

## Save combined rca data subsetted to those with CellID in the benthos data (smaller file size) ----
## Commented out because only need to save once
setwd(paste0(pwd, "/data_rca"))
write.table(
    x = rca_ts_inputs,
    file = "rca_ts_1986-2015_2_biomassCellID.csv",
    sep = ",",
    col.names = TRUE,
    row.names = FALSE
)
setwd(pwd)

## Change Date column name to SAMPLE_DATE so can be combined with the benthos data
colnames(rca_ts_inputs)[which(colnames(rca_ts_inputs) == "Date")] = "SAMPLE_DATE"

## Pick a benthos taxa/functional group before adding in, and lagging, the environmental variables in rca_ts
benthos_names <- benthos_biomass_grp$Aggregate_grp %>% unique() %>% sort()

## Loop across benthic functional groups ----
for (name in benthos_names) {
    biomass_data <- benthos_biomass_grp_cast %>% dplyr::select(all_of(c("STATION", "CellID", "SAMPLE_DATE", name)))

    ## Remove NAs
    biomass_data %<>% dplyr::filter(!is.na(get(name)), !is.na(CellID))

    ## log transform (specified at the top)
    if (log_biomass) {
        biomass_data[[name]] = log(1 + biomass_data[[name]])
    }

    ## Only consider benthos with at least 20 observations
    if (nrow(biomass_data) >= 20) {

        ## Loop across biomass_data adding corresponding rca data
        for (row in 1:nrow(biomass_data)) {
            print(row)

            ## Focal observation for each loop
            data_addition <- biomass_data[row,]

            ## Lags, as dates
            lagged_dates <- as.Date(data_addition$SAMPLE_DATE) - c(0, lag_gap*(1:num_lags))

            ## Add rca data at corresponding lagged dates
            rca_row <- rca_ts_inputs %>% 
                dplyr::filter(
                    CellID == data_addition$CellID,
                    SAMPLE_DATE %in% lagged_dates
                ) %>% 
                reshape2::melt() %>%
                as_tibble()

            ## Add lag names (*_t1, *_t2, etc)
            rca_row$variable = paste0(rca_row$variable, rep(paste0("_t", 0:num_lags), ncol(rca_ts_inputs)))

            ## Convert to a single row matrix
            rca_addition <- matrix(rca_row$value, nrow = 1) %>% as_tibble()
            colnames(rca_addition) = rca_row$variable

            ## Combine newly lagged rca data
            if (row == 1) {
                rca_combined <- rca_addition
            } else {
                rca_combined %<>% dplyr::bind_rows(rca_addition)
            }
        }

        ## We no longer need CellID or SAMPLE_DATE columns
        rca_combined %<>% dplyr::select(-which(grepl(pattern = c("CellID|SAMPLE_DATE"), x = colnames(rca_combined))))

        ## Combine benthos and rca data
        data_combined <- dplyr::bind_cols(biomass_data, rca_combined)

        ## Create neural net input ----
        df_input <- data_combined

        ## Save plot of pairwise correlations
        ## Commented out because these only need to be saved once, and can be rather large
        # setwd(paste0(pwd, "/images/correlations"))
        # # grDevices::cairo_pdf(filename = paste0("Pairs_", name, "_", tu_id, ".pdf"), width = 30, height = 30)
        # grDevices::cairo_pdf(filename = paste0("Correlations_", name, "_", tu_id, ".pdf"), width = 30, height = 30)
        #     # fig <- pairs(df_input %>% dplyr::select(-STATION))
        #     fig <- df_input %>% dplyr::select(-c("STATION", "CellID", "SAMPLE_DATE")) %>% reshape2::melt(id.vars = name) %>% as_tibble() %>% ggplot() + theme_classic() + geom_point(aes(x = .data[[name]], y = value)) + facet_wrap(~variable)
        #     print(fig)
        # dev.off()
        # setwd(pwd)

        ## Unit scale data for neural net
        df_input[,sapply(df_input, class) == "numeric"] = sapply(df_input[,sapply(df_input, class) == "numeric"], function(x){x/max(x)}) %>% as_tibble()

        ## Split off time variable from rest of rca data
        data_x <- df_input[, grepl(paste0(c("_avg", "_sd"), collapse = "|"), colnames(df_input))]
        data_time <- df_input[, "SAMPLE_DATE"]

        ## Variable being predicted
        data_y <- df_input[, name]

        ## Dither
        if (dither) {
            data_x_raw <- data_x
            data_x <- ((data_x) + (runif(n = nrow(data_x) * ncol(data_x)) * 1e-5)) %>% as_tibble()
            data_y_raw <- data_y
            data_y <- ((data_y) + (runif(n = nrow(data_y) * ncol(data_y)) * 1e-5)) %>% as_tibble()
        }

        ## Keep data where the prediced variable is non-zero
        ids_keep <- which(unlist(data_y) != 0)
        data_x <- data_x[ids_keep,]
        data_y <- data_y[ids_keep,]
        data_time <- data_time[ids_keep,]

        ## Training/testing split
        train_split <- round(train_percent * nrow(data_y))
        obs <- 1:train_split

        x_train <- data_x[obs,] %>% as.matrix(ncol = ncol(data_x))
        x_test <-  data_x[-obs,] %>% as.matrix(ncol = ncol(data_x))
        y_train <- data_y[obs,] %>% as.matrix(ncol = 1)
        y_test <- data_y[-obs,] %>% as.matrix(ncol = 1)
        time_train <- data_time[obs,] %>% as.matrix(ncol = 1)
        time_test <- data_time[-obs,] %>% as.matrix(ncol = 1)

        ## Neural network ----
        model <- keras_model_sequential() %>%
            layer_dense(
                input_shape = ncol(data_x),
                units = 10,
                activation = "relu"
            ) %>%
            layer_dense(
                units = 20,
                activation = "relu"
            ) %>%
            layer_dense(
                units = 50,
                activation = "relu"
            ) %>%
            layer_dropout(rate = 0.1) %>%
            layer_dense(
                units = 20,
                activation = "relu"
            ) %>%
            layer_dense(
                units = 5,
                activation = "linear"
            )

        ## Compile the model
        model %>% compile(
            optimizer = optimizer_adam(lr = 1e-5, decay = 1e-7),
            loss = "mse", #'binary_crossentropy',
            optimizer = "adam"
        )
        
        ## Fit the model
        model %>% fit(
            verbose = 0,
            x = x_train,
            y = y_train,
            validation_split = 0.2, ## FROM SLAVA
            epochs = 100
        ) 

        ## Predict testing data
        y_pred <- model %>% predict(x_test)

        ## Assemble model output for plotting
        ggdata <- tibble(
            Time = time_test[,1],
            ID = 1:nrow(time_test),
            Observed = y_test[,1],
            Predicted = y_pred[,1]
        )
        ggdata$Time = as.Date(ggdata$Time)

        ## Un-log the data (if logged above)
        if (log_biomass) {
            ggdata$Observed = exp(ggdata$Observed) - 1
            ggdata$Predicted = exp(ggdata$Predicted) - 1
        }

        ## Save (rounded) loss to include in plot
        loss <- model %>% 
            evaluate(
                x_test,
                y_test,
                verbose = 0
            ) %>%
            as.numeric() %>%
            round(2)

        ## Plot model fits ----
        ## Observed vs predicted
        setwd(paste0(pwd, "/images/neuralnet_predictions"))
            grDevices::cairo_pdf(filename = paste0("ObservedVsPredicted_", name, "_", tu_id, ".pdf"), width = 10, height = 10)
                fig <- ggdata %>%
                    ggplot() +
                        cowplot::theme_minimal_grid() +
                        geom_point(
                            aes(
                                x= Observed,
                                y = Predicted
                            ),
                            size = 10
                        ) +
                        geom_point(
                            aes(
                                x= Observed,
                                y = Predicted
                            ),
                            size = 9,
                            color = "white"
                        ) +
                        geom_point(
                            aes(
                                x= Observed,
                                y = Predicted
                            ),
                            size = 9,
                            alpha = 0.1
                        ) +
                        geom_abline(
                            slope = 1,
                            intercept = 0,
                            color = "slateblue"
                        ) +
                        labs(
                            title = paste0(name , " | Correlation = ", cor(ggdata$Observed, ggdata$Predicted) %>% round(2), " | Loss = ", loss)
                        )
                print(fig)
            dev.off()

            ## Observed vs predicted time series
            grDevices::cairo_pdf(filename = paste0("TS_ObservedAndPredicted_", name, "_", tu_id, ".pdf"), width = 10, height = 5)   
                fig <- ggdata %>% 
                    reshape2::melt(
                        id.vars = c("Time", "ID")
                    ) %>%
                    as_tibble() %>%
                    ggplot() +
                        cowplot::theme_minimal_grid() +
                        geom_line(
                            aes(
                                x= ID,
                                y = value,
                                color = variable
                            )
                        ) +
                        scale_color_brewer(palette = "Set1") +
                        labs(
                            y = "",
                            title = paste0(name , " | Correlation = ", cor(ggdata$Observed, ggdata$Predicted) %>% round(2), " | Loss = ", loss)
                        )
                print(fig)
            dev.off()
        setwd(pwd)

    } ## if (nrow(biomass_data) >= 20) {
} ## for (name in benthos_names) {
