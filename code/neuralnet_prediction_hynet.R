## Primitives ----
cat("\014")
rm(list = ls())
gc()
set.seed(42)

## Primary working directory ----
# pwd <- "/Users/ryan/Windows/Documents/Post UCB/Research/MD Sea Grant/Grant_2/ForageHypoxia"
pwd <- "/home/ryan/ForageHypoxia"

## Libraries ----
library("readr")
library("dplyr")
library("reshape2")
library("magrittr")
library("nonlinearTseries")
library("keras")
library("tensorflow")
library("ggplot2")
library("gridExtra")
library("randomForest")
library("cowplot")
library("rEDM")

## Parallelization
library("foreach")
library("parallel")
library("doParallel")

## Variables
num_lags <- 1 ## Will include one extra time step for the present
log_biomass <- TRUE
dither <- FALSE
dither_amount <- 1e-5
train_percent <- 0.6
remove_zeros <- FALSE ## FALSE to save all data for Slava, TRUE for model fitting

## UUID (1 second resolution) for each model run so we do not accidentally save over previous model fit output
tu_id <- gsub(x = Sys.time(), pattern = "-|:|[.]| ", replacement = "")

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
# setwd(paste0(pwd, "/data_rca"))
# write.table(
#     x = rca_ts,
#     file = "rca_ts_1986-2015_2.csv",
#     sep = ",",
#     col.names = TRUE,
#     row.names = FALSE
# )
# setwd(pwd)


## Model data prep ----
benthos_biomass_grp <- dplyr::left_join(benthos_biomass, benthos_taxa, by = "LBL") %>%
    dplyr::group_by(STATION, SAMPLE_DATE, Aggregate_grp) %>%
    dplyr::summarize(VALUE = sum(VALUE, na.rm = TRUE)) %>%
    ungroup()

benthos_biomass_grp %<>% dplyr::left_join(dplyr::select(stations_cells, STATION, CellID), by = "STATION")

## Remove NAs and . grouping names (check with Woodland about what the periods mean)
benthos_biomass_grp %<>% dplyr::filter(!is.na(Aggregate_grp), Aggregate_grp != ".")

## Remove trouble CellID but NOTE: Slava thought was in the Atlantic (see following lines for his cutoffs) but it does not seem to be
# ## Remove Atlantic cells in the way Slava did
# Atlantic =
#     (LATITUDE < 37.22 & LONGITUDE > -75.96) |
#     (LATITUDE < 37.5 & LONGITUDE > -75.8) |
#     (LATITUDE < 38.0 & LONGITUDE > -75.6) |
#     (LATITUDE < 36.96 & LONGITUDE > -75.99)

benthos_biomass_grp %<>% dplyr::filter(!(CellID %in% c(1582)))

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
# setwd(paste0(pwd, "/data_rca"))
# write.table(
#     x = rca_ts_inputs,
#     file = "rca_ts_1986-2015_2_biomassCellID.csv",
#     sep = ",",
#     col.names = TRUE,
#     row.names = FALSE
# )
# setwd(pwd)

## Change Date column name to SAMPLE_DATE so can be combined with the benthos data
colnames(rca_ts_inputs)[which(colnames(rca_ts_inputs) == "Date")] = "SAMPLE_DATE"


## Add hynet embeddings
hynet <- readr::read_csv(
    paste0(
        pwd,
        "/data_hynet/",
        "embedding_ASM.csv"
    )
)

## And also add "manual" hynet embeddings
hynet_manual <- readr::read_csv(
    paste0(
        pwd,
        "/data_hynet/",
        "embedding_Manual.csv"
    )
)

## Add year to rca_ts_inputs so can join on it
rca_ts_inputs %<>% dplyr::mutate(Year = substr(rca_ts_inputs$SAMPLE_DATE, 1, 4))

## Remove x prefix to CellID in hynet data so same as in rca_ts_inputs
hynet$CellID = substr(hynet$CellID, 2, nchar(hynet$CellID))
hynet_manual$CellID = substr(hynet_manual$CellID, 2, nchar(hynet_manual$CellID))

## Convert CellID to a string so can join with hynet
rca_ts_inputs$CellID = as.character(rca_ts_inputs$CellID)

## Convert Year to a double so can join with hynet
rca_ts_inputs$Year = as.numeric(rca_ts_inputs$Year)

## Remove Year and Cell because will join on CellID and SAMPLE_DATE
hynet_manual %<>% dplyr::select(-c("Year", "Cell"))

## Join hynet data to rca_ts_inputs
rca_ts_inputs_nohynet <- rca_ts_inputs
rca_ts_inputs %<>% dplyr::left_join(hynet, by = c("CellID", "Year"))
rca_ts_inputs %<>% dplyr::left_join(hynet_manual, by = c("CellID", "SAMPLE_DATE"))

## Remove Year which is no longer needed because SAMPLE_DATE contains the same information
rca_ts_inputs %<>% dplyr::select(-c("Year"))



## Detrend rca_ts_inputs
cols_to_detrend <- which(grepl(pattern="WTEMP|SAL|TPOC|CHLAVEG|DOAVEG|ASM", names(rca_ts_inputs)))
cols_to_detrend_nohynet <- which(grepl(pattern="WTEMP|SAL|TPOC|CHLAVEG|DOAVEG", names(rca_ts_inputs)))


rca_detrend <- rca_ts_inputs[, cols_to_detrend_nohynet]
# rca_ts_inputs_raw <- rca_ts_inputs
rca_ts_inputs[, cols_to_detrend_nohynet] = apply(rca_detrend, 2, \(x){stats::decompose(x = x %>% ts(deltat = 1/365))$random} )



# lag_gap <- sapply(unique(rca_ts_inputs$CellID), \(x){sapply(cols_to_detrend, \(y){z <- rca_ts_inputs %>% dplyr::filter(CellID == x) %>% dplyr::select(y) %>% unlist() %>% as.numeric(); z = z[!is.na(z)]; nonlinearTseries::timeLag(time.series = z, do.plot = FALSE)  })  })

# lag_gap <- apply(rca_ts_inputs[, cols_to_detrend], 2, \(x){nonlinearTseries::timeLag(time.series = x)})

unique_cells <- unique(rca_ts_inputs$CellID)
unique_cells = unique_cells[-which(unique_cells == "1582")] ## Slava wrote on Slack: "The cell ID 1582 matches my checks and is classified as Atlantic cell (and is removed), but the cell has been matched with a station 18M01 in 2011 that has just a bit different coordinates that differ in 3rd and 4th decimals from the "Atlantic" checks so the station is not classified as Atlantic automatically."


## CAREFUL NOT TO USE TOO MANY CORES
cluster <- makeCluster(100, outfile = "")
registerDoParallel(cluster)

cell_loop <- foreach (cell = seq_along(unique_cells), .packages = c("dplyr", "nonlinearTseries")) %dopar% {
    print(cell)

    lag_gap_matrix_addition <- rep(NA, length(cols_to_detrend))
    embed_dim_matrix_addition <- rep(NA, length(cols_to_detrend))

    for (col_id in seq_along(cols_to_detrend)) {
        col = cols_to_detrend[col_id]
        print(c(cell, col))

        z <- rca_ts_inputs %>% dplyr::filter(CellID == unique_cells[cell]) %>% dplyr::select(SAMPLE_DATE, all_of(col))
        z = z[!is.na(z[[2]]),]
        lag_gap_addition = nonlinearTseries::timeLag(time.series = z[[2]], do.plot = FALSE) #, technique = "ami")
        lag_gap_matrix_addition[col_id] = lag_gap_addition

        dim <- nonlinearTseries::estimateEmbeddingDim(
            z[[2]],
            time.lag = lag_gap_addition,
            max.embedding.dim = 20,
            do.plot = FALSE
        )

        embed_dim_matrix_addition[col_id] = dim
    }

    return(list(
        lag_gap_matrix = lag_gap_matrix_addition,
        embed_dim_matrix = embed_dim_matrix_addition
    ))
}


## Shut down parallel backend and free up resources
parallel::stopCluster(cluster)
# closeAllConnections() ## CAREFUL: This will close other jobs
gc()

## Extract output from parallelized loop
lag_gap_matrix <- do.call(rbind, lapply(cell_loop, \(x){x$lag_gap_matrix})) %>% as_tibble()
embed_dim_matrix <- do.call(rbind, lapply(cell_loop, \(x){x$embed_dim_matrix})) %>% as_tibble()


setwd(pwd)


lag_gap_vector = apply(lag_gap_matrix, 2, \(x){mean(x, na.rm = TRUE)}) %>% round()
embed_dim_vector = apply(embed_dim_matrix, 2, \(x){mean(x, na.rm = TRUE)})  %>% round()

names(lag_gap_vector) = colnames(rca_ts_inputs)[cols_to_detrend]
names(embed_dim_vector) = colnames(rca_ts_inputs)[cols_to_detrend]

## Plot number of lags from this approach
lag_gap_vector %>% enframe %>% ggplot() + geom_col(aes(x = name, y = value), fill = "black") + theme_minimal() + labs( x = "", y = "Lag (number of days)")
embed_dim_vector %>% enframe %>% ggplot() + geom_col(aes(x = name, y = value), fill = "black") + theme_minimal() + labs( x = "", y = "Number of Lags to Use") + scale_y_continuous(breaks = c(0:13))
(lag_gap_vector*embed_dim_vector) %>% enframe %>% ggplot() + geom_col(aes(x = name, y = value), fill = "black") + theme_minimal() + labs( x = "", y = "Total Days into the Past Required") + scale_y_continuous(breaks = 365*c(0:8))


## Pick a benthos taxa/functional group before adding in, and lagging, the environmental variables in rca_ts
benthos_names <- benthos_biomass_grp$Aggregate_grp %>% unique() %>% sort()


## Loop across benthic functional groups ----
for (name_id in seq_along(benthos_names)) {
    lag_gap = lag_gap_vector[name_id]
    name = benthos_names[name_id]

    biomass_data <- benthos_biomass_grp_cast %>% dplyr::select(all_of(c("STATION", "CellID", "SAMPLE_DATE", name)))

    ## Remove NAs
    ## NOTE: This is where the different number of observations for each taxa comes from
    biomass_data %<>% dplyr::filter(!is.na(get(name)), !is.na(CellID))


    ## log transform (specified at the top)
    if (log_biomass) {
        # biomass_data[[name]] = log(1 + biomass_data[[name]]) ## Seems ineffective at stretching data with a magnitude less than one
        biomass_data[[name]] = log(biomass_data[[name]])
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
            # rca_row <- rca_ts_inputs %>%
            #     dplyr::filter(
            #         CellID == data_addition$CellID,
            #         SAMPLE_DATE %in% lagged_dates
            #     ) %>%
            #     reshape2::melt(
            #         id.vars = NULL ## Melt all variables
            #     ) %>%
            #     as_tibble()
            rca_row <- rca_ts_inputs[(rca_ts_inputs$CellID == data_addition$CellID) & ((rca_ts_inputs$SAMPLE_DATE == lagged_dates[1]) | (rca_ts_inputs$SAMPLE_DATE == lagged_dates[2])),] #! Only works if num_lags = 1, but MUCH faster than the dplyr::filter approach above
            rca_row %<>% reshape2::melt(
                    id.vars = NULL ## Melt all variables
                ) %>%
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
        rca_combined_celldate <- rca_combined
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

        ## Normalize data for neural net
        cols_to_convert <- which(grepl(pattern="WTEMP|SAL|TPOC|CHLAVEG|DOAVEG|ASM|Emb", names(df_input)))
        df_input[, cols_to_convert] = apply(df_input[, cols_to_convert], 2, as.numeric)
        ## Unit normalization
        # df_input[,sapply(df_input, class) == "numeric"] = sapply(df_input[,sapply(df_input, class) == "numeric"], function(x){x/max(x)}) %>% as_tibble()
        ## Z-score normalization
        df_input[,sapply(df_input, class) == "numeric"] = sapply(df_input[,sapply(df_input, class) == "numeric"], function(x){scale(x, center = TRUE, scale = TRUE)}) %>% as_tibble()

        ## Split off time variable from rest of rca data
        ## Remember that the ASM hynet variables are constant within a year, so no need to include whatever lag was calculated (which seems to be close to a year anyway)
        data_x <- df_input[, grepl(paste0(c("_avg", "_sd", "ASM_t0", "Emb"), collapse = "|"), colnames(df_input))]
        data_time <- df_input[, "SAMPLE_DATE"]

        ## Variable being predicted
        data_y <- df_input[, name]

        ## Dither
        if (dither) {
            data_x_raw <- data_x
            data_x <- ((data_x) + (runif(n = nrow(data_x) * ncol(data_x)) * dither_amount)) %>% as_tibble()
            data_y_raw <- data_y
            data_y <- ((data_y) + (runif(n = nrow(data_y) * ncol(data_y)) * dither_amount)) %>% as_tibble()
        }

        ## Keep data where the predicted variable is non-zero
        if (remove_zeros) {
            ids_keep <- which(unlist(data_y) != 0)
            data_x <- data_x[ids_keep,]
            data_y <- data_y[ids_keep,]
            data_time <- data_time[ids_keep,]
        }

        # ## Save data for Slava
        # setwd(paste0(pwd, "/data_neuralnet"))
        #     write.table(
        #         x = cbind(data_time, tibble(CellID = rca_combined_celldate$CellID_t0), data_y, data_x) %>% as_tibble(),
        #         file = paste0("NeuralNetData_", name,  "_", tu_id, ".csv"),
        #         sep = ",",
        #         col.names = TRUE,
        #         row.names = FALSE
        #     )
        # setwd(pwd)


        ## Training/testing split
        train_split <- round(train_percent * nrow(data_y))

        ## Training and testing data are contiguous
        # obs <- 1:train_split

        ## Training and testing data are randomly chosen
        obs <- round(sample(x = 1:nrow(data_y), size = train_split, replace = FALSE))
        obs = sort(obs)

        x_train <- data_x[obs,] %>% as.matrix(ncol = ncol(data_x))
        x_test <-  data_x[-obs,] %>% as.matrix(ncol = ncol(data_x))

        y_train <- data_y[obs,] %>% as.matrix(ncol = 1)
        y_test <- data_y[-obs,] %>% as.matrix(ncol = 1)

        time_train <- data_time[obs,] %>% as.matrix(ncol = 1)
        time_test <- data_time[-obs,] %>% as.matrix(ncol = 1)

        ## Add a required batch_size
        dim(x_train) <- c(1, dim(x_train))
        dim(x_test) <- c(1, dim(x_test))

        dim(y_train) <- c(1, dim(y_train))
        dim(y_test) <- c(1, dim(y_test))

        ## Neural network
        model <- keras_model_sequential() %>%
            layer_dense(units = 2) %>%
            layer_dense(units = 10, kernel_regularizer = regularizer_l1_l2(l1 = 0.1, l2 = 1)) %>%
            layer_dropout(rate = 0.2) %>%
            layer_dense(units = 100, bias_regularizer = regularizer_l2(0.01)) %>%
            layer_dense(units = 10) %>%
            layer_simple_rnn(units = 100, return_sequences = TRUE) %>%
            layer_simple_rnn(units = 10, return_sequences = TRUE) %>%
            layer_dropout(rate = 0.1) %>%
            layer_simple_rnn(units = 10, return_sequences = TRUE) %>%
            layer_simple_rnn(units = 10, return_sequences = TRUE) %>%
            layer_dropout(rate = 0.1) %>%
            layer_simple_rnn(units = 10, return_sequences = TRUE) %>%
            layer_simple_rnn(units = 100, return_sequences = TRUE) %>%
            layer_dropout(rate = 0.2) %>%
            layer_simple_rnn(units = 100, return_sequences = TRUE) %>%
            layer_dropout(rate = 0.2) %>%
            layer_simple_rnn(units = 100, return_sequences = TRUE) %>%
            layer_dropout(rate = 0.2) %>%
            layer_simple_rnn(units = 10, return_sequences = TRUE) %>%
            layer_simple_rnn(units = 100, return_sequences = TRUE) %>%
            layer_simple_rnn(units = 10, return_sequences = TRUE) %>%
            layer_dropout(rate = 0.1) %>%
            layer_dense(units = 10, bias_regularizer = regularizer_l2(10.01)) %>%
            layer_dropout(rate = 0.2) %>%
            layer_dense(units = 10, bias_regularizer = regularizer_l2(10.01)) %>%
            layer_dense(units = 100, bias_regularizer = regularizer_l2(0.01)) %>%
            layer_dropout(rate = 0.2) %>%
            layer_dense(units = 10) %>%
            layer_dropout(rate = 0.2) %>%
            layer_dense(units = 10) %>%
            layer_dropout(rate = 0.1) %>%
            layer_dense(units = 10, kernel_regularizer = regularizer_l1_l2(l1 = 0.1, l2 = 1)) %>%
            layer_dense(units = 1)

        ## Compile the model
        model %>% compile(
            loss = "mse", #"msle", #"mse", #'binary_crossentropy',
            optimizer = "adam" #optimizer_adam(learning_rate = 1e-5, decay = 1e-7) #"adam"
        )

        ## Fit the model
        model_history <- model %>% fit(
            verbose = 0,
            x = x_train,
            y = y_train,
            # validation_split = 0.2, ## FROM SLAVA
            epochs = 1000
        )

        ## Save convergence plot for the combined multiplanel plot
        fig_model_history <- plot(model_history) + theme_minimal() + labs(title = "Model Training History") + ylim(0, NA)

        ## Predict testing data
        y_pred_train <- model %>% predict(x_train)
        y_pred_test <- model %>% predict(x_test)

        ## Assemble model output for plotting
        ggdata_train <- tibble(
            Time = time_train[,1],
            ID = 1:nrow(time_train),
            Observed = y_train[,1:dim(y_train)[2],], #y_test[,1],
            Predicted = y_pred_train[,1:dim(y_pred_train)[2],] #y_pred[,1]
        )
        ggdata_train$Time = as.Date(ggdata_train$Time)

        ggdata_test <- tibble(
            Time = time_test[,1],
            ID = 1:nrow(time_test),
            Observed = y_test[,1:dim(y_test)[2],], #y_test[,1],
            Predicted = y_pred_test[,1:dim(y_pred_test)[2],] #y_pred[,1]
        )
        ggdata_test$Time = as.Date(ggdata_test$Time)

        ## Save (rounded) loss to include in plot
        loss_train <- model %>%
            evaluate(
                x_train,
                y_train,
                verbose = 0
            ) %>%
            as.numeric() %>%
            round(2)

        loss_test <- model %>%
            evaluate(
                x_test,
                y_test,
                verbose = 0
            ) %>%
            as.numeric() %>%
            round(2)

        ggdata <- dplyr::bind_rows(
            ggdata_train %>% dplyr::mutate(Kind = "Train"),
            ggdata_test %>% dplyr::mutate(Kind = "Test")
        )
        ggdata$Kind = factor(ggdata$Kind, levels = c("Train", "Test"))

        fig_OBSvsPRED_train <- ggdata_train %>%
                            ggplot() +
                                cowplot::theme_minimal_grid() +
                                geom_point(
                                    aes(
                                        x= Observed,
                                        y = Predicted
                                    ),
                                    size = 2 #10
                                ) +
                                geom_abline(
                                    slope = 1,
                                    intercept = 0,
                                    color = "slateblue"
                                ) +
                                labs(
                                    title = paste0("Training | ", name , " | Correlation = ", cor(ggdata_train$Observed, ggdata_train$Predicted) %>% round(2), " | Loss = ", loss_train)
                                )

        fig_OBSvsPRED_test <- ggdata_test %>%
                            ggplot() +
                                cowplot::theme_minimal_grid() +
                                geom_point(
                                    aes(
                                        x= Observed,
                                        y = Predicted
                                    ),
                                    size = 2 #10
                                ) +
                                geom_abline(
                                    slope = 1,
                                    intercept = 0,
                                    color = "slateblue"
                                ) +
                                labs(
                                    title = paste0("Testing | ", name , " | Correlation = ", cor(ggdata_test$Observed, ggdata_test$Predicted) %>% round(2), " | Loss = ", loss_test)
                                )

            ## Plot model fits ----
            setwd(paste0(pwd, "/images/neuralnet_predictions"))
                grDevices::cairo_pdf(filename = paste0("Convergence-TrainingFit-TestingFit_", name, "_", tu_id, ".pdf"), width = 20, height = 10)
                    grid.arrange(fig_model_history, fig_OBSvsPRED_train, fig_OBSvsPRED_test, nrow = 1)
                dev.off()
            setwd(pwd)

    } ## if (nrow(biomass_data) >= 20) {
} ## for (name in benthos_names) {
