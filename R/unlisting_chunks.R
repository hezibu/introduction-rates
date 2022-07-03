library(future.apply)

n_files <- length(list.files("/data_list/", full.names = T))

plan("multicore", workers = n_files)

full_list <- future_lapply(list.files("/data_list/", full.names = T), readRDS)

plan("sequential")

cli::cli_alert_success("Successfully read files")

full_list <- unlist(full_list, recursive = FALSE)


cli::cli_alert_success("Successfully unlisted list")

saveRDS(full_list,"/data_list/time_series_length_full.RDS")

cli::cli_alert_success("Done.")