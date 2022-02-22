source("R/00_required_packages.R")
# Rcpp::sourceCpp("R/01_solow_costello_cpp_functions.cpp")
source("R/02_global_functions.R")

sim_params_ts_length <- set_params(c(-1.1, 0.014, -1.46, 0.00001, 0.0000004))


# full_results_list <- lapply(paste0("Full Result Lists/Time Series Length/",list.files("Full Result Lists/Time Series Length/")[-1]),
#                             readRDS)
# 
# full_results_list <- unlist(full_results_list, recursive = FALSE)
# 
# saveRDS(full_results_list,"Full Result Lists/Time Series Length/time series length simulation unlisted results")

tsl_results <- readRDS("Full Result Lists/Time Series Length/time series length simulation unlisted results")

biases <- pblapply(tsl_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_bias)))
})

biases <- biases[sapply(biases, function(x) all(class(x) != "try-error"))] %>% bind_rows()

variances <- pblapply(tsl_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_variance)))
})

variances <- variances[sapply(variances, function(x) all(class(x) != "try-error"))] %>% bind_rows()

ci_widths <-  pblapply(tsl_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_ci_width)))
})

ci_widths <- ci_widths[sapply(ci_widths, function(x) all(class(x) != "try-error"))] %>% bind_rows()

mses <-  pblapply(tsl_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_mse)))
})

mses <- mses[sapply(mses, function(x) all(class(x) != "try-error"))] %>% bind_rows()

coverages <-  pblapply(tsl_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_coverage)))
})

coverages <- coverages[sapply(coverages, function(x) all(class(x) != "try-error"))] %>% bind_rows()

rm(tsl_results)

summed_result <- list(bias = biases,
                      variance = variances,
                      ci_width = ci_widths, 
                      mse = mses, 
                      coverage = coverages)


# dir.create("Results/Time Series Length", recursive = TRUE)
saveRDS(summed_result, "Results/Time Series Length/tsl results February 21st")
