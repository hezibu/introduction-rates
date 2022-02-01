source("R/00_required_packages.R")
# Rcpp::sourceCpp("R/01_solow_costello_cpp_functions.cpp")
source("R/02_global_functions.R")

sim_params_ts_length <- set_params(c(-1.1, 0.014, -1.46, 0.00001, 0.0000004))


# tsl_results_parts_list <- readRDS("Results/Prior Invasion/prior invasion simulation results 31012022")
# 
# tsl_results_parts_list <- unlist(tsl_results_parts_list, recursive = FALSE,use.names = TRUE)
# 
# saveRDS(tsl_results_parts_list, "Results/Prior Invasion/prior invasion simulation unlisted results 31012022")

tsl_results <- readRDS("Full Result Lists/Prior Invasion/prior invasion simulation unlisted results 31012022")

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

dir.create("Results/Prior Introduction")

rm(tsl_results)

summed_result <- list(bias = biases, variance = variances, ci_width = ci_widths, mse = mses, coverage = coverages)


dir.create("Results/Prior Introduction", recursive = TRUE)
saveRDS(summed_result, "Results/Prior Introduction/tsl results February 1st")