source("R/00_required_packages.R")
# Rcpp::sourceCpp("R/01_solow_costello_cpp_functions.cpp")
source("R/02_global_functions.R")


full_results_list <- readRDS("Full Result Lists/Increasing Rate/simulations 04012022")
            
params_range_results <- unlist(full_results_list, recursive = FALSE)

biases <- pblapply(params_range_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_bias)))
})

biases <- biases[sapply(biases, function(x) all(class(x) != "try-error"))] %>% bind_rows()

variances <- pblapply(params_range_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_variance)))
})

variances <- variances[sapply(variances, function(x) all(class(x) != "try-error"))] %>% bind_rows()

ci_widths <-  pblapply(params_range_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_ci_width)))
})

ci_widths <- ci_widths[sapply(ci_widths, function(x) all(class(x) != "try-error"))] %>% bind_rows()

mses <-  pblapply(params_range_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_mse)))
})

mses <- mses[sapply(mses, function(x) all(class(x) != "try-error"))] %>% bind_rows()

coverages <-  pblapply(params_range_results, function(trial) {
  return(try(summarize_trials(trial, 4, report_coverage)))
})

coverages <- coverages[sapply(coverages, function(x) all(class(x) != "try-error"))] %>% bind_rows()

rm(params_range_results)

summed_result <- list(bias = biases,
                      variance = variances,
                      ci_width = ci_widths, 
                      mse = mses, 
                      coverage = coverages)

