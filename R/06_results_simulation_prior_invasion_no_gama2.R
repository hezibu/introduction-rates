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


# Plotting:

summed_result <- readRDS("Results/Prior Introduction/tsl results February 1st")

summed_result[["coverage"]] %>%
  ggplot()+
  aes(x = as.factor(actual_length), y = as.factor(beta1), fill = beta1_fn >= 950)+
  geom_tile()+
  coord_equal()+
  facet_wrap(~beta0)

summed_result[["coverage"]] %>%
  ggplot()+
  aes(x = actual_length, y = (beta1_fn/1000), group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_se, geom = "line")+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  geom_hline(yintercept = 0.95, color = "red", linetype = 2)+
  facet_wrap(~beta1)+
  labs(x = "Actual introduction Length", y = "Proportion of CI coverage")+
  theme(legend.position = "none")

summed_result[["ci_width"]] %>%
  ggplot()+
  aes(x = actual_length, y = (beta1_fn), group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_se, geom = "line")+
  #scale_y_log10()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1,scales = "free_y")+
  labs(x = "Actual introduction Length", y = "CI Width")+
  theme(legend.position = "none")

summed_result[["bias"]] %>%
  ggplot()+
  aes(x = actual_length, y = (beta1_fn), group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_se, geom = "line")+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1, scales = "free_y")+
  labs(x = "Actual introduction Length", y = "Bias")+
  theme(legend.position = "none")

summed_result[["variance"]] %>%
  ggplot()+
  aes(x = actual_length, y = (beta1_fn), group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_se, geom = "line")+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1, scales = "free_y")+
  labs(x = "Actual introduction Length", y = "Variance")+
  theme(legend.position = "none")


summed_result[["bias"]] %>%
  ggplot()+
  aes(x = as.factor(actual_length), y = as.factor(beta1), fill = beta1_fn >= 0)+
  geom_tile()+
  coord_equal()+
  facet_wrap(~beta0)
