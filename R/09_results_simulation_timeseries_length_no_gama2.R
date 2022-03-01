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
# saveRDS(summed_result, "Results/Time Series Length/tsl results February 21st")

# Plotting:

?mean_cl_normal

summed_result <- readRDS("Results/Time Series Length/tsl results February 21st")

summed_result[["coverage"]] %>%
  ggplot()+
  aes(x = length, y = (beta1_fn/1000), group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_cl_normal, geom = "line")+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  geom_hline(yintercept = 0.95, color = "red", linetype = 2)+
  facet_wrap(~beta1)+
  labs(x = "Time Series Length", y = "Proportion of CI coverage")+
  theme_demo+
  theme(legend.position = "none")

summed_result[["ci_width"]] %>%
  ggplot()+
  aes(x = length, y = (beta1_fn), group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_cl_normal, geom = "line")+
  #scale_y_log10()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1)+
  coord_cartesian(ylim = c(0,0.01))+
  labs(x = "Length", y = "CI Width")+
  theme_demo+
  theme(legend.position = "none")

ggplot()+
  aes(x = length, group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  # stat_summary(data = summed_result[["coverage"]], aes(y = beta1_fn), fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  # stat_summary(data = summed_result[["coverage"]], aes(y = beta1_fn), fun.data = mean_se, geom = "line")+
  stat_summary(data = summed_result[["ci_width"]], aes(y = beta1_fn), fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  stat_summary(data = summed_result[["ci_width"]], aes(y = beta1_fn), fun.data = mean_se, geom = "line")+
  scale_y_continuous(sec.axis = sec_axis(trans = function(x) sqrt(x), name = "CI Width"))+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1)

summed_result[["coverage"]] %>%
  ggplot()+
  aes(x = length, y = as.factor(beta1), fill = (beta1_fn/1000) >= 0.95)+
  geom_tile()+
  facet_wrap(~beta0)


summed_result[["mse"]] %>%
  ggplot()+
  aes(x = length, y = beta1_fn, group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_se, geom = "line")+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  scale_y_log10()+
  geom_line(mapping = aes(y = beta1), color = "red", linetype = 2)+
  facet_wrap(~beta1)


summed_result[["bias"]] %>%
  # mutate(bias_percentage = beta1_)
  ggplot()+
  aes(x = length, y = 100*beta1_fn/beta1, group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_cl_normal, geom = "line")+
  coord_cartesian(ylim = c(-100, 1))+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1)+
  labs(x = "Length", y = "Bias %")+
  theme_demo+
  theme(legend.position = "none")

summed_result[["variance"]] %>%
  ggplot()+
  aes(x = length, y = (beta1_fn), group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2)+
  stat_summary(fun.data = mean_se, geom = "line")+
  coord_cartesian(ylim = c(0, 0.01))+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1)+
  labs(x = "Length", y = "Variance")+
  theme(legend.position = "none")


summed_result[["bias"]] %>%
  ggplot()+
  aes(x = as.factor(length), y = as.factor(beta1), fill = beta1_fn >= 0)+
  geom_tile()+
  coord_equal()+
  facet_wrap(~beta0)


left_join(summed_result[["ci_width"]], summed_result[["coverage"]], by = c('beta0', 'beta1', 'length')) %>% 
  filter(beta1_fn.x < 1) %>% 
  ggplot()+
  aes(x = beta1_fn.x, y = beta1_fn.y/1000, group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  geom_point()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  scale_x_log10()+
  facet_wrap(~beta1)+
  labs(x = "CI Width", y = "Proportion of CI coverage")+
  theme_demo+
  theme(legend.position = "none")

left_join(summed_result[["bias"]], summed_result[["ci_width"]], by = c('beta0', 'beta1', 'length')) %>%
  filter(beta1_fn.y < 1) %>% 
  ggplot()+
  aes(y = abs(beta1_fn.x), x = beta1_fn.y, group = as.factor(beta1), color = as.factor(beta1), fill = as.factor(beta1))+
  geom_point()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  facet_wrap(~beta1)+
  scale_x_log10()+scale_y_log10()+
  labs(y = "Bias", x = "CI Width")+
  theme(legend.position = "none")
