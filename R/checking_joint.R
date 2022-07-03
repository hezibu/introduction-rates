library(tidyverse)

summarize_stan_trials <- function(a_single_run, parameters, fn){
  real_values_tbl <- bind_rows(a_single_run[[1]])
  relevant_real_values <- select(real_values_tbl, any_of(parameters))
  real_values <- as.numeric(relevant_real_values)
  names(real_values) <- names(relevant_real_values)
  randomizations <- a_single_run[[2]]
  out <- bind_rows(fn(real_values, parameters, randomizations)) %>%
    rename_with(.fn = function(x) str_glue("{x}_fn"))
  return(bind_cols(real_values_tbl, out))
}

report_stan_bias <- function(real_values, parameters, randomizations){
  all_b <- sapply(seq_len(length(randomizations)),
                  FUN = function(i) {
                    pars <- randomizations[[i]]$draws %>% 
                      filter(variable %in% parameters) %>% 
                      pull(mean)
                    out <- pars - real_values # ˆd − d
                    return(out)
                  }
  )
  
  if (is.null(dim(all_b))){
    mean <- mean(all_b) # E[ˆd − d]
    names(mean) <- names(real_values)
    return(mean)
  } else {
    return(rowMeans(all_b)) # E[ˆd − d]
  }
}

report_stan_variance <- function(real_values, parameters, randomizations){
  all_d <- sapply(seq_len(length(randomizations)),
                  FUN = function(i) {
                    pars <- randomizations[[i]]$draws %>% 
                      filter(variable %in% parameters) %>% 
                      pull(mean)
                    return(pars)
                  }
  )
  
  if (is.null(dim(all_d))){
    e_all_d <- mean(all_d) # E[ˆd]
    variance <- (all_d - e_all_d)^2
    mean_variance <- mean(variance)
    names(mean_variance) <- names(real_values)
    return(mean_variance)
    
  } else {
    e_all_d <- rowMeans(all_d) # E[ˆd]
    names(e_all_d) <- names(real_values)
    variance <- apply(all_d, MARGIN = 2, function(set)
      (set - e_all_d)^2) # (ˆd − E[ˆd])^2
    return(rowMeans(variance)) #E[(ˆd − E[ˆd])^2]
  }

  
}

report_stan_mse <- function(real_values, parameters, randomizations){
  bias <- report_stan_bias(real_values, parameters, randomizations)
  variance <- report_stan_variance(real_values, parameters, randomizations)
  return(bias^2 + variance)
}

# test <- readRDS("/data_list/time_series_length_full.RDS")
# saveRDS(full_list, "/results_tables/time_series_length_full_beta0_0.RDS")


test_bias <- pbapply::pblapply(full_list, function(trial) {
  return(try(summarize_stan_trials(trial, c("beta1","b1"), report_stan_bias)))
}) %>% bind_rows()


test_variance <- pbapply::pblapply(full_list, function(trial) {
  return(try(summarize_stan_trials(trial, c("beta1","b1"), report_stan_variance)))
}) %>% bind_rows()

test_mse <- pbapply::pblapply(full_list, function(trial) {
  return(try(summarize_stan_trials(trial, c("beta1","b1"), report_stan_mse)))
}) %>% bind_rows()


theme_demo <- theme_classic()+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))

ggplot(test_bias)+
  aes(x = length, y = 100*beta1_fn/beta1, color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.5) + 
  stat_summary(fun = mean, geom = "point")+
  stat_summary(fun = mean, geom = "line")+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  #facet_wrap(~beta1)+
  labs(x = "Time Series Length", y = "Bias %")+
  theme_demo+ # facet_wrap(~beta0)+
  theme(legend.position = c(0.9, 0.15))

ggplot(test_bias)+
  aes(x = length, y = 100*beta1_fn/beta1, color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.5) + 
  stat_summary(fun = mean, geom = "point")+
  stat_summary(fun = mean, geom = "line")+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  #facet_wrap(~beta1)+
  labs(x = "Time Series Length", y = "Bias %")+
  theme_demo+ # facet_wrap(~beta0)+
  theme(legend.position = c(0.9, 0.15))


ggplot(test_variance)+
  aes(x = length, y = beta1_fn, color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.5) + 
  stat_summary(fun = mean, geom = "point")+
  stat_summary(fun = mean, geom = "line")+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  #facet_wrap(~beta1)+
  labs(x = "Time Series Length", y = "Variance")+
  theme_demo+
  # facet_wrap(~beta0)+
  theme(legend.position = c(0.1, 0.8))


ggplot(test_mse)+
  aes(x = length, y = beta1_fn, color = as.factor(beta1), fill = as.factor(beta1))+
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.5) + 
  stat_summary(fun = mean, geom = "point")+
  stat_summary(fun = mean, geom = "line")+
  geom_hline(yintercept = 0, color = "red", linetype = 2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  #facet_wrap(~beta1)+
  labs(x = "Time Series Length", y = "MSE")+
  theme_demo+
  # facet_wrap(~beta0)+
  theme(legend.position = c(0.1, 0.8))

what <- readRDS("/results_tables/time_series_length_part_1_of_20.RDS")

what_2 <- readRDS("/results_tables/time_series_length_22062022.RDS")

saveRDS(list(test_bias, test_variance, test_mse), "/results_tables/time_series_length_22062022.RDS")
