source("R/00_required_packages.R")
# Rcpp::sourceCpp("R/01_solow_costello_cpp_functions.cpp")
source("R/02_global_functions.R")

sim_params_ts_length <- set_params(c(-1.1, 0.014, -1.46, 0.00001, 0.0000004))

list_tsl <- list(
  beta0 = seq(-1.5,1.5, length.out = 20),
  beta1 = seq(0, 0.03, length.out = 20),
  length = seq(10, 400, by = 15)
)

df_tsl <- expand.grid(list_tsl)

parts_list_tsl <- list(
  `1` = filter(df_tsl, length < 86),
  `2` = filter(df_tsl, length > 86, length < 166),
  `3` = filter(df_tsl, length > 166, length < 244),
  `4` = filter(df_tsl, length > 244, length < 322),
  `5` = filter(df_tsl, length > 322, length <= 400)
)

tsl_results_parts_list <- vector(mode = "list", length = length(parts_list_tsl))


tsl_results_parts_list[[2]] <-  pbapply::pbapply(parts_list_tsl$`2`, MARGIN = 1, function(row){
  
  sim_params_ts_length[1:2] <- as.numeric(row[1:2]) 
  
  length <- as.numeric(row[3])
  
  result <- replicate(1000, expr = {
    
    simulation <- sim(length, params = sim_params_ts_length)
    
    years <- seq_along(simulation)-1
    
    # We'll use the coefficients from this model as starting parameters:
    # beta0 = intercept; beta1 = slope
    simple_model <- lm(c(log((simulation) + 1)) ~ years)
    
    guess <- set_params(c(simple_model$coefficients[1],
                          simple_model$coefficients[2],
                          0, 0, 0))[1:4]
    
    estimates <- try(optim(count_log_like_no_gama2, par = guess, constant = 0, method = "BFGS",
                           first_record_data = simulation, hessian = TRUE), silent = T)
    
    return(list(timeseries = simulation, estimates = estimates))
  }, simplify = FALSE)
  
  return(list(sim_params_tsl = c(sim_params_ts_length, length = length), simulation = result))
}, cl = 12)

*dir.create("Full Result List/Time Series Length/", recursive = TRUE)

saveRDS(tsl_results_parts_list,"Full Result List/Time Series Length/time series length simulation results 03022022")
saveRDS(parts_list_tsl,"Full Result List/Time Series Length/time series length parameters 03022022")
