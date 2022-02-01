source("R/00_required_packages.R")
# Rcpp::sourceCpp("R/01_solow_costello_cpp_functions.cpp")
source("R/02_global_functions.R")

# simulation: increasing parts:

list_params <- list(
  beta0 = seq(-1.5,1.5, length.out = 10),
  beta1 = seq(0, 0.1, length.out = 10),
  gama0 = seq(-1,1, length.out = 10),
  gama1 = seq(0, 0.01, length.out = 10),
  gama2 = seq(0, 0.0001, length.out = 10)
)


df_params <- expand.grid(list_params)

parts_list <- split(sample(seq_len(nrow(df_params))), 1:10)

parts_list <- lapply(parts_list, function(rows){
  return(df_params[rows,])
})

# We'll create list for partial results:

results_parts_list <- vector(mode = "list", length = length(parts_list))

# Setting up the iteration:

results_parts_list[[1]] <- pbapply::pbapply(parts_list$`1`, MARGIN = 1, function(row){
  sim_params <- set_params(as.numeric(row))
  
  result <- replicate(1, expr = {
    
    simulation <- sim(150, params = sim_params)
    
    # We'll use the coefficients from this model as starting parameters:
    # beta0 = intercept; beta1 = slope
    simple_model <- lm(c(log((simulation) + 1)) ~ c(0:149))
    
    
    guess <- set_params(c(simple_model$coefficients[1],
                          simple_model$coefficients[2],
                          0, 0, 0))
    
    estimates <- optim(count_log_like_new, par = guess, constant = 0, method = "BFGS",
                       first_record_data = simulation, hessian = TRUE)
    
    return(list(timeseries = simulation, estimates = estimates))
  }, simplify = FALSE)
  
  return(list(sim_params = sim_params, simulation = result))},
  cl = 12)
