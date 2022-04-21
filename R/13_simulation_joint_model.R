# Loading libraries:
library(tidyverse)
library(rstan)
library(loo)
library(shinystan)
library(brms)
library(bayesplot)
library(mobsim)
source("R/02_global_functions.R")

# loading the model:

stan_model <- readRDS("Stan/stan_exp_model.rds")

# simulation: increasing parts:

list_params <- list(
  beta0 = seq(-1.5,1.5, length.out = 5),
  beta1 = c(0, 0.002, 0.005, 0.01, 0.015, 0.02)
)

df_params <- expand.grid(list_params)

results <- vector(mode = "list", length = nrow(df_params))

results <- pbapply::pbapply(df_params, MARGIN = 1, function(row) {
  
  # Defining variables and parameters:
  time_span <- 150      # time series length.
  b0 <- row[1]          # parameters defining u as a function of t.
  b1 <- row[2]          # parameters defining u as a function of t.
  species_pool <- 1:500 # array with species names as integers (e.g, 1 is species a, 2 is species b, etc.) 
  n_enter <- NULL       # Poisson random variate, drawn from distribution with mean u, number of species entering in t
  
  u <- exp(b0 + b1*seq_len(time_span))
  
  cli::cli_alert_info("Using {.var {b0}} and {.var {b1}} for simulation")
  
  result <- replicate(4, expr = {
    
    
    # Creating an introduction
    
    simulation_data <- generate_introduction(time_span, u, species_pool)
    
    # Simulating Discoveries:
    # Assumption: probability of sampling invasives and natives = their proportion in total species pool
    # Note to self: Proportion in richness, not estimated abundance...
    
    M <- 6000             # Estimated number of native species
    mean_effort <- 50     # Mean of poisson distribution - number of species sampled each sampling
    sampling_times <- 150 # Sampling effort in terms of number of sampling done throughout the time series
    
    sampling_data <- simulate_discoveries(simulation_data, M, mean_effort, sampling_times)
    
    data_for_stan <- list(
      M = M,                                       # Estimated number of native species
      N = nrow(sampling_data),                     # Number of samples - number of rows for sampling data
      dI = sampling_data$n_new_invasives,          # Number of newly discovered invasive species
      d_Nativ = sampling_data$n_new_natives,       # Number of newly discovered native species
      dsps =  sampling_data$n_new_species,         # Number of newly discovered species (natives + invasives)
      t = sampling_data$time_span,                 # Sampling events time
      n_Inv   = sampling_data$total_new_invasives, # Total number of invasives discovered
      n_Nativ =  sampling_data$total_new_natives   # Total number of natives discovered
    )
    
    stan_samples <- sampling(object = stan_model,
                         data = data_for_stan,
                         chains = 3,     # number of Markov chains
                         cores  = 3,     # number of cores   
                         warmup = 10000,   # number of warm-up iterations per chain
                         iter = 20000,     # total number of iterations per chain
                         refresh = 0,    # show progress every 'refresh' iterations
                         thin = 2,       # take sample every 2 steps in the chain
                         control = list(adapt_delta = 0.99)
    )
    
    simulation <- sampling_data$n_new_invasives
    
    years <- seq_along(simulation)-1
    
    # We'll use the coefficients from this model as starting parameters:
    # beta0 = intercept; beta1 = slope
    simple_model <- lm(c(log((simulation) + 1)) ~ years)
    
    guess <- set_params(c(simple_model$coefficients[1],
                          simple_model$coefficients[2],
                          0, 0, 0))
    
    sc_model <- try(optim(count_log_like, par = guess, constant = 0, method = "BFGS",
                          first_record_data = simulation, hessian = TRUE), silent = T)
    
    return(list(sampling_simulation = sampling_data, stan_samples = stan_samples, sc_optim = sc_model))
  }, simplify = FALSE)
  return(list(sim_params = row, simulation = result))
})
