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
  beta1 = seq(0, 0.1, length.out = 5)
)

df_params <- expand.grid(list_params)

results <- vector(mode = "list", length = nrow(df_params))

for (row in 7:nrow(df_params)){
  
  # Defining variables and parameters:
  time_span <- 150      # time series length.
  b0 <- df_params[row,1]          # parameters defining u as a function of t.
  b1 <- df_params[row,2]          # parameters defining u as a function of t.
  species_pool <- 1:500 # array with species names as integers (e.g, 1 is species a, 2 is species b, etc.) 
  n_enter <- NULL       # Poisson random variate, drawn from distribution with mean u, number of species entering in t
  
  u <- exp(b0 + b1*seq_len(time_span))
  
  cli::cli_alert_info("Using {.var {b0}} and {.var {b1}} for simulation")
  
  result <- replicate(20, expr = {
    
    
    # Creating an introduction
    
    simulation_data <- generate_introduction(time_span, u, species_pool)
    
    # Simulating Discoveries:
    # Assumption: probability of sampling invasives and natives = their proportion in total species pool
    # Note to self: Proportion in richness, not estimated abundance...
    
    M <- 6000             # Estimated number of native species
    mean_effort <- 50    # Mean of poisson distribution - number of species sampled each sampling
    sampling_times <- 61 # Sampling effort in terms of number of sampling done throughout the time series
    
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
                         iter = 40000,     # total number of iterations per chain
                         refresh = 0,    # show progress every 'refresh' iterations
                         thin = 2,       # take sample every 2 steps in the chain
                         control = list(adapt_delta = 0.99)
    )
    
    return(list(sampling_simulation = sampling_data, samples = stan_samples))
  }, simplify = FALSE)
  
  results[[row]] <- list(sim_params = df_params[row,], simulation = result)
}

# dir.create("Full Result Lists/Stan Model/Range", recursive = TRUE)
saveRDS(results, "Full Result Lists/Stan Model/Range/paramerter range 05042022")
