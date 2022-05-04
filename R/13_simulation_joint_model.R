# Loading libraries:
library(tidyverse)
library(rstan)
library(loo)
library(shinystan)
library(brms)
library(bayesplot)
library(mobsim)
library(snc2022)
library(cmdstanr)
source("R/02_global_functions.R")

# loading the model:

mod <- cmdstan_model("Stan/stan_exp_improved.stan")

# simulation: increasing parts:

list_params <- list(
  beta50 = seq(-1.5,1.5, length.out = 5),
  beta1 = c(0, 0.005, 0.01, 0.015, 0.02)
)

df_params <- expand.grid(list_params)

results <- vector(mode = "list", length = nrow(df_params))
replicates_per_parameters <- 1000

library(future)
library(progressr)
handlers("txtprogressbar")
plan(multisession(workers = 12))

with_progress({
  timer <- 1
  p <- progressor(along = 1:(replicates_per_parameters*nrow(df_params)))
  results <- future.apply::future_apply(df_params, MARGIN = 1, FUN = function(row) {
    
    # Defining variables and parameters:
    time_span <- 150      # time series length.
    b0 <- row[1]           # parameters defining u as a function of t.
    b1 <- row[2]           # parameters defining u as a function of t.
    species_pool <- 1:500 # array with species names as integers (e.g, 1 is species a, 2 is species b, etc.) 
    n_enter <- NULL       # Poisson random variate, drawn from distribution with mean u, number of species entering in t
    
    u <- exp(b0 + b1*seq_len(time_span))
    
    
    result <- replicate(replicates_per_parameters, expr = {
      
      
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
      
      data_for_stan$posterior_predictive <- 0
      samples <- mod$sample(
        data = data_for_stan,
        chains = 3,     # number of Markov chains
        parallel_chains = 3,     # number of cores
        refresh = 0, iter_warmup = 2000, iter_sampling = 8000,
        show_messages = FALSE)
      timer <- timer + 1
      p(sprintf("timer=%g", timer))
      
      return(list(sampling_simulation = sampling_data, stan_samples = samples))
    }, simplify = FALSE)
   
    return(list(sim_params = row, simulation = result))
  }, future.seed = TRUE)
})
