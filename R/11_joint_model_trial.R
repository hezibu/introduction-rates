# Understanding Moshe's model

# Loading libraries:
library(tidyverse)
library(rstan)
library(loo)
library(shinystan)
library(brms)
library(bayesplot)
library(mobsim)
source("R/02_global_functions.R")
rstan_options("auto_write" = TRUE)

exp_model <- stan_model("Stan/stan_exp_model.stan", model_name = "exp_model")

# Defining variables and parameters:
time_span <- 150      # time series length.
b0 <- 0            # parameters defining u as a function of t.
b1 <- 0.02          # parameters defining u as a function of t.
species_pool <- 1:500 # array with species names as integers (e.g, 1 is species a, 2 is species b, etc.) 
n_enter <- NULL       # Poisson random variate, drawn from distribution with mean u, number of species entering in t

# Creating an introduction:

u <- exp(b0 + b1*seq_len(time_span))

simulation_data <- generate_introduction(time_span, u, species_pool)

# Visualizing the introduction:

ggplot(simulation_data)+
  aes(x = time_span, y = n_total)+
  geom_point()

# Simulating Discoveries:
# Assumption: probability of sampling invasives and natives = their proportion in total species pool
# Note to self: Proportion in richness, not estimated abundance...

M <- 6000             # Estimated number of native species
mean_effort <- 50    # Mean of poisson distribution - number of species sampled each sampling
sampling_times <- 61 # Sampling effort in terms of number of sampling done throughout the time series

sampling_data <- simulate_discoveries(simulation_data, M, mean_effort, sampling_times)

sampling_data %>% 
  ggplot()  + 
  geom_point(aes(x = time_span, y = n_total), col='grey') +
  xlim(0, last(sampling_data$time_span)) + 
  ylim(0,(last(sampling_data$total_new_natives)+last(sampling_data$total_new_invasives))) + 
  ylab("No. of Species") +
  geom_point(aes(x = time_span, y = total_new_invasives), color='red') +
  geom_point(aes(x = time_span, y = total_new_natives), color='blue')

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

fit1.test <- sampling(object = exp_model,
                  data = data_for_stan,
                  chains = 3,     # number of Markov chains
                  cores  = 3,     # number of cores   
                  warmup = 2000,   # number of warm-up iterations per chain
                  iter = 8000,     # total number of iterations per chain
                  refresh = 0,    # show progress every 'refresh' iterations
                  thin = 2,       # take sample every 2 steps in the chain
                  control = list(adapt_delta = 0.99)
)

# ALL OF THE BELOW NEED TO BE CHANGED:

fitted.data = rstan::extract(fit1.test, permuted=TRUE)
sampling_data$av.Itot=apply(fitted.data$Itot, 2, median)
sampling_data$av.dI=apply(fitted.data$dI_rep, 2, median)
sampling_data$av.discov_t=apply(fitted.data$discov_t, 2, median)

sampling_data$quant.01.Itot=apply(fitted.data$Itot, 2, 
                             function (x) quantile(x, probs=c(0.01)))
sampling_data$quant.99.Itot=apply(fitted.data$Itot, 2, 
                             function (x) quantile(x, probs=c(0.99)))                   
sampling_data$quant.01.dI=apply(fitted.data$dI_rep, 2, 
                           function (x) quantile(x, probs=c(0.01)))
sampling_data$quant.99.dI=apply(fitted.data$dI_rep, 2, 
                           function (x) quantile(x, probs=c(0.99)))   
sampling_data$quant.01.discov_t=apply(fitted.data$discov_t, 2, 
                                 function (x) quantile(x, probs=c(0.01)))
sampling_data$quant.99.discov_t=apply(fitted.data$discov_t, 2, 
                                 function (x) quantile(x, probs=c(0.99))) 

f.Itot = ggplot(data = sampling_data)  + 
  geom_point(aes(x = time_span, y = total_new_invasives), col='black') +
  geom_line(aes(x = time_span, y = av.Itot), col='red') +
  xlim(1,150) #+ ylim(0,length(species_pool)+20) + ylab("Introduced Species") 

f.Itot = f.Itot + geom_ribbon(data = sampling_data, 
                              aes(x = time_span, ymin=quant.01.Itot, ymax=quant.99.Itot), 
                              alpha=0.2, fill='red') 
f.Itot

f.discov_t = ggplot(data = sampling_data)  + 
  geom_point(aes(x = time_span, y = total_new_invasives), col='black') +
  geom_line(aes(x = time_span, y = av.discov_t), col='red') +
  xlim(1,150) + ylab("Discovered Inv Species") 

f.discov_t = f.discov_t + geom_ribbon(data = sampling_data, 
                                      aes(x = time_span, ymin=quant.01.discov_t, ymax=quant.99.discov_t), 
                                      alpha=0.2, fill='red') 
f.discov_t

shinystan::launch_shinystan(object = fit1.test)
