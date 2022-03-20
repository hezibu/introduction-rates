# Understanding Moshe's model

# Loading libraries:
library(tidyverse)
library(rstan)
library(loo)
library(shinystan)
library(brms)
library(bayesplot)
library(mobsim)

# Defining variables and parameters:
time_span <- 1:150            # time series.
b0 <- 0            # parameters defining u as a function of t.
b1 <- 0.014           # parameters defining u as a function of t.
species_pool <- 1:200 # array with species names as integers (e.g, 1 is species a, 2 is species b, etc.) 
n_enter <- NULL       # Poisson random variate, drawn from distribution with mean u, number of species entering in t

# Defining functions:

# function returns names of new species
generate_introduction <- function(time_span, n_enter, species_pool){
  
  species_remaining <- species_pool # At the start of introduction, all species in species pool are candidates
  new_sps <- NULL
  out <- list()
  
  for (t in time_span){
    sps_enter <- sample(species_pool, n_enter[t], replace = F)     # We sample from the entire species pool
    new_species_id <- species_remaining %in% sps_enter             # get indices of new species in the remaining species pool
    new_sps <- sort(c(new_sps, species_remaining[new_species_id])) # Add this t's new species to pool of already introduced 
    species_remaining <- species_remaining[!new_species_id]        # Remove newly introduced species from remaining species pool
    
    out[[t]] <- list(introduced = new_sps,          # Identity of all species 
                     n_total = length(new_sps),     # Number of introduced species so far
                     n_new = sum(new_species_id),   # Number of introduced species in t
                     remaining = species_remaining) # Species remaining in pool
  }
  
  return(out) 
}

# Creating an introduction:

u <- exp(b0 + b1*time_span)

simulation_data <- tibble(time_span, u) %>% 
  mutate(y = sapply(u, function(x) rpois(1, x)),
         simulation = generate_introduction(time_span, y, species_pool)) %>% 
  mutate(n_new = map_int(simulation, "n_new"),
         n_total = map_int(simulation, "n_total"))

# Visualizing the introduction:

ggplot(simulation_data)+
  aes(x = time_span, y = n_total)+
  geom_point()

# Simulating Discoveries:
# Assumption: probality of sampling invasives and natives = their proportion in total species pool
# Note to self: Proportion in richness, not estimated abundance...

M <- 300 # Estimated number of unrecorded native species

sampling_times <- 61 # Sampling effort in terms of number of sampling done throughout the time series
sample_events <- sort(sample(time_span[-1], sampling_times-1, replace = F))
sampling_data <- simulation_data[c(1, sample_events), ] %>% 
  mutate(species_entered = map(simulation, "introduced")) %>% 
  select(time_span, n_total,  species_entered)

sampling_effort <- rpois(1, 15)
discovered_natives <- NULL
discovered_invasives <- NULL


sampling_data <- sampling_data %>% 
  mutate(sad = map(n_total, function(n_total) {
    sim_sad(s_pool = M + n_total, 
            n_sim = sampling_effort, 
            sad_type = "lnorm",
            sad_coef = list("meanlog" = 5, "sdlog" = 0.5))
  })) %>% 
  mutate(prob_invasive = n_total/(M + n_total),
         n_species_sample = map2(.x = sad, .y = prob_invasive,
                                 .f = function(sad, prob_invasive) 
                                   as.numeric(rmultinom(1, size = length(sad), 
                                                        prob = c(prob_invasive, 1-prob_invasive)))))

sampling_data <- cbind(sampling_data, t(data.frame(sampling_data$n_species_sample))) %>% as_tibble

sampling_data <- sampling_data %>% 
  mutate(invasives = map(`1`, function(x) sample(species_entered, x, replace = F)),
         natives = map(`2`, function(x) sample(1:M, x, replace = F))) 

for (t in seq_len(sampling_times)){
  new_species <- sampling_data$natives[t][!sampling_data$natives[t] %in% discovered_natives]
  discovered_natives <- c(discovered_natives, sampling_data$natives[t])
}
