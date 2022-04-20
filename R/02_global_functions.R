# Functions:

# Functions based off of Solow & Costello, 2004:
# Rcpp::sourceCpp("R/solow_costello_rcpp_functions.cpp")

# Functions:

# Function to name our paramaters.
#  Input: numeric values, and names if specifying the parameters
# Output: Named numeric vector
set_params <- function(numeric_vector,parameters = c("beta0","beta1","gama0","gama1","gama2")){
  return(set_names(x = numeric_vector,nm = parameters))
}

# Load functions for precision and accuracy reporting:
source("R/03_precision_and_accuracy_functions.R")

# Function to create a Solow and Costello time series simulation.
#  Input: N- number of years and simulation parameters.
# Output: A vector of length N with yearly number of new introduced species

sim <- function(N, params){
  constant <- 0
  lambda <- vector(mode = "numeric",length = N)
  Am = count_m(0:(N-1),params,constant) 
  Ap = rev(count_p_no_gama2(0:(N-1),params,constant))
  for(i in 0:(N-1)){
    S=c(0:i)
    Yt = rbinom(Am[S],1,rev(Ap[S]))
    lambda[[i+1]] = round(sum(Am[S] * Yt))
  }
  
  return(lambda)
}

# Function to create a sequence from its center value
#  Input: center value, proportion of lower and upper values in relation to center,
#         and number of steps from the center.
# Output: A sequence of length steps_from_center + 1

seq_from_center <- function(center, range, steps_from_center){
  first_part <- seq(from = center, to = center + center * range, length.out = steps_from_center + 1)
  second_part <- seq(from = center - (center * range), to =  center , length.out = steps_from_center + 1)
  
  return(sort(unique(c(first_part,second_part))))
}

# Function to summarize te simulation results given a reporting function.
#  Input: A single element of a list, with a parameters and a simulations object
# Output: a row of a tibble giving the real parameters and the reported values

summarize_trials <- function(a_single_run, n_params, fn){
  real_values <- a_single_run[[1]][1:n_params]
  real_values_tbl <- bind_rows(a_single_run[[1]])
  randomizations <- a_single_run[[2]]
  out <- bind_rows(fn(real_values, randomizations)) %>% 
    rename_with(.fn = function(x) str_glue("{x}_fn"))
  return(bind_cols(real_values_tbl, out))
}


# Function to sample normally if n < length(x), and take all of x if n >= length(x)
#  Input: similar to `base::sample`: vector, sample size, replace (T/F), prob
# Output: a sample of size n if n < length(x), and x otherwise

sample_from_pool <- function(x, size, replace = FALSE, prob = NULL){
  if (size < length(x)){
    return(sample(x, size, replace = F))  
  } else {
    return(x)
  }
}

# Function to generate an introduction of invasive species from alien species pool
#  Input: time_span = number of t's of introduction as integer,
#         u = yearly mean introduction rate - vector of length time_span,
#         species_pool = alien species pool, numeric vector where each number represents a species
# Output: a tibble containing the identities of species entered, number of new species each t, 
#         total number of introduced species so far, along with the input parameters.
#   Note: u represents number of species introduces, not necessarily number of NEW species introduced.

generate_introduction <- function(time_span, u, species_pool){
  simulation_data <- tibble(time_span = seq_len(time_span), u) %>% 
    mutate(species_pool = list(species_pool)) %>% 
    mutate(y = sapply(u, function(x) rpois(1, x))) %>% 
    mutate(introduced = map2(.x = species_pool, .y = y, 
                             .f = function(x,y) sample_from_pool(x, y)))
  
  new_invasives_introduced <- vector(mode = "list", length = time_span)
  new_invasives_introduced[[1]] <- simulation_data$introduced[[1]]
  for (i in 2:length(new_invasives_introduced)){
    new_invasives_introduced[[i]] <- unique(sort(c(simulation_data$introduced[[i]], new_invasives_introduced[[i-1]])))
  }
  
  simulation_data <- simulation_data %>% 
    mutate(species_entered = new_invasives_introduced,
           n_total = map_int(new_invasives_introduced, length),
           n_new = ifelse(time_span == 1, n_total, n_total - lag(n_total)))
  
  return(simulation_data)
}

# Function to simulate the discoveries given a simulated introduction.
#  Input: simulation_data = Simulation tibble from `generate_introduction`, 
#         M = number of natives, 
#         mean_effort = Mean of poisson distribution - number of species sampled each sampling
#         sampling_times = Sampling effort in terms of number of sampling done throughout the time series
# Output: A tibble containing many simulated traits. 
#         Mainly the number of newly discovered species (natives, invasives, total), their identities,
#         and their corresponding cumulative sums.

simulate_discoveries <- function(simulation_data, M, mean_effort, sampling_times){
  
  # Step 1: simulate sampling events:
  
  sample_events <- sort(sample(seq_len(nrow(simulation_data))[-1], sampling_times-1, replace = F)) # Ts of sampling
  
  sampling_data <- simulation_data[c(1, sample_events), ] %>% 
    mutate(M_unknown = M) %>% 
    select(time_span, n_total,  species_entered, M_unknown) %>% 
    mutate(sampling_effort = rpois(sampling_times, mean_effort))
  
  # Step 2: Use SAD to simulate species discoveries, and use proportion of invasives
  #         for probability of discovery. 
  
  sampling_data <- sampling_data %>% 
    mutate(sad = pmap(.l = list(n_total, M_unknown, sampling_effort),
                      .f = function(n_total, M, sampling_effort) {
                        sim_sad(s_pool = M + n_total, 
                                n_sim = sampling_effort, 
                                sad_type = "lnorm",
                                sad_coef = list("meanlog" = 5, "sdlog" = 0.5))
                      })) %>% 
    mutate(prob_invasive = n_total/(M_unknown + n_total),
           n_species_sample = map2(.x = sad, .y = prob_invasive,
                                   .f = function(sad, prob_invasive) 
                                     as.numeric(rmultinom(1, size = length(sad), 
                                                          prob = c(prob_invasive, 1-prob_invasive)))))
  
  sampling_data <- cbind(sampling_data, t(data.frame(sampling_data$n_species_sample))) %>% as_tibble
  
  sampling_data <- sampling_data %>% 
    mutate(invasives = map2(.x = `1`, .y = species_entered, function(n, pool)  sample_from_pool(pool, n)),
           natives = map(`2`, function(x) sample_from_pool(1:M, x)))
  
  # Step 3: From the discovered species, simulate which are new and which were already discovered in previous t's.
  
  new_invasives_discovered <- vector(mode = "list", length = length(sampling_data$invasives))
  if (is.null(sampling_data$invasives[[1]])){
    new_invasives_discovered[1] <- list(integer(0L))
  } else {
    new_invasives_discovered[1] <- sampling_data$invasives[1]
  }
  for (i in 2:length(new_invasives_discovered)){
    new_invasives_discovered[[i]] <- unique(sort(c(sampling_data$invasives[[i]], new_invasives_discovered[[i-1]])))
  }
  
  new_natives_discovered <- vector(mode = "list", length = length(sampling_data$natives))
  if (is.null(sampling_data$natives[[1]])){
    new_natives_discovered[1] <- list(integer(0L))
  } else {
    new_natives_discovered[1] <- sampling_data$natives[1]
  }
  for (i in 2:length(new_natives_discovered)){
    new_natives_discovered[[i]] <- unique(sort(c(sampling_data$natives[[i]], new_natives_discovered[[i-1]])))
  }
  
  # Step 4: Join all the simulations in a tibble
  
  sampling_data <- sampling_data %>% 
    mutate(new_invasives_discovered,
           new_natives_discovered,
           total_new_invasives = map_int(new_invasives_discovered, length),
           total_new_natives = map_int(new_natives_discovered, length),
           n_new_invasives = ifelse(time_span == 1, total_new_invasives, total_new_invasives - lag(total_new_invasives)),
           n_new_natives = ifelse(time_span == 1, total_new_natives, total_new_natives - lag(total_new_natives))) %>% 
    mutate(n_new_species = n_new_natives + n_new_invasives)
  
  # Step 5: Shift discoveries so total number of discovered species at t+1 equals to number discovered at t
  
  sampling_data <- sampling_data %>% 
    mutate(total_new_invasives = lag(total_new_invasives),
           total_new_natives = lag(total_new_natives)) %>% 
    replace_na(list(total_new_invasives = 0, total_new_natives = 0))

  return(sampling_data)
}

