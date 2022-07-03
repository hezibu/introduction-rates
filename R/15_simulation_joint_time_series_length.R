# Loading libraries:
library(tidyverse)
library(mobsim)
library(cmdstanr)
source("/home/hezi/introduction-rates/R/02_global_functions.R")
# cmdstanr::install_cmdstan()

# Creating a folder on the NMVe drive(s):
if (!file.exists("/data_1/Stan_csv_files/"))
  dir.create("/data_1/Stan_csv_files/")
if (!file.exists("/data_2/Stan_csv_files/"))
  dir.create("/data_2/Stan_csv_files/")

# loading the model:

mod <- cmdstan_model("/home/hezi/introduction-rates/Stan/stan_exp_improved.stan")

# simulation: increasing parts:

replicates_per_parameters <- 1000

list_params <- list(
  replicat = 1:replicates_per_parameters,
  beta0 = seq(-1.5,1.5, length.out = 3),
  beta1 = c(0, 0.01, 0.02),
  length = seq(10, 400, by = 15)
)

df_params <- expand.grid(list_params)

parts_list <-  df_params %>% 
  group_by(cut(length, 5)) %>% 
  nest() %>%
  pull(data)

names(parts_list) <- 1:5

results <- vector(mode = "list", length = nrow(df_params))



library(future)
library(future.apply)
library(progressr)
handlers("txtprogressbar")

plan("multisession", workers = 32)

cli::cli_alert_success("So far, so good... Starting the simulations.")

for (i in 4:5){
  
  cli::cli_text("Now running chunk {.var {i}}")
  
  with_progress({
    p <- progressor(along = 1:nrow(parts_list[[i]]))
    results <- future.apply::future_apply(parts_list[[i]], MARGIN = 1, FUN = function(row) {
      
      # Defining variables and parameters:
      replicate <- row[1]    # replicate for this set of parameters
      b0 <- row[2]           # parameters defining u as a function of t.
      b1 <- row[3]           # parameters defining u as a function of t.
      time_span <- row[4]       # time series length.
      
      species_pool <- 1:500 # array with species names as integers (e.g, 1 is species a, 2 is species b, etc.) 
      n_enter <- NULL       # Poisson random variate, drawn from distribution with mean u, number of species entering in t
      
      u <- exp(b0 + b1*seq_len(time_span))
      
      # Creating an introduction
      
      simulation_data <- generate_introduction(time_span, u, species_pool)
      
      # Simulating Discoveries:
      # Assumption: probability of sampling invasives and natives = their proportion in total species pool
      # Note to self: Proportion in richness, not estimated abundance...
      
      M <- 6000                # Estimated number of native species
      mean_effort <- 50        # Mean of poisson distribution - number of species sampled each sampling
      sampling_times <- row[4] # Sampling effort in terms of number of sampling done throughout the time series
      
      sampling_data <- simulate_discoveries(simulation_data, M, mean_effort, sampling_times)
      
      data_for_stan <- list(
        posterior_predictive = 0,                    # To sample from posterior
        M = M,                                       # Estimated number of native species
        N = nrow(sampling_data),                     # Number of samples - number of rows for sampling data
        dI = sampling_data$n_new_invasives,          # Number of newly discovered invasive species
        d_Nativ = sampling_data$n_new_natives,       # Number of newly discovered native species
        dsps =  sampling_data$n_new_species,         # Number of newly discovered species (natives + invasives)
        t = sampling_data$time_span,                 # Sampling events time
        n_Inv   = sampling_data$total_new_invasives, # Total number of invasives discovered
        n_Nativ =  sampling_data$total_new_natives   # Total number of natives discovered
      )
      
      json_file <- tempfile(fileext = ".json", tmpdir = "/json_files")
      write_stan_json(data_for_stan, json_file)
      
      # sample between data and data_2 for the Stan output csv
      
      stan_output_file <- sample(x = c("/data_1/Stan_csv_files/",
                                       "/data_2/Stan_csv_files/"),
                                 size = 1)
      
      samples <- suppressMessages(
        mod$sample(
          data = json_file,
          chains = 3,         # number of Markov chains
          parallel_chains = 3,# number of cores
          refresh = 0, iter_warmup = 2000, iter_sampling = 8000,
          show_messages = FALSE, output_dir = stan_output_file)
      )
      
      summary <- samples$summary()
      
      # call_csv_files <- list(draws = samples$draws(),
      #                        diagnostics = samples$sampler_diagnostics())
      
      p(sprintf("row=%g", row))
      
      
      return(list(sim_params = row, sampling_simulation = data_for_stan, stan_samples = summary))
    }, future.seed = TRUE)
  })
  
  cli::cli_alert_success("Done, now saving RDS file..")
  saveRDS("saving", "/data_list/success.RDS")
  saveRDS(results, stringr::str_glue("/data_list/joint_run_part_{i}_of_5.RDS"))
  
  
  ##############################
  
  cli::cli_alert_success("Done, now manipulating data for saving...")
  
  data_params <- map(results, function(item) item$sim_params) %>% bind_rows()
  # 
  data_tbl <- data_params %>%
    mutate(stan = map(results, function(item) item$stan_samples[1:3,])) %>%
    unnest(stan)
  
  cli::cli_alert_success("Done, now saving table in RDS..")
  
  saveRDS(data_tbl, stringr::str_glue("/results_tables/time_series_length_part_{i}_of_5.RDS"))
}

cli::cli_alert_success("Done!")