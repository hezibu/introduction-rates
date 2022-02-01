source("R/00_required_packages.R")
# Rcpp::sourceCpp("R/01_solow_costello_cpp_functions.cpp")
source("R/02_global_functions.R")

sim_params_ts_length <- set_params(c(-1.1, 0.014, -1.46, 0.00001, 0.0000004))


# Run the simulation/underlying curve using sim/count_lambda respectively
set.seed(333) # setting seed before `sim` guarantees that curve looks identical for each set of parameters  
test_tsl <- sim(N = 200, params = sim_params_ts_length)

# visualize `test_tsl`:

qplot(x = 0:199, cumsum(test_tsl), geom = "line")+
  geom_line(aes(x= 50:199, y = cumsum(test_tsl[50:199])), color = "red")

list_tsl <- list(
  beta0 = seq(-1.5,1.5, length.out = 20),
  beta1 = seq(0, 0.03, length.out = 20),
  actual_length = seq(150, 400, by = 10)
)

df_tsl <- expand.grid(list_tsl)

parts_list_tsl <- split(sample(seq_len(nrow(df_tsl))), 1:5)

parts_list_tsl <- lapply(parts_list_tsl, function(rows){
  return(df_tsl[rows,])
})

tsl_results_parts_list <- vector(mode = "list", length = length(parts_list_tsl))

parts_list_tsl <- readRDS("Results/Prior Invasion/prior invasion parameters 28012022")


cl <- parallel::makeCluster(detectCores())
plan(cluster, workers = cl)
handlers(handler_progress(format="[:bar] :percent :eta :message"))
with_progress({
  p <- progressor(along = 1:nrow(parts_list_tsl$`5`))
  tsl_results_parts_list[[5]] <-  future_apply(parts_list_tsl$`5`, MARGIN = 1, function(row){
    
    p()
    
    sim_params_ts_length[1:2] <- as.numeric(row[1:2]) 
    
    actual_length <- as.numeric(row[3])
    
    result <- replicate(1000, expr = {
      
      simulation <- sim(actual_length, params = sim_params_ts_length)
      
      years_to_cut <- actual_length - 150
      
      if (actual_length > 150) {
        simulation <- simulation[-(seq_len(years_to_cut))] # return length 150
      }
      
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
    
    return(list(sim_params_tsl = c(sim_params_ts_length, actual_length = actual_length), simulation = result))
  }, future.seed = NULL, future.chunk.size = 1)
  

})

# dir.create("Results/Prior Invasion/", recursive = TRUE)

saveRDS(tsl_results_parts_list,"Results/Prior Invasion/prior invasion simulation results 31012022")
saveRDS(parts_list_tsl,"Results/Prior Invasion/prior invasion parameters 31012022")
