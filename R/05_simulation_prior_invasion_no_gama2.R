source(here::here("R","02_global_functions.R"))

source("R/02_global_functions.R")
Rcpp::sourceCpp("R/01_solow_costello_cpp_functions.cpp")

sim_params_ts_length <- set_params(c(-1.1, 0.014, -1.46, 0.00001, 0.0000004))


# Run the simulation/underlying curve using sim/count_lambda respectively
set.seed(333) # setting seed before `sim` guarantees that curve looks identical for each set of parameters  
test_tsl <- sim(N = 200, params = sim_params_ts_length)

# visualize `test_tsl`:

ggplot2::qplot(x = 0:199, cumsum(test_tsl), geom = "line")+
  ggplot2::geom_line(ggplot2::aes(x= 50:199, y = cumsum(test_tsl[50:199])), color = "red")

list_tsl <- list(
  beta0 = seq(-1.5,1.5, length.out = 20),
  beta1 = seq(0, 0.03, length.out = 20),
  actual_length = seq(150, 400, by = 10)
)

df_tsl <- expand.grid(list_tsl)

parts_list_tsl <- split(sample(seq_len(nrow(df_tsl))), 1:1000)

parts_list_tsl <- lapply(parts_list_tsl, function(rows){
  return(df_tsl[rows,])
})


# 
# # Set cluster type
# 
# clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
# 
# # Create cluster
# 
# clust <- snow::makeCluster(getOption("cl.cores", 7), type = clusterType)
# 
# # Load the params, and function to the cluster
# 
# 
# clusterEvalQ(clust, source("R/functions.R"))
# clusterEvalQ(clust, source("R/simulations functions_modified.R"))
# clusterEvalQ(clust, Rcpp::sourceCpp("R/creating_new_p.cpp"))
# clusterExport(clust, deparse(substitute(sim_params_ts_length)), envir = environment())



# We'll create list for partial results:

tsl_results_parts_list <- vector(mode = "list", length = nrow(df_tsl))



library(future.apply)
cl <- parallel::makeCluster(detectCores())
plan(cluster, workers = cl)

library(progressr)
handlers(handler_progress(format="[:bar] :percent :eta :message"))

with_progress({
  p <- progressor(steps = nrow(df_tsl))
  tsl_results_parts_list <-  future_apply(df_tsl, MARGIN = 1, function(row){
    
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

saveRDS(tsl_results_parts_list,"Results/Time series length no gama 2/prior invasion simulation results 20012022")
saveRDS(parts_list_tsl,"Results/Time series length no gama 2/prior invasion parameters 20012022")




tsl_results_parts_list <- tsl_results_parts_list[1:3]

test_func_tsl <- function(list_of_results, i){
  list_of_hessians <- map(map(list_of_results[[2]], 2), 6) # for the hessians
  
  list_of_trials <- lapply(1:100, function(j){
    try(
      100*(map(list_of_results[[2]], 2)[[j]]$par - list_of_results$sim_params_tsl[1:4])/list_of_results$sim_params_tsl[1:4],
      silent = TRUE)
  })
  
  indices <- sapply(list_of_trials,function(x) all(class(x) != "try-error"))
  
  list_of_trials <- list_of_trials[indices]
  
  if (length(list_of_trials) == 0) {
    mean_bias <- rep(NA,4)
    names(mean_bias) <- c("beta0_bias","beta1_bias","gama0_bias","gama1_bias")
    mean_bias <- bind_rows(mean_bias)
  }else{
    mean_bias <- list_of_trials %>% bind_rows %>% colMeans() %>% bind_rows() %>% 
      rename_with(.fn = function(x) str_glue("{x}_bias"))
  }
  
  row <- bind_cols(bind_rows(list_of_results$sim_params_tsl),mean_bias)
  
  return(row)
  
}



clusterExport(clust, deparse(substitute(test_func_tsl)), envir = environment())
clusterExport(clust, deparse(substitute(calculate_percentage)), envir = environment())


mmm <- unlist(tsl_results_parts_list, recursive = FALSE,use.names = TRUE)

aaa <- pbapply::pblapply(mmm, function(row){
  return(test_func_tsl(row, names(row)))
})



tsl_success <- bind_rows(aaa, .id = "trial") %>% 
  dplyr::select(-c(gama0, gama1, gama2)) %>% 
  pivot_longer(cols = c(beta0_bias, beta1_bias, gama0_bias,  gama1_bias), names_to = "parameter", values_to = "bias") %>% 
  rename("years_cut" = years_cut.trim)


tsl_success %>% 
  filter(years_cut < 100) %>% 
  ggplot()+
  aes(x = years_cut, y = bias, group = parameter, color = parameter)+
  stat_summary(geom = 'line')+
  scale_color_viridis_d()+
  facet_wrap(~parameter,scales = "free")


this <- mmm[[1]]

(this[2]) %>% class()

did_optim_succeed(map(this[[2]],2)[[1]], this$sim_params_gap[1:4])


