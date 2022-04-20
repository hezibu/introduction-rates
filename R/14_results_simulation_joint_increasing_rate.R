
library(pbapply)


rm(results)

#Functions:
{
  summarize_stan_trials <- function(a_single_run, n_params, fn){
  real_values <- as.numeric(a_single_run[[1]][1:n_params])
  real_values_tbl <- bind_rows(a_single_run[[1]])
  randomizations <- a_single_run[[2]]
  out <- bind_rows(fn(real_values, randomizations)) %>%
    rename_with(.fn = function(x) str_glue("{x}_fn"))
  return(bind_cols(real_values_tbl, out))
}

report_stan_bias <- function(real_values, randomizations){
  all_b <- sapply(seq_len(length(randomizations)),
                  FUN = function(i) {
                    pars <- randomizations[[i]]$samples %>% summary() %>% .$summary %>% .[,1]
                    out <- pars[-length(pars)] - real_values # ˆd − d
                    return(out)
                    }
  )
  return(rowMeans(all_b)) # E[ˆd − d]
}

report_stan_variance <- function(real_values, randomizations){
  all_d <- sapply(seq_len(length(randomizations)),
                  FUN = function(i) {
                    pars <- randomizations[[i]]$samples %>% summary() %>% .$summary %>% .[,1] # ˆd
                    pars <- pars[-length(pars)]
                    return(pars)
                  }
  )
  
  e_all_d <- rowMeans(all_d) # E[ˆd]
  variance <- apply(all_d, MARGIN = 2, function(set)
    (set - e_all_d)^2) # (ˆd − E[ˆd])^2
  return(rowMeans(variance)) #E[(ˆd − E[ˆd])^2]
}

report_stan_ci_width <- function(real_values, randomizations){
  all_widths <- sapply(seq_len(length(randomizations)),
                       FUN = function(i) {
                         ci_low <- randomizations[[i]]$samples %>% summary() %>% .$summary %>% .[,1]
                         ci_high <- randomizations[[i]]$samples %>% summary() %>% .$summary %>% .[,8]
                         width <- ci_high - ci_low
                         return(width[-length(width)])
                       }
  )
  return(rowMeans(all_widths))
}

report_stan_mse <- function(real_values, randomizations){
  bias <- report_stan_bias(real_values, randomizations)
  variance <- report_stan_variance(real_values, randomizations)
  return(bias^2 + variance)
}

report_stan_coverage <- function(real_values, randomizations){
  all_widths <- sapply(seq_len(length(randomizations)),
                       FUN = function(i) {
                         summary_table <- randomizations[[i]]$samples %>% summary() %>% .$summary
                         summary_table <- summary_table[-nrow(summary_table),]
                         d_low <- summary_table[,1]
                         d_high <- summary_table[,8]
                         return(d_low < real_values & real_values < d_high)
                       }
  )
  return(rowSums(all_widths)/ncol(all_widths)) #change later when no differing number of randomizations
}

report_stan_rhat <-  function(real_values, randomizations){
  all_b <- sapply(seq_len(length(randomizations)),
                  FUN = function(i) {
                    pars <- randomizations[[i]]$samples %>% summary() %>% .$summary %>% .[,ncol(.)]
                    out <- pars[-length(pars)] # ˆd − d
                    return(out)
                  }
  )
return(apply(all_b,MARGIN = 1, mean)) # E[ˆd − d]
}

results <- readRDS("Full Result Lists/Stan Model/Range/paramerter range 07042022")

biases_1 <- pbapply::pblapply(results, function(trial) {
  return(try(summarize_stan_trials(trial, 2, report_stan_bias)))
})

biases_1 <- biases_1[sapply(biases_1, function(x) all(class(x) != "try-error"))] %>% bind_rows()

variances_1 <- pblapply(results, function(trial) {
  return(try(summarize_stan_trials(trial, 2, report_stan_variance)))
})

variances_1 <- variances_1[sapply(variances_1, function(x) all(class(x) != "try-error"))] %>% bind_rows()

ci_widths_1 <-  pblapply(results, function(trial) {
  return(try(summarize_stan_trials(trial, 2, report_stan_ci_width)))
})

ci_widths_1 <- ci_widths_1[sapply(ci_widths_1, function(x) all(class(x) != "try-error"))] %>% bind_rows()

mses_1 <-  pblapply(results, function(trial) {
  return(try(summarize_stan_trials(trial, 2, report_stan_mse)))
})

mses_1 <- mses_1[sapply(mses_1, function(x) all(class(x) != "try-error"))] %>% bind_rows()

coverages_1 <-  pblapply(results, function(trial) {
  return(try(summarize_stan_trials(trial, 2, report_stan_coverage)))
})

coverages_1 <- coverages_1[sapply(coverages_1, function(x) all(class(x) != "try-error"))] %>% bind_rows()

rhat_1 <- pbapply::pblapply(results, function(trial) {
  return(try(summarize_stan_trials(trial, 2, report_stan_rhat)))
})


biases <- bind_rows(biases_2, biases_1)
variances <- bind_rows(variances_2, variances_1)
ci_widths <- bind_rows(ci_widths_2, ci_widths_1)
mses <- bind_rows(mses_2, mses_1)
coverages <- bind_rows(coverages_2, coverages_1)

summed_result <- list(bias = biases,
                      variance = variances,
                      ci_width = ci_widths, 
                      mse = mses, 
                      coverage = coverages)

ggplot(biases)+
  aes(x = as.factor(beta0), y = as.factor(beta1), fill = b1_fn)+
  geom_tile(color = "black")+
  scale_fill_viridis_b()

ggplot(variances)+
  aes(x = as.factor(beta0), y = as.factor(beta1), fill = b1_fn)+
  geom_tile(color = "black")+
  scale_fill_viridis_b()

ggplot(mses)+
  aes(x = as.factor(beta0), y = as.factor(beta1), fill = b1_fn)+
  geom_tile(color = "black")+
  scale_fill_viridis_b()

ggplot(ci_widths)+
  aes(x = as.factor(beta0), y = as.factor(beta1), fill = b1_fn)+
  geom_tile(color = "black")+
  scale_fill_viridis_b()

ggplot(coverages)+
  aes(x = as.factor(beta0), y = as.factor(beta1), fill = b1_fn)+
  geom_tile(color = "black")+
  scale_fill_viridis_b()


