# These functions look at different metrics as described by Bolker, 2008:
# https://ms.mcmaster.ca/~bolker/emdbook/book.pdf

# Bias (accuracy): The expected difference between the estimate
# and the true value of the parameter. If you run a large number of
# simulations with a true value of d and estimate a value of ˆd for each
# one, then the bias is E[ˆd − d].

report_bias <- function(real_values, randomizations){
  all_b <- sapply(seq_len(length(randomizations)),
                  FUN = function(i)
                    randomizations[[i]]$estimates$par - real_values # ˆd − d
  )
  return(rowMeans(all_b)) # E[ˆd − d]
}

# Variance (precision): variance, or E[(ˆd − E[ˆd])^2], measures the
# variability of the point estimates (ˆd) around their mean value. Just as
# an accurate but imprecise answer is worthless, unbiased answers are
# worthless if they have high variance. With low bias we know that we
# get the right answer on average, but high variability means that any
# particular estimate could be way off. With real data, we never know
# which estimates are right and which are wrong.

report_variance <- function(real_values, randomizations){
  all_d <- sapply(seq_len(length(randomizations)),
                  FUN = function(i)
                    randomizations[[i]]$estimates$par # ˆd
  )
  
  e_all_d <- rowMeans(all_d) # E[ˆd]
  variance <- apply(all_d, MARGIN = 2, function(set)
    (set - e_all_d)^2) # (ˆd − E[ˆd])^2
  return(rowMeans(variance)) #E[(ˆd − E[ˆd])^2]
}

# Confidence interval width (precision): the width of the confidence 
# intervals, either in absolute terms or as a proportion of the estimated
# value, provides useful information on the precision of your estimate.
# If the confidence interval is estimated correctly (see coverage, below)
# then the confidence interval should be related to the variance among
# estimates.

report_ci_width <- function(real_values, randomizations){
  all_widths <- sapply(seq_len(length(randomizations)),
                       FUN = function(i) {
                         opt <- randomizations[[i]]$estimates
                         hessian <- opt$hessian
                         parameter_se <- HelpersMG::SEfromHessian(hessian) 
                         width <- 2 * 1.96 * parameter_se
                       }
  )
  return(rowMeans(all_widths))
}

# Mean squared error (MSE: accuracy and precision) combines bias and
# variance as (bias^2+variance). It represents the total variation around
# the true value, rather than the average estimated value (E[d − ˆd])^2 +
# E[(ˆd − E[ˆd])^2] = E[(ˆd − d)^2]. MSE gives an overall sense of the quality
# of the estimator.

report_mse <- function(real_values, randomizations){
  bias <- report_bias(real_values, randomizations)
  variance <- report_variance(real_values, randomizations)
  return(bias^2 + variance)
}

# Coverage (accuracy): when we sample data and estimate parameters,
# we try to estimate the uncertainty in those parameters. Coverage 
# describes how accurate those confidence intervals are, and (once again)
# can only be estimated via simulation. If the confidence intervals (for
# a given confidence level 1 − α) are d_low and d_high, then the coverage
# describes the proportion or percentage of simulations in which the
# confidence intervals actually include the true value 
# (Prob(d_low < d < d_high)). Ideally, the observed coverage should equal the 
# nominal coverage of 1 − α; values that are too high are pessimistic, 
# overstating the level of uncertainty, while values that are too low are 
# optimistic. (It often takes several hundred simulations to get a reasonably 
# precise estimate of the coverage, especially when estimating the coverage for
# 95% confidence intervals.)

report_coverage <- function(real_values, randomizations){
  all_widths <- sapply(seq_len(length(randomizations)),
                       FUN = function(i) {
                         opt <- randomizations[[i]]$estimates
                         hessian <- opt$hessian
                         parameter_se <- HelpersMG::SEfromHessian(hessian) 
                         d_low = opt$par - 1.96 * parameter_se
                         d_high = opt$par + 1.96 * parameter_se
                         return(d_low < real_values & real_values < d_high)
                       }
  )
  return(rowSums(all_widths))
}

# Power (precision): The narrow-sense power gives the probability
# of correctly rejecting the null hypothesis, or in other words the
# fraction of the times that the null-hypothesis value d_0 will be outside
# of the confidence limits: (Prob(d_0 < d_low or d_0 > d_high)). In frequentist
# language, it is 1−β, where β is the probability of making a type II
# error. 
# Not applicable here?




