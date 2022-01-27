# This scripts loads the necessary libraries and functions

source("R/required packages.R")

# Functions:

# Functions based off of Solow & Costello, 2004:
Rcpp::sourceCpp("R/solow_costello_rcpp_functions.cpp")

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
  constant <- c()
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

