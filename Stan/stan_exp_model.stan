data{
    int<lower=1> N;
    int<lower=1> M; // no. of med species
    int dI[N];
    int dsps[N];
    int n_Inv [N]; 
    int n_Nativ [N];
    real t[N];
}

parameters {
  real b0;
  real b1;
}

transformed parameters{
}

model{
    real It [N]; 
    real P [N];
    real It0;
    
 //priors
  b0 ~ normal(0, 2);
  b1 ~ normal(0, 4);
    
    for (i in 1:N) {
      It0 = exp(b0 + b1 * t[i]);
      It[i] = fmax(It0 , (n_Inv[i]+dI[N]));
      //P[i] = (0.1+It[i]-n_Inv[i]) / ((0.1+It[i]-n_Inv[i])+(M-n_Nativ[i]));
      P[i] = (It[i]-n_Inv[i]) / ((It[i]-n_Inv[i])+(M-n_Nativ[i]));
      dI[i] ~ binomial( dsps[i] , P[i] );
    }
}

generated quantities{
    real It [N]; 
    real It0;
    real P [N];
    real Itot [N];
    real log_lik [N];
    real discov_t [N];
    int dI_rep [N];
    
    for (i in 1:N) {
      It0 = exp(b0 + b1 * t[i]);
      It[i] = fmax(It0 , (n_Inv[i]+dI[N]));
      //P[i] = (0.1+It[i]-n_Inv[i]) / ((0.1+It[i]-n_Inv[i])+(M-n_Nativ[i]));
      P[i] = (It[i]-n_Inv[i]) / ((It[i]-n_Inv[i])+(M-n_Nativ[i]));
    }
      
    for (i in 1:N){
      dI_rep[i] = binomial_rng(dsps[i], P[i]);
      discov_t[i] = sum(dI_rep[1:i]);
      Itot[i] = n_Inv[i] + (P[i]/(1-P[i]))*(M-n_Nativ[i]);
      log_lik[i] = binomial_lpmf(dI[i] | dsps[i], P[i]); 
    }
    
}
