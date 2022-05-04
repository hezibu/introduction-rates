data{
  int<lower=1> N;
  int<lower=1> M; // no. of med species
  int dI[N];
  int dsps[N];
  vector[N] n_Inv; 
  vector[N] n_Nativ;
  vector[N] t;
  int posterior_predictive;
}

transformed data {
  vector[N] deltaMt;
  
  deltaMt = M - n_Nativ;
}

parameters {
  real b0;
  real b1;
}

transformed parameters {
  vector[N] u;
  
  u = exp(b0 + b1*t);
}

model{
  
  vector[N] It;
  vector[N] P;
  
  //priors
  b0 ~ normal(0, 1);
  b1 ~ normal(0, 0.02);
  
  if (!posterior_predictive) {
    It = cumulative_sum(exp(b0 + b1 * t));
    P = (It - n_Inv)./((It - n_Inv) + deltaMt);
    dI ~ binomial(dsps, P);
  }
  
}

generated quantities{
  vector[N] It;
  vector[N] P;
  vector[N] invasives;
  
  if (!posterior_predictive) {
    It = cumulative_sum(exp(b0 + b1 * t));
    P = (It - n_Inv)./((It - n_Inv) + deltaMt);
    for (i in 1:N){
      invasives[i] = binomial_rng(dsps[i], P[i]);
    }
  }
}
