data {
  int N; // number of data points
  int C[N]; // case series
  int wmax; // max day of generation time distribution
  row_vector[wmax] w; // generation time distribution
  int N_covariates;
  matrix[N, N_covariates] X; // matrix's first column is 1s
}

parameters {
  real<lower=0> sigma;
  vector[N_covariates] beta;
  real<lower=0> epsilon[N];
}

transformed parameters {

  real E_cases[N];
  real<lower=0> R[N];

  {
    vector[wmax] I_temp;
    for(t in 1:N) {
      R[t] = exp(X[t] * beta + epsilon[t]);
      if(t == 1) {
        I_temp = rep_vector(0, wmax);
      } else if(t <= wmax) {
        int kk = wmax - t + 1;
        for(i in 1:(t - 1))
          I_temp[i] = C[t - i]; // needs to be lagged by one time point
        for(i in 1:kk)
          I_temp[i + t - 1] = 0;
      } else {
        for(i in 1:wmax)
          I_temp[i] = C[t - i]; // needs to be lagged by one time point
      }
      E_cases[t] = R[t] * w * I_temp;
    }
  }
}


model {
  for(t in 2:N) {
    C[t] ~ poisson(E_cases[t]);
    epsilon[t] ~ normal(epsilon[t - 1], sigma);
  }

  epsilon[1] ~ cauchy(0, 2);
  sigma ~ cauchy(0, 1);
  beta ~ normal(0, 5);
}

generated quantities {
  vector[N] log_likelihood;
  for(t in 1:N)
    log_likelihood[t] = poisson_lpmf(C[t]|E_cases[t]);
}
