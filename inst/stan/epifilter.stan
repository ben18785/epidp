data {
  int N; // number of data points
  int C[N]; // case series
  int wmax; // max day of generation time distribution
  row_vector[wmax] w;
}

parameters {
  real<lower=0> R[N];
  real<lower=0> sigma;
}

transformed parameters {

  real E_cases[N];

  {
    vector[wmax] I_temp;
    for(t in 1:N) {
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
    R[t] ~ normal(R[t - 1], sigma * sqrt(R[t]));
  }

  R[1] ~ uniform(0, 50);
  sigma ~ cauchy(0, 1);
}

generated quantities {
  vector[N] log_likelihood;
  for(t in 1:N)
    log_likelihood[t] = poisson_lpmf(C[t]|E_cases[t]);
}
