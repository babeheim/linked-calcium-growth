data {
  int<lower=1> N_ind;
  int sex[N_ind];
  int<lower=0> N_obs;
  int pid[N_obs];
  real age_su[N_obs];
  real<lower=0> y[N_obs];
  // add priors as input
  real l_pop_mu;
  real<lower=0> l_pop_sigma;
  real<lower=0> l_sex_sigma_rate;
  real<lower=0> b_mu;
  real<lower=0> b_sigma;
}
parameters {
  real l_pop;
  real<lower=0> l_sex_sigma;
  vector[2] l_sex;
  real<lower=0> b;
}
model {
  vector[N_ind] l;
  vector[N_obs] theta;
  l_pop ~ normal(l_pop_mu, l_pop_sigma);
  l_sex_sigma ~ exponential(l_sex_sigma_rate);
  l_sex ~ std_normal();
  b ~ normal(b_mu, b_sigma) T[0, ];
  for (j in 1:N_ind) {
    l[j] = l_pop +
      l_sex[sex[j]] * l_sex_sigma;
  }
  for (i in 1:N_obs) {
    theta[i] = b * (age_su[i] - l[pid[i]]);
    if (y[i] == 0) {
      target += log1m(inv_logit(theta[i]));
    } else {
      target += log(inv_logit(theta[i]));
    }  
  }
}
generated quantities {
  vector[N_ind] l;
  for (j in 1:N_ind) {
    l[j] = l_pop +
      l_sex[sex[j]] * l_sex_sigma;
  }
}
