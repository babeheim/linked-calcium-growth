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
  real<lower=0> l_ind_sigma_rate;
  real<lower=0> l_sex_sigma_rate;
  real<lower=0> b_mu;
  real<lower=0> b_sigma;
  real<lower=0> k_pop_mu;
  real<lower=0> k_pop_sigma;
  real<lower=0> k_ind_sigma_rate;
  real<lower=0> k_sex_sigma_rate;
  real<lower=0> s_rate;
}
parameters {
  real l_pop;
  real<lower=0> l_sex_sigma;
  vector[2] l_sex;
  real<lower=0> l_ind_sigma;
  vector[N_ind] l_ind;
  real<lower=0> b;
  real t0_pop;
  real<lower=0> t0_sex_sigma;
  vector[2] t0_sex;
  real<lower=0> t0_ind_sigma;
  vector[N_ind] t0_ind;
  real k_pop;
  real<lower=0> k_sex_sigma;
  vector[2] k_sex;
  real<lower=0> k_ind_sigma;
  vector[N_ind] k_ind;
  real<lower=0> s;
}
model {
  vector[N_ind] l;
  vector[N_ind] t0;
  vector[N_obs] theta;
  vector[N_ind] k;
  vector[N_obs] m;
  l_pop ~ normal(l_pop_mu, l_pop_sigma);
  l_ind_sigma ~ exponential(l_ind_sigma_rate);
  l_ind ~ std_normal();
  l_sex_sigma ~ exponential(l_sex_sigma_rate);
  l_sex ~ std_normal();
  b ~ normal(b_mu, b_sigma) T[0, ];
  t0_pop ~ normal(l_pop_mu, l_pop_sigma);
  t0_ind_sigma ~ exponential(l_ind_sigma_rate);
  t0_ind ~ std_normal();
  t0_sex_sigma ~ exponential(l_sex_sigma_rate);
  t0_sex ~ std_normal();
  k_pop ~ normal(k_pop_mu, k_pop_sigma);
  k_ind_sigma ~ exponential(k_ind_sigma_rate);
  k_ind ~ std_normal();
  k_sex_sigma ~ exponential(k_sex_sigma_rate);
  k_sex ~ std_normal();
  s ~ exponential(s_rate);
  for (j in 1:N_ind) {
    l[j] = l_pop +
      l_sex[sex[j]] * l_sex_sigma +
      l_ind[j] * l_ind_sigma;
    t0[j] = t0_pop +
      t0_sex[sex[j]] * t0_sex_sigma +
      t0_ind[j] * t0_ind_sigma;
    k[j] = exp(k_pop +
      k_sex[sex[j]] * k_sex_sigma +
      k_ind[j] * k_ind_sigma);
  }
  for (i in 1:N_obs) {
    theta[i] = b * (age_su[i] - l[pid[i]]);
    if (y[i] == 0) {
      target += log1m(inv_logit(theta[i]));
    } else {
      m[i] = k[pid[i]] * (age_su[i] - t0[pid[i]]);
      target += log(inv_logit(theta[i])) + lognormal_lpdf(y[i] | m[i], s);
    }  
  }
}
generated quantities {
  vector[N_ind] l;
  vector[N_ind] t0;
  vector[N_ind] k;
  for (j in 1:N_ind) {
    l[j] = l_pop +
      l_sex[sex[j]] * l_sex_sigma +
      l_ind[j] * l_ind_sigma;
    t0[j] = t0_pop +
      t0_sex[sex[j]] * t0_sex_sigma +
      t0_ind[j] * t0_ind_sigma;
    k[j] = exp(k_pop +
      k_sex[sex[j]] * k_sex_sigma +
      k_ind[j] * k_ind_sigma);
  }
}
