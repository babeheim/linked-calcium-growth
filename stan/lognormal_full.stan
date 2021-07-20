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
  real d_pop_mu;
  real<lower=0> d_pop_sigma;
  real<lower=0> d_ind_sigma_rate;
  real<lower=0> d_sex_sigma_rate;
  real<lower=0> s_rate;
}
parameters {
  real t0_pop;
  real<lower=0> t0_sex_sigma;
  vector[2] t0_sex;
  real<lower=0> t0_ind_sigma;
  vector[N_ind] t0_ind;
  real d_pop;
  real<lower=0> d_sex_sigma;
  vector[2] d_sex;
  real<lower=0> d_ind_sigma;
  vector[N_ind] d_ind;
  real<lower=0> s;
}
model {
  vector[N_ind] t0;
  vector[N_ind] d;
  vector[N_ind] k;
  vector[N_obs] m;
  t0_pop ~ normal(l_pop_mu, l_pop_sigma);
  t0_ind_sigma ~ exponential(l_ind_sigma_rate);
  t0_ind ~ std_normal();
  t0_sex_sigma ~ exponential(l_sex_sigma_rate);
  t0_sex ~ std_normal();
  d_pop ~ normal(d_pop_mu, d_pop_sigma);
  d_ind_sigma ~ exponential(d_ind_sigma_rate);
  d_ind ~ std_normal();
  d_sex_sigma ~ exponential(d_sex_sigma_rate);
  d_sex ~ std_normal();
  s ~ exponential(s_rate);
  for (j in 1:N_ind) {
    t0[j] = t0_pop +
      t0_sex[sex[j]] * t0_sex_sigma +
      t0_ind[j] * t0_ind_sigma;
    d[j] = exp(d_pop +
      d_sex[sex[j]] * d_sex_sigma +
      d_ind[j] * d_ind_sigma);
    k[j] = log(2)/d[j];
  }
  for (i in 1:N_obs) {
    if (y[i] != 0) {
      m[i] = k[pid[i]] * (age_su[i] - t0[pid[i]]);
      target += lognormal_lpdf(y[i] | m[i], s);
    }
  }
}
generated quantities {
  vector[N_ind] t0;
  vector[N_ind] d;
  vector[N_ind] k;
  for (j in 1:N_ind) {
    t0[j] = t0_pop +
      t0_sex[sex[j]] * t0_sex_sigma +
      t0_ind[j] * t0_ind_sigma;
    d[j] = exp(d_pop +
      d_sex[sex[j]] * d_sex_sigma +
      d_ind[j] * d_ind_sigma);
    k[j] = log(2)/d[j];
  }
}
