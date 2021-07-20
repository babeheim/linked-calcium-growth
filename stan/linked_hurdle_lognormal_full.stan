data {
  int<lower=0> N_ind;
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
  real d_pop_mu;
  real<lower=0> d_pop_sigma;
  real<lower=0> d_ind_sigma_rate;
  real<lower=0> d_sex_sigma_rate;
  real<lower=0> s_rate;
}
parameters {
  real l_pop;
  real<lower=0> l_sex_sigma;
  vector[2] l_sex;
  real<lower=0> l_ind_sigma;
  vector[N_ind] l_ind;
  real<lower=0> b;
  vector[N_ind] t0;
  real d_pop;
  real<lower=0> d_sex_sigma;
  vector[2] d_sex;
  real<lower=0> d_ind_sigma;
  vector[N_ind] d_ind;
  real<lower=0> s;
}
model {
  vector[N_ind] l;
  vector[N_obs] theta;
  vector[N_ind] d;
  vector[N_ind] k;
  vector[N_obs] m;
  l_pop ~ normal(l_pop_mu, l_pop_sigma);
  l_ind_sigma ~ exponential(l_ind_sigma_rate);
  l_ind ~ std_normal();
  l_sex_sigma ~ exponential(l_sex_sigma_rate);
  l_sex ~ std_normal();
  b ~ normal(b_mu, b_sigma) T[0, ];
  d_pop ~ normal(d_pop_mu, d_pop_sigma);        // on the log scale
  d_ind_sigma ~ exponential(d_ind_sigma_rate);  // on the log scale
  d_ind ~ std_normal();                         // on the log scale
  d_sex_sigma ~ exponential(d_sex_sigma_rate);  // on the log scale
  d_sex ~ std_normal();                         // on the log scale
  s ~ exponential(s_rate);
  for (j in 1:N_ind) {
    l[j] = l_pop +
      l_sex[sex[j]] * l_sex_sigma +
      l_ind[j] * l_ind_sigma;
    d[j] = exp(d_pop +
      d_sex[sex[j]] * d_sex_sigma +
      d_ind[j] * d_ind_sigma);
    k[j] = log(2)/d[j];
  }
  t0 ~ logistic(l, 1/b);
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
  vector[N_ind] d;
  vector[N_ind] k;
  for (j in 1:N_ind) {
    l[j] = l_pop +
      l_sex[sex[j]] * l_sex_sigma +
      l_ind[j] * l_ind_sigma;
    d[j] = exp(d_pop +
      d_sex[sex[j]] * d_sex_sigma +
      d_ind[j] * d_ind_sigma);
    k[j] = log(2)/d[j];
  }
}
