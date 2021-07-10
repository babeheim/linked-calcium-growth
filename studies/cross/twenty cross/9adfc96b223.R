list(
  sex_effect = "l on, k on",
  priors = "perfect",
  ind_variation = "lots",
  sampling = "normal",
  seed = NA,
  n_runs = 10,
  sim = list(
    longitudinal = FALSE,
    N_ind        =    20,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =   0.1,
    l_sex_sigma  =   1.0,
    b_mu         = 1/1.0,
    b_sigma      =   0.5,
    k_pop_mu     =   1.0,
    k_pop_sigma  =   0.1,
    k_ind_sigma  =  0.24,
    k_sex_sigma  =   0.5,
    s            =   0.1
  ),
  stan = list(
    # prior model parameters
    l_pop_mu          = 0,
    l_pop_sigma       = 0.18,
    l_ind_sigma_rate  = 1/0.1,
    l_sex_sigma_rate  = 1/1.0,
    b_mu              = 1.0,
    b_sigma           = 0.5,
    k_pop_mu          = 1,
    k_pop_sigma       = 0.1,
    k_ind_sigma_rate  = 1/0.24,
    k_sex_sigma_rate  = 1/0.5,
    s_rate            = 1/0.1,
    # sampling parameters
    n_chains = 4,
    n_cores = 4,
    n_iter = 5000,
    adapt_delta = 0.90,
    max_treedepth = 10
  )
)
