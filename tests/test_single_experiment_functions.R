
test_that("calc_key_estimates works", {

  pars <- list(
    seed = 1,
    sim = list(
      longitudinal =  TRUE,
      N_ind        =    10,
      l_pop_mu     =     0,
      l_pop_sigma  =  0.18,
      l_ind_sigma  =   0.5,
      l_sex_sigma  =   1.0,
      b_mu         =   1.5, # b = 1 / v, v is a scale term in decades
      b_sigma      =   0.5,
      d_pop_mu     = -0.07, # on the log scale
      d_pop_sigma  =  0.01,
      d_ind_sigma  =  0.04,
      d_sex_sigma  =  0.01,
      s            =   0.1
    ),
    stan = list(
      # prior model parameters
      l_pop_mu          = 0,
      l_pop_sigma       = 0.18,
      l_ind_sigma_rate  = 1 / 0.5,
      l_sex_sigma_rate  = 1 / 1.0,
      b_mu              = 1.0,
      b_sigma           = 0.5,
      d_pop_mu          = -0.7,
      d_pop_sigma       = 0.01,
      d_ind_sigma_rate  = 1 / 0.04,
      d_sex_sigma_rate  = 1 / 0.01,
      s_rate            = 1 / 0.1,
      # sampling parameters
      n_chains = 3,
      n_cores = 3,
      n_iter = 200,
      adapt_delta = 0.8,
      max_treedepth = 15
    )
  )

  pars$start_time <- Sys.time()

  ######### simulate data #########

  pars$sim$seed <- pars$seed
  sim <- sim_cac(pars$sim)
  pars$sim <- sim$pars # save updated simulation pars

  ######### extract / prep data for analysis #########

  data <- list()
  data$obs <- sim$obs
  data$people <- sim$people

  data$obs$age_su <- (data$obs$age - 50) / 10

  ######### fit stan models #########

  stan_data <- c(
    pars$stan,
    list(
      N_ind = nrow(data$people),
      sex = as.array(data$people$sex),
      N_obs = nrow(data$obs),
      pid = as.array(data$obs$pid),
      age_su = as.array(data$obs$age_su),
      y = as.array(data$obs$cac)
    )
  )

  estimates <- list()
  fitted_models <- list()
  summaries <- list()
  diagnostics <- list()

  stan_pars <- c("l_pop", "l_sex_sigma", "t0_pop", "t0_sex_sigma",
    "t0_sex", "k_pop", "k_sex_sigma", "k_sex", "s")
  init <- vector("list", pars$stan$n_chains)
  init_one <- pars$sim[names(pars$sim) %in% stan_pars]
  for (i in seq_along(init)) init[[i]] <- init_one

  # fit cmdstan models
  model_name <- "hurdle_full"
  fit <- models[[model_name]]$sample(
    parallel_chains = pars$stan$n_cores,
    chains = pars$stan$n_chains,
    iter_warmup = floor(pars$stan$n_iter / 2),
    iter_sampling = pars$stan$n_iter,
    adapt_delta = pars$stan$adapt_delta,
    max_treedepth = pars$stan$max_treedepth,
    data = stan_data,
    init = init,
    step_size = 0.01
  )
  est <- calc_key_estimates(fit, data$people)
  expect_true(identical(names(est), c("b_mean", "l_mean", "l_male",
    "l_female", "l_male_effect")))
  expect_true(length(est$b_mean) == 1)
  expect_true(length(est$l_mean) == pars$sim$N_ind)

  model_name <- "hurdle_full"
  fit <- models[[model_name]]$sample(
    parallel_chains = pars$stan$n_cores,
    chains = pars$stan$n_chains,
    iter_warmup = floor(pars$stan$n_iter / 2),
    iter_sampling = pars$stan$n_iter,
    adapt_delta = pars$stan$adapt_delta,
    max_treedepth = pars$stan$max_treedepth,
    data = stan_data,
    init = init,
    step_size = 0.01
  )
  est <- calc_key_estimates(fit, data$people)
  expect_true(identical(names(est), c("b_mean", "l_mean", "l_male",
    "l_female", "l_male_effect")))
  expect_true(length(est$b_mean) == 1)
  expect_true(length(est$l_mean) == pars$sim$N_ind)

  model_name <- "lognormal_full"
  fit <- models[[model_name]]$sample(
    parallel_chains = pars$stan$n_cores,
    chains = pars$stan$n_chains,
    iter_warmup = floor(pars$stan$n_iter / 2),
    iter_sampling = pars$stan$n_iter,
    adapt_delta = pars$stan$adapt_delta,
    max_treedepth = pars$stan$max_treedepth,
    data = stan_data,
    init = init,
    step_size = 0.01
  )
  est <- calc_key_estimates(fit, data$people)
  expect_true(identical(names(est),
    c(
      "t0_mean",        "t0_male",       "t0_female",
      "t0_male_effect", "d_mean",        "d_male",
      "d_female",       "d_male_effect", "k_mean",
      "k_male",         "k_female",      "k_male_effect"
    )
  ))
  expect_true(length(est$t0_mean) == pars$sim$N_ind)
  expect_true(length(est$k_mean) == pars$sim$N_ind)
  expect_true(length(est$d_mean) == pars$sim$N_ind)

  model_name <- "lognormal_full"
  fit <- models[[model_name]]$sample(
    parallel_chains = pars$stan$n_cores,
    chains = pars$stan$n_chains,
    iter_warmup = floor(pars$stan$n_iter / 2),
    iter_sampling = pars$stan$n_iter,
    adapt_delta = pars$stan$adapt_delta,
    max_treedepth = pars$stan$max_treedepth,
    data = stan_data,
    init = init,
    step_size = 0.01
  )
  est <- calc_key_estimates(fit, data$people)
  expect_true(identical(names(est),
    c(
      "t0_mean",        "t0_male",       "t0_female",
      "t0_male_effect", "d_mean",        "d_male",
      "d_female",       "d_male_effect", "k_mean",
      "k_male",         "k_female",      "k_male_effect"
    )
  ))
  expect_true(length(est$t0_mean) == pars$sim$N_ind)
  expect_true(length(est$k_mean) == pars$sim$N_ind)
  expect_true(length(est$d_mean) == pars$sim$N_ind)

})


test_that("run_experiment works on well-formed parameter lists", {

  pars <- list(
    seed = 1,
    sim = list(
      longitudinal =  TRUE,
      N_ind        =    10,
      l_pop_mu     =     0,
      l_pop_sigma  =  0.18,
      l_ind_sigma  =   0.5,
      l_sex_sigma  =   1.0,
      b_mu         =   1.5, # b = 1 / v, v is a scale term in decades
      b_sigma      =   0.5,
      d_pop_mu     = -0.07, # on the log scale
      d_pop_sigma  =  0.01,
      d_ind_sigma  =  0.04,
      d_sex_sigma  =  0.01,
      s            =   0.1
    ),
    stan = list(
      # prior model parameters
      l_pop_mu          = 0,
      l_pop_sigma       = 0.18,
      l_ind_sigma_rate  = 1 / 0.5,
      l_sex_sigma_rate  = 1 / 1.0,
      b_mu              = 1.0,
      b_sigma           = 0.5,
      d_pop_mu          = -0.7,
      d_pop_sigma       = 0.01,
      d_ind_sigma_rate  = 1 / 0.04,
      d_sex_sigma_rate  = 1 / 0.01,
      s_rate            = 1 / 0.1,
      # sampling parameters
      n_chains = 3,
      n_cores = 3,
      n_iter = 200,
      adapt_delta = 0.8,
      max_treedepth = 15
    )
  )

  # check when fit_models = FALSE
  expect_error(capture.output(x <- run_experiment(pars, fit_models = FALSE)), NA)
  expect_true(setequal(names(x), c("seed", "sim", "stan", "start_time",
  "stop_time", "run_time", "data", "estimates", "summaries", "diagnostics")))
  # note fitted_models is removed even if we dont fit
  expect_true(length(x$estimates) == 0)
  expect_true(length(x$summaries) == 0)
  expect_true(length(x$diagnostics) == 0)

  # check when fit_models = TRUE
  expect_error(capture.output(x <- run_experiment(pars, fit_models = TRUE)), NA)
  expect_true(setequal(names(x), c("seed", "sim", "stan", "start_time",
    "stop_time", "run_time", "data", "estimates", "summaries", "diagnostics")))
  expect_true(length(x$estimates) == length(models))
  expect_true(length(x$summaries) == length(models))
  expect_true(length(x$diagnostics) == length(models))

})
