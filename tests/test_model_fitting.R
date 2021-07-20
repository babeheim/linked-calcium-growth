

test_that("all models sample without error", {

  pars <- list()

  pars$sim <- list(
    longitudinal =  TRUE,
    seed         =     1,
    N_ind        =     7,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =   0.2,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     = -0.07, # on the log scale
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  sim <- sim_cac(pars$sim)

  obs <- sim$obs
  people <- sim$people

  obs$age_su <- (obs$age - 50) / 10
  obs$sex <- people$sex[match(obs$pid, people$pid)]

  pm <- people
  dm <- obs

  pars$stan <- list(
    l_pop_mu          = 0,
    l_pop_sigma       = 0.18,
    l_ind_sigma_rate  = 1 / 0.2,
    l_sex_sigma_rate  = 1 / 0.32,
    b_mu              = 1.5,
    b_sigma           = 0.5,
    d_pop_mu          = -0.7,
    d_pop_sigma       = 0.01,
    d_ind_sigma_rate  = 1 / 0.04,
    d_sex_sigma_rate  = 1 / 0.01,
    s_rate            = 1 / 0.1
  )

  stan_data <- c(
    pars$stan,
    list(
      N_ind = nrow(pm),
      sex = as.array(pm$sex),
      N_obs = nrow(dm),
      pid = as.array(dm$pid),
      age_su = as.array(dm$age_su),
      y = as.array(dm$cac)
    )
  )

  init <- list(sim$pars, sim$pars, sim$pars)

  for (i in seq_along(models)) {
    model_name <- names(models)[i]
    expect_error(
      capture.output(
        models[[model_name]]$sample(
          seed = 1, # set for consistency
          refresh = 0, # disables output
          parallel_chains = 3,
          chains = 3,
          iter_warmup = 50,
          iter_sampling = 100,
          adapt_delta = 0.8,
          max_treedepth = 15,
          data = stan_data,
          init = init,
          step_size = 0.01,
          show_messages = FALSE
        )
      ),
      NA # the NA here means expect_NO_error
    )
  }

})




test_that("linked hurdle lognormal fits longitudinal accurately", {

  pars <- list()

  pars$sim <- list(
    longitudinal =  TRUE,
    seed         =     1,
    N_ind        =     7,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =   0.2,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     = -0.07, # on the log scale
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  sim <- sim_cac(pars$sim)

  obs <- sim$obs
  people <- sim$people

  obs$age_su <- (obs$age - 50) / 10
  obs$sex <- people$sex[match(obs$pid, people$pid)]

  pm <- people
  dm <- obs

  pars$stan <- list(
    l_pop_mu          = 0,
    l_pop_sigma       = 0.18,
    l_ind_sigma_rate  = 1 / 0.2,
    l_sex_sigma_rate  = 1 / 0.32,
    b_mu              = 1.5,
    b_sigma           = 0.5,
    d_pop_mu          = -0.7,
    d_pop_sigma       = 0.01,
    d_ind_sigma_rate  = 1 / 0.04,
    d_sex_sigma_rate  = 1 / 0.01,
    s_rate            = 1 / 0.1
  )

  stan_data <- c(
    pars$stan,
    list(
      N_ind = nrow(pm),
      sex = as.array(pm$sex),
      N_obs = nrow(dm),
      pid = as.array(dm$pid),
      age_su = as.array(dm$age_su),
      y = as.array(dm$cac)
    )
  )

  init <- list(sim$pars, sim$pars, sim$pars, sim$pars)

  model_name <- "linked_hurdle_lognormal_full"

  fit <- models[[model_name]]$sample(
    seed = 1, # set for consistency
    refresh = 0, # disables output
    parallel_chains = 4,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    data = stan_data,
    init = init,
    step_size = 0.01,
    show_messages = FALSE
  )

  diagnostics <- fit$cmdstan_diagnose()$stdout
  summaries <- fit$cmdstan_summary()$stdout
  estimates <- calc_key_estimates(fit, people)

  rhat_good <- grepl("Split R-hat values satisfactory all parameters.", diagnostics)
  expect_true(rhat_good)
  ess_good <- grepl("Effective sample size satisfactory.", diagnostics)
  expect_true(ess_good)
  energy_good <- grepl("E-BFMI satisfactory for all transitions.", diagnostics)
  expect_true(energy_good)

  true_b <- sim$pars$b
  true_t0 <- sim$pars$t0 # one entry per person
  true_t0_mean <- mean(sim$pars$t0)
  true_t0_male_avg <- mean(true_t0[people$sex == 1])
  true_t0_female_avg <- mean(true_t0[people$sex == 2])
  true_t0_male_effect <- true_t0_male_avg - true_t0_female_avg

  true_k <- sim$pars$k # one entry per person
  true_k_mean <- mean(true_k)
  true_k_male_effect <- mean(true_k[people$sex == 1]) - mean(true_k[people$sex == 2])

  true_d <- sim$pars$d # one entry per person
  true_d_mean <- mean(true_d)
  true_d_male_effect <- mean(true_d[people$sex == 1]) - mean(true_d[people$sex == 2])

  # calculate key estimates
  report <- list()
  x <- estimates

  if ("b_mean" %in% names(x)) report$b_mean <- x$b_mean

  if ("l_mean" %in% names(x)) {
    report$l_mean <- mean(x$l_mean)
    if (length(x$l_mean) > 1) {
      report$l_lb <- HPDI(x$l_mean)[1]
      report$l_ub <- HPDI(x$l_mean)[2]
    }
    expect_true(length(x$l_mean) == length(true_t0))
    report$l_mean_bias <- mean(x$l_mean - true_t0)
  }

  if ("l_male_effect" %in% names(x)) {
    report$l_male_effect_mean <- mean(x$l_male_effect)
    report$l_male_effect_lb <- HPDI(x$l_male_effect)[1]
    report$l_male_effect_ub <- HPDI(x$l_male_effect)[2]
    report$l_male_effect_bias <- mean(x$l_male_effect - true_t0_male_effect)
  }

  if ("t0_mean" %in% names(x)) {
    report$t0_mean <- mean(x$t0_mean)
    if (length(x$t0_mean) > 1) {
      report$t0_lb <- HPDI(x$t0_mean)[1]
      report$t0_ub <- HPDI(x$t0_mean)[2]
    }
    expect_true(length(x$t0_mean) == length(true_t0))
    report$t0_mean_bias <- mean(x$t0_mean - true_t0)
  }

  if ("t0_male_effect" %in% names(x)) {
    report$t0_male_effect_mean <- mean(x$t0_male_effect)
    report$t0_male_effect_lb <- HPDI(x$t0_male_effect)[1]
    report$t0_male_effect_ub <- HPDI(x$t0_male_effect)[2]
    report$t0_male_effect_bias <- mean(x$t0_male_effect - true_t0_male_effect)
  }

  if ("k_mean" %in% names(x)) {
    report$k_mean <- mean(x$k_mean)
    if (length(x$k_mean) > 1) {
      report$k_lb <- HPDI(x$k_mean)[1]
      report$k_ub <- HPDI(x$k_mean)[2]
    }
    expect_true(length(x$k_mean) == length(true_t0))
    report$k_mean_bias <- mean(x$k_mean - true_k)
  }

  if ("k_male_effect" %in% names(x)) {
    report$k_male_effect_mean <- mean(x$k_male_effect)
    report$k_male_effect_lb <- HPDI(x$k_male_effect)[1]
    report$k_male_effect_ub <- HPDI(x$k_male_effect)[2]
    report$k_male_effect_bias <- mean(x$k_male_effect - true_k_male_effect)
  }

  if ("d_mean" %in% names(x)) {
    report$d_mean <- mean(x$d_mean)
    if (length(x$d_mean) > 1) {
      report$d_lb <- HPDI(x$d_mean)[1]
      report$d_ub <- HPDI(x$d_mean)[2]
    }
    expect_true(length(x$d_mean) == length(true_t0))
    report$d_mean_bias <- mean(x$d_mean - true_k)
  }

  if ("d_male_effect" %in% names(x)) {
    report$d_male_effect_mean <- mean(x$d_male_effect)
    report$d_male_effect_lb <- HPDI(x$d_male_effect)[1]
    report$d_male_effect_ub <- HPDI(x$d_male_effect)[2]
    report$d_male_effect_bias <- mean(x$d_male_effect - true_d_male_effect)
  }

  expect_true(abs(report$l_mean_bias - (-0.004)) < 0.005)
  expect_true(abs(report$l_male_effect_bias - (0.007)) < 0.005)
  expect_true(abs(report$t0_mean_bias - (-0.009)) < 0.005)
  expect_true(abs(report$t0_male_effect_bias - (-0.022)) < 0.005)
  expect_true(abs(report$k_mean_bias - (0.015)) < 0.005)
  expect_true(abs(report$k_male_effect_bias - 0.036) < 0.005)
  expect_true(abs(report$d_mean_bias - (0.191)) < 0.005)
  expect_true(abs(report$d_male_effect_bias - (-0.005)) < 0.005)

})

test_that("linked hurdle lognormal fits cross-sectional accurately", {

  pars <- list()

  pars$sim <- list(
    longitudinal = FALSE,
    seed         =     1,
    N_ind        =   100,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =   0.2,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     = -0.07, # on the log scale
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  sim <- sim_cac(pars$sim)

  obs <- sim$obs
  people <- sim$people

  obs$age_su <- (obs$age - 50) / 10
  obs$sex <- people$sex[match(obs$pid, people$pid)]

  pm <- people
  dm <- obs

  pars$stan <- list(
    l_pop_mu          = 0,
    l_pop_sigma       = 0.18,
    l_ind_sigma_rate  = 1 / 0.2,
    l_sex_sigma_rate  = 1 / 0.32,
    b_mu              = 1.5,
    b_sigma           = 0.5,
    d_pop_mu          = -0.7,
    d_pop_sigma       = 0.01,
    d_ind_sigma_rate  = 1 / 0.04,
    d_sex_sigma_rate  = 1 / 0.01,
    s_rate            = 1 / 0.1
  )

  stan_data <- c(
    pars$stan,
    list(
      N_ind = nrow(pm),
      sex = as.array(pm$sex),
      N_obs = nrow(dm),
      pid = as.array(dm$pid),
      age_su = as.array(dm$age_su),
      y = as.array(dm$cac)
    )
  )

  init <- list(sim$pars, sim$pars, sim$pars, sim$pars)

  model_name <- "linked_hurdle_lognormal_full"

  fit <- models[[model_name]]$sample(
    seed = 1, # set for consistency
    refresh = 0, # disables output
    parallel_chains = 4,
    chains = 4,
    iter_warmup = 6000,
    iter_sampling = 12000,
    adapt_delta = 0.9,
    max_treedepth = 15,
    data = stan_data,
    init = init,
    step_size = 0.01,
    show_messages = FALSE
  )

  diagnostics <- fit$cmdstan_diagnose()$stdout
  summaries <- fit$cmdstan_summary()$stdout
  estimates <- calc_key_estimates(fit, people)

  rhat_good <- grepl("Split R-hat values satisfactory all parameters.", diagnostics)
  expect_true(rhat_good)
  ess_good <- grepl("Effective sample size satisfactory.", diagnostics)
  expect_true(ess_good)

  true_b <- sim$pars$b
  true_t0 <- sim$pars$t0 # one entry per person
  true_t0_mean <- mean(sim$pars$t0)
  true_t0_male_avg <- mean(true_t0[people$sex == 1])
  true_t0_female_avg <- mean(true_t0[people$sex == 2])
  true_t0_male_effect <- true_t0_male_avg - true_t0_female_avg

  true_k <- sim$pars$k # one entry per person
  true_k_mean <- mean(true_k)
  true_k_male_effect <- mean(true_k[people$sex == 1]) - mean(true_k[people$sex == 2])

  true_d <- sim$pars$d # one entry per person
  true_d_mean <- mean(true_d)
  true_d_male_effect <- mean(true_d[people$sex == 1]) - mean(true_d[people$sex == 2])

  # calculate key estimates
  report <- list()
  x <- estimates

  if ("b_mean" %in% names(x)) report$b_mean <- x$b_mean

  if ("l_mean" %in% names(x)) {
    report$l_mean <- mean(x$l_mean)
    if (length(x$l_mean) > 1) {
      report$l_lb <- HPDI(x$l_mean)[1]
      report$l_ub <- HPDI(x$l_mean)[2]
    }
    expect_true(length(x$l_mean) == length(true_t0))
    report$l_mean_bias <- mean(x$l_mean - true_t0)
  }

  if ("l_male_effect" %in% names(x)) {
    report$l_male_effect_mean <- mean(x$l_male_effect)
    report$l_male_effect_lb <- HPDI(x$l_male_effect)[1]
    report$l_male_effect_ub <- HPDI(x$l_male_effect)[2]
    report$l_male_effect_bias <- mean(x$l_male_effect - true_t0_male_effect)
  }

  if ("t0_mean" %in% names(x)) {
    report$t0_mean <- mean(x$t0_mean)
    if (length(x$t0_mean) > 1) {
      report$t0_lb <- HPDI(x$t0_mean)[1]
      report$t0_ub <- HPDI(x$t0_mean)[2]
    }
    expect_true(length(x$t0_mean) == length(true_t0))
    report$t0_mean_bias <- mean(x$t0_mean - true_t0)
  }

  if ("t0_male_effect" %in% names(x)) {
    report$t0_male_effect_mean <- mean(x$t0_male_effect)
    report$t0_male_effect_lb <- HPDI(x$t0_male_effect)[1]
    report$t0_male_effect_ub <- HPDI(x$t0_male_effect)[2]
    report$t0_male_effect_bias <- mean(x$t0_male_effect - true_t0_male_effect)
  }

  if ("k_mean" %in% names(x)) {
    report$k_mean <- mean(x$k_mean)
    if (length(x$k_mean) > 1) {
      report$k_lb <- HPDI(x$k_mean)[1]
      report$k_ub <- HPDI(x$k_mean)[2]
    }
    expect_true(length(x$k_mean) == length(true_t0))
    report$k_mean_bias <- mean(x$k_mean - true_k)
  }

  if ("k_male_effect" %in% names(x)) {
    report$k_male_effect_mean <- mean(x$k_male_effect)
    report$k_male_effect_lb <- HPDI(x$k_male_effect)[1]
    report$k_male_effect_ub <- HPDI(x$k_male_effect)[2]
    report$k_male_effect_bias <- mean(x$k_male_effect - true_k_male_effect)
  }

  if ("d_mean" %in% names(x)) {
    report$d_mean <- mean(x$d_mean)
    if (length(x$d_mean) > 1) {
      report$d_lb <- HPDI(x$d_mean)[1]
      report$d_ub <- HPDI(x$d_mean)[2]
    }
    expect_true(length(x$d_mean) == length(true_t0))
    report$d_mean_bias <- mean(x$d_mean - true_k)
  }

  if ("d_male_effect" %in% names(x)) {
    report$d_male_effect_mean <- mean(x$d_male_effect)
    report$d_male_effect_lb <- HPDI(x$d_male_effect)[1]
    report$d_male_effect_ub <- HPDI(x$d_male_effect)[2]
    report$d_male_effect_bias <- mean(x$d_male_effect - true_d_male_effect)
  }

  expect_true(abs(report$l_mean_bias - (0.499)) < 0.005) # its off by 4.9 years
  expect_true(abs(report$l_male_effect_bias - (-0.163)) < 0.005)
  expect_true(abs(report$t0_mean_bias - (0.541)) < 0.005) # its off by 5.4 years
  expect_true(abs(report$t0_male_effect_bias - (-0.139)) < 0.005)
  expect_true(abs(report$k_mean_bias - (0.609)) < 0.005)
  expect_true(abs(report$k_male_effect_bias - (0.0043)) < 0.005)
  expect_true(abs(report$d_mean_bias - (-0.236)) < 0.005) # off by 2.3 years
  expect_true(abs(report$d_male_effect_bias - (-0.0058)) < 0.005)

})
