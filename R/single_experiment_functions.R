

calc_key_estimates <- function(fit, people) {

  out <- list()

  post_array <- fit$draws() # iterations x chains x variable
  post <- as.data.frame(as_draws_df(post_array)) # (iterations x chains) x variable

  # hurdle component

  if ("b" %in% colnames(post)) out$b_mean <- mean(post$b)

  if (any(grepl("^l\\[", colnames(post)))) {
    post_l <- as.matrix(post[, grepl("^l\\[", colnames(post))])
    out$l_mean <- as.numeric(apply(post_l, 2, mean))
    if (any(people$sex == 1)) {
      post_l_male <- as.matrix(post_l[, people$sex == 1])
      out$l_male <- as.numeric(apply(post_l_male, 1, mean))
    }
    if (any(people$sex == 2)) {
      post_l_female <- as.matrix(post_l[, people$sex == 2])
      out$l_female <- as.numeric(apply(post_l_female, 1, mean))
    }
    if (all(c("l_female", "l_male") %in% names(out))) {
      out$l_male_effect <- as.numeric(out$l_male - out$l_female)
    }
  }

  # lognormal component

  if (any(grepl("^t0\\[", colnames(post)))) {
    post_t0 <- as.matrix(post[, grepl("^t0\\[", colnames(post))])
    out$t0_mean <- as.numeric(apply(post_t0, 2, mean))
    if (any(people$sex == 1)) {
      post_t0_male <- as.matrix(post_t0[, people$sex == 1])
      out$t0_male <- as.numeric(apply(post_t0_male, 1, mean))
    }
    if (any(people$sex == 2)) {
      post_t0_female <- as.matrix(post_t0[, people$sex == 2])
      out$t0_female <- as.numeric(apply(post_t0_female, 1, mean))
    }
    if (all(c("t0_female", "t0_male") %in% names(out))) {
      out$t0_male_effect <- as.numeric(out$t0_male - out$t0_female)
    }
  }

  if (any(grepl("^d\\[", colnames(post)))) {
    post_d <- as.matrix(post[, grepl("^d\\[", colnames(post))])
    out$d_mean <- as.numeric(apply(post_d, 2, mean))
    if (any(people$sex == 1)) {
      post_d_male <- as.matrix(post_d[, people$sex == 1])
      out$d_male <- as.numeric(apply(post_d_male, 1, mean))
    }
    if (any(people$sex == 2)) {
      post_d_female <- as.matrix(post_d[, people$sex == 2])
      out$d_female <- as.numeric(apply(post_d_female, 1, mean))
    }
    if (all(c("d_female", "d_male") %in% names(out))) {
      out$d_male_effect <- as.numeric(out$d_male - out$d_female)
    }
  }

  if (any(grepl("^k\\[", colnames(post)))) {
    post_k <- as.matrix(post[, grepl("^k\\[", colnames(post))])
    out$k_mean <- as.numeric(apply(post_k, 2, mean))
    if (any(people$sex == 1)) {
      post_k_male <- as.matrix(post_k[, people$sex == 1])
      out$k_male <- as.numeric(apply(post_k_male, 1, mean))
    }
    if (any(people$sex == 2)) {
      post_k_female <- as.matrix(post_k[, people$sex == 2])
      out$k_female <- as.numeric(apply(post_k_female, 1, mean))
    }
    if (all(c("k_female", "k_male") %in% names(out))) {
      out$k_male_effect <- as.numeric(out$k_male - out$k_female)
    }
  }

  return(out)

}


run_experiment <- function(pars, fit_models = TRUE) {

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
  if (fit_models) {
    for (i in seq_along(models)) {
      model_name <- names(models)[i]
      fitted_models[[model_name]] <- models[[model_name]]$sample(
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
      diagnostics[[model_name]] <- fitted_models[[model_name]]$cmdstan_diagnose()$stdout
      summaries[[model_name]] <- fitted_models[[model_name]]$cmdstan_summary()$stdout
      estimates[[model_name]] <- calc_key_estimates(fitted_models[[model_name]], data$people)
    }
  }

  ######### wrap up and save experiment output #########

  pars$stop_time <- Sys.time()
  pars$run_time <- round(difftime(pars$stop_time, pars$start_time), 1)

  results <- c(
    pars,
    list(
      data = data,
      estimates = estimates,
      summaries = summaries,
      diagnostics = diagnostics
    )
  )

  return(results)

}
