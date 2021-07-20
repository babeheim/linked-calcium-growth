
calc_lm_lognormal_plus1 <- function(res) {

  out <- list()

  obs <- res$data$obs
  people <- res$data$people
  obs$sex <- people$sex[match(obs$pid, people$pid)]

  obs$male <- as.numeric(obs$sex == 1)
  obs$lncac_plus1 <- log(obs$cac + 1)

  true_b <- res$sim$b
  true_t0 <- res$sim$t0 # one t0 per person; t0 is in decades
  true_t0_male_avg <- mean(true_t0[people$sex == 1])
  true_t0_female_avg <- mean(true_t0[people$sex == 2])
  true_t0_male_effect <- true_t0_male_avg - true_t0_female_avg

  true_k <- res$sim$k # one k per person; k is on the natural (non-log) scale
  true_k_male_effect <- mean(true_k[people$sex == 1]) - mean(true_k[people$sex == 2])

  true_d <- res$sim$d # one d per person; d is in years
  true_d_male_effect <- mean(true_d[people$sex == 1]) - mean(true_d[people$sex == 2])

  out$model <- "lm_lognormal_plus1"
  dm <- obs # all observations get used because its log(cac + 1)

  has_both_sexes <- var(dm$male) > 0.1
  if (var(dm$age_su) > 0 & has_both_sexes & nrow(dm) >= 20) {
    m_ln <- lm(lncac_plus1 ~ male*age_su, data = dm)
    post_ln <- extract.samples(m_ln)
    a_post <- ifelse(
      people$sex == 1,
      post_ln[, "Intercept"] + post_ln[, "male"],
      post_ln[, "Intercept"]
    )
    b_post <- ifelse(
      people$sex == 1,
      post_ln[, "age_su"] + post_ln[, "male:age_su"],
      post_ln[, "age_su"]
    )
    t0_post <- (-a_post / b_post)
    out$t0_mean <- mean(t0_post)
    out$t0_lb <- HPDI(t0_post)[1]
    out$t0_ub <- HPDI(t0_post)[2]
    out$t0_mean_bias <- out$t0_mean - mean(true_t0)
    a_female <- post_ln[, "Intercept"]
    b_female <- post_ln[, "age_su"]
    a_male <- post_ln[, "Intercept"] + post_ln[, "male"]
    b_male <- post_ln[, "age_su"] + post_ln[, "male:age_su"]
    t0_male_effect <- -(a_male / b_male) - (-a_female / b_female)
    out$t0_male_effect_mean <- mean(t0_male_effect)
    out$t0_male_effect_lb <- HPDI(t0_male_effect)[1]
    out$t0_male_effect_ub <- HPDI(t0_male_effect)[2]
    out$t0_male_effect_bias <- out$t0_male_effect_mean - true_t0_male_effect
    k_post <- ifelse(
      people$sex == 1,
      post_ln[, "age_su"] + post_ln[, "male:age_su"],
      post_ln[, "age_su"]
    )
    out$k_mean <- mean(k_post)
    out$k_lb <- HPDI(k_post)[1]
    out$k_ub <- HPDI(k_post)[2]
    out$k_mean_bias <- out$k_mean - mean(true_k)
    k_male_effect <- post_ln[, "male:age_su"]
    out$k_male_effect_mean <- mean(k_male_effect)
    out$k_male_effect_lb <- HPDI(k_male_effect)[1]
    out$k_male_effect_ub <- HPDI(k_male_effect)[2]
    out$k_male_effect_bias <- out$k_male_effect_mean - true_k_male_effect
    d_post <- log(2) / k_post
    out$d_mean <- mean(d_post)
    out$d_lb <- HPDI(d_post)[1]
    out$d_ub <- HPDI(d_post)[2]
    out$d_mean_bias <- out$d_mean - mean(true_d)
    d_male_effect <- post_ln[, "male:age_su"] / (b_female * b_male)
    out$d_male_effect_mean <- mean(d_male_effect)
    out$d_male_effect_lb <- HPDI(d_male_effect)[1]
    out$d_male_effect_ub <- HPDI(d_male_effect)[2]
    out$d_male_effect_bias <- out$d_male_effect_mean - true_d_male_effect
  }
  if (var(dm$age_su) > 0 & !has_both_sexes & nrow(dm) >= 20) {
    m_ln <- lm(lncac_plus1 ~ age_su, data = dm)
    post_ln <- extract.samples(m_ln)
    a_post <- post_ln[, "Intercept"]
    b_post <- post_ln[, "age_su"]
    t0_post <- (-a_post / b_post)
    out$t0_mean <- mean(t0_post)
    out$t0_lb <- HPDI(t0_post)[1]
    out$t0_ub <- HPDI(t0_post)[2]
    out$t0_mean_bias <- out$t0_mean - mean(true_t0)
    k_post <- b_post
    out$k_mean <- mean(k_post)
    out$k_lb <- HPDI(k_post)[1]
    out$k_ub <- HPDI(k_post)[2]
    out$k_mean_bias <- out$k_mean - mean(true_k)
    d_post <- log(2) / k_post
    out$d_mean <- mean(d_post)
    out$d_lb <- HPDI(d_post)[1]
    out$d_ub <- HPDI(d_post)[2]
    out$d_mean_bias <- out$d_mean - mean(true_d)
  }

  return(out)

}

calc_lm_lognormal <- function(res) {

  out <- list()

  obs <- res$data$obs
  people <- res$data$people
  obs$sex <- people$sex[match(obs$pid, people$pid)]

  obs$male <- as.numeric(obs$sex == 1)
  obs$lncac <- log(obs$cac)

  true_b <- res$sim$b
  true_t0 <- res$sim$t0 # one t0 per person; t0 is in decades
  true_t0_male_avg <- mean(true_t0[people$sex == 1])
  true_t0_female_avg <- mean(true_t0[people$sex == 2])
  true_t0_male_effect <- true_t0_male_avg - true_t0_female_avg

  true_k <- res$sim$k # one k per person; k is on the natural (non-log) scale
  true_k_male_effect <- mean(true_k[people$sex == 1]) - mean(true_k[people$sex == 2])

  true_d <- res$sim$d # one d per person; d is in years
  true_d_male_effect <- mean(true_d[people$sex == 1]) - mean(true_d[people$sex == 2])

  out$model <- "lm_lognormal"
  dm <- obs[which(obs$cac > 0),]

  has_both_sexes <- var(dm$male) > 0.1
  if (var(dm$age_su) > 0 & has_both_sexes & nrow(dm) >= 20) {
    m_ln <- lm(lncac ~ male*age_su, data = dm)
    post_ln <- extract.samples(m_ln)
    a_post <- ifelse(
      people$sex == 1,
      post_ln[, "Intercept"] + post_ln[, "male"],
      post_ln[, "Intercept"]
    )
    b_post <- ifelse(
      people$sex == 1,
      post_ln[, "age_su"] + post_ln[, "male:age_su"],
      post_ln[, "age_su"]
    )
    t0_post <- (-a_post / b_post)
    out$t0_mean <- mean(t0_post)
    out$t0_lb <- HPDI(t0_post)[1]
    out$t0_ub <- HPDI(t0_post)[2]
    out$t0_mean_bias <- out$t0_mean - mean(true_t0)
    a_female <- post_ln[, "Intercept"]
    b_female <- post_ln[, "age_su"]
    a_male <- post_ln[, "Intercept"] + post_ln[, "male"]
    b_male <- post_ln[, "age_su"] + post_ln[, "male:age_su"]
    t0_male_effect <- -(a_male / b_male) - (-a_female / b_female)
    out$t0_male_effect_mean <- mean(t0_male_effect)
    out$t0_male_effect_lb <- HPDI(t0_male_effect)[1]
    out$t0_male_effect_ub <- HPDI(t0_male_effect)[2]
    out$t0_male_effect_bias <- out$t0_male_effect_mean - true_t0_male_effect
    k_post <- ifelse(
      people$sex == 1,
      post_ln[, "age_su"] + post_ln[, "male:age_su"],
      post_ln[, "age_su"]
    )
    out$k_mean <- mean(k_post)
    out$k_lb <- HPDI(k_post)[1]
    out$k_ub <- HPDI(k_post)[2]
    out$k_mean_bias <- out$k_mean - mean(true_k)
    k_male_effect <- post_ln[, "male:age_su"]
    out$k_male_effect_mean <- mean(k_male_effect)
    out$k_male_effect_lb <- HPDI(k_male_effect)[1]
    out$k_male_effect_ub <- HPDI(k_male_effect)[2]
    out$k_male_effect_bias <- out$k_male_effect_mean - true_k_male_effect
    d_post <- log(2) / k_post
    out$d_mean <- mean(d_post)
    out$d_lb <- HPDI(d_post)[1]
    out$d_ub <- HPDI(d_post)[2]
    out$d_mean_bias <- out$d_mean - mean(true_d)
    d_male_effect <- post_ln[, "male:age_su"] / (b_female * b_male)
    out$d_male_effect_mean <- mean(d_male_effect)
    out$d_male_effect_lb <- HPDI(d_male_effect)[1]
    out$d_male_effect_ub <- HPDI(d_male_effect)[2]
    out$d_male_effect_bias <- out$d_male_effect_mean - true_d_male_effect
  }
  if (var(dm$age_su) > 0 & !has_both_sexes & nrow(dm) >= 20) {
    m_ln <- lm(lncac ~ age_su, data = dm)
    post_ln <- extract.samples(m_ln)
    a_post <- post_ln[, "Intercept"]
    b_post <- post_ln[, "age_su"]
    t0_post <- (-a_post / b_post)
    out$t0_mean <- mean(t0_post)
    out$t0_lb <- HPDI(t0_post)[1]
    out$t0_ub <- HPDI(t0_post)[2]
    out$t0_mean_bias <- out$t0_mean - mean(true_t0)
    k_post <- b_post
    out$k_mean <- mean(k_post)
    out$k_lb <- HPDI(k_post)[1]
    out$k_ub <- HPDI(k_post)[2]
    out$k_mean_bias <- out$k_mean - mean(true_k)
    d_post <- log(2) / k_post
    out$d_mean <- mean(d_post)
    out$d_lb <- HPDI(d_post)[1]
    out$d_ub <- HPDI(d_post)[2]
    out$d_mean_bias <- out$d_mean - mean(true_d)
  }

  return(out)

}

calc_glm_hurdle <- function(res) {

  out <- list()

  obs <- res$data$obs
  people <- res$data$people
  obs$sex <- people$sex[match(obs$pid, people$pid)]

  obs$male <- as.numeric(obs$sex == 1)
  obs$has_cac <- as.numeric(obs$cac > 0)

  true_b <- res$sim$b
  true_t0 <- res$sim$t0 # one t0 per person; t0 is in decades
  true_t0_male_avg <- mean(true_t0[people$sex == 1])
  true_t0_female_avg <- mean(true_t0[people$sex == 2])
  true_t0_male_effect <- true_t0_male_avg - true_t0_female_avg

  out$model <- "glm_hurdle"
  dm <- obs

  has_cac_var <- all(c(0, 1) %in% dm$has_cac)
  has_age_var <- var(dm$age_su) > 0
  has_both_sexes <- var(dm$male) > 0.1
  if (has_cac_var & has_age_var & has_both_sexes & nrow(dm) >= 20) {
    m_h <- glm(has_cac ~ male + age_su, data = dm, family = "binomial")
    post_h <- extract.samples(m_h)
    a_post <- ifelse(people$sex == 1,
      post_h[, "Intercept"] + post_h[, "male"],
      post_h[, "Intercept"]
    )
    b_post <- post_h[, "age_su"]
    # note: both men and women are assumed to have the same logistic b
    l_post <- -(a_post / b_post)
    out$l_mean <- mean(l_post)
    out$l_lb <- HPDI(l_post)[1]
    out$l_ub <- HPDI(l_post)[2]
    out$l_mean_bias <- out$l_mean - mean(true_t0)
    l_male_effect <- (-1 / b_post) * (post_h[, "male"])
    out$l_male_effect_mean <- mean(l_male_effect)
    out$l_male_effect_lb <- HPDI(l_male_effect)[1]
    out$l_male_effect_ub <- HPDI(l_male_effect)[2]
    out$l_male_effect_bias <- out$l_male_effect_mean - true_t0_male_effect
  }
  if (has_cac_var & has_age_var & !has_both_sexes & nrow(dm) >= 20) {
    m_h <- glm(has_cac ~ age_su, data = dm, family = "binomial")
    post_h <- extract.samples(m_h)
    a_post <- post_h[, "Intercept"]
    b_post <- post_h[, "age_su"]
    l_post <- -(a_post / b_post)
    out$l_mean <- mean(l_post)
    out$l_lb <- HPDI(l_post)[1]
    out$l_ub <- HPDI(l_post)[2]
    out$l_mean_bias <- out$l_mean - mean(true_t0)
  }

  return(out)

}



summarize_experiment <- function(res) {

  # output is one-row-per-model

  # pull out simulation data

  people <- res$data$people

  # pull out true parameter values
  # (need these here to calculate model estimate biases)

  true_b <- res$sim$b
  true_t0 <- res$sim$t0 # one entry per person
  true_t0_mean <- mean(res$sim$t0)
  true_t0_male_avg <- mean(true_t0[people$sex == 1])
  true_t0_female_avg <- mean(true_t0[people$sex == 2])
  true_t0_male_effect <- true_t0_male_avg - true_t0_female_avg

  true_k <- res$sim$k # one entry per person
  true_k_mean <- mean(true_k)
  true_k_male_effect <- mean(true_k[people$sex == 1]) - mean(true_k[people$sex == 2])

  true_d <- res$sim$d # one entry per person
  true_d_mean <- mean(true_d)
  true_d_male_effect <- mean(true_d[people$sex == 1]) - mean(true_d[people$sex == 2])

  # evaluate each model and make a report

  model_names <- names(res$summaries)
  n_models <- length(model_names)

  report <- vector("list", n_models + 3) # plus 3 for the three lm / glm models to include

  for (i in seq_along(model_names)) {

    report[[i]]$model <- model_names[i]

    # assess HMC sampling diagnostics
    x <- res$diagnostics[[model_names[i]]]
    report[[i]]$rhat_good <- grepl("Split R-hat values satisfactory all parameters.", x)
    report[[i]]$ess_good <- grepl("Effective sample size satisfactory.", x)
    report[[i]]$energy_good <- grepl("E-BFMI satisfactory for all transitions.", x)
    report[[i]]$divergence_good <- grepl("No divergent transitions found.", x)

    # extract model timings
    x <- res$summaries[[model_names[i]]]
    pattern <- "Sampling took.*\\n"
    m <- gregexpr(pattern, x, perl = TRUE)
    report[[i]]$sampling_time <- regmatches(x, m)[[1]]
    report[[i]]$sampling_time <- gsub("\\n$", "", report[[i]]$sampling_time)
    report[[i]]$sampling_time_min <- extract_sampling_time(report[[i]]$sampling_time)

    # calculate key estimates
    x <- res$estimates[[model_names[i]]]
    if ("b_mean" %in% names(x)) report[[i]]$b_mean <- x$b_mean

    if ("l_mean" %in% names(x)) {
      report[[i]]$l_mean <- mean(x$l_mean)
      if (length(x$l_mean) > 1) {
        report[[i]]$l_lb <- HPDI(x$l_mean)[1]
        report[[i]]$l_ub <- HPDI(x$l_mean)[2]
      }
      expect_true(length(x$l_mean) == length(true_t0))
      report[[i]]$l_mean_bias <- mean(x$l_mean - true_t0) # need to confirm this is the right dimensions...
    }
    if ("l_male_effect" %in% names(x)) {
      report[[i]]$l_male_effect_mean <- mean(x$l_male_effect)
      report[[i]]$l_male_effect_lb <- HPDI(x$l_male_effect)[1]
      report[[i]]$l_male_effect_ub <- HPDI(x$l_male_effect)[2]
      report[[i]]$l_male_effect_bias <- mean(x$l_male_effect - true_t0_male_effect)
    }
    if ("t0_mean" %in% names(x)) {
      report[[i]]$t0_mean <- mean(x$t0_mean)
      if (length(x$t0_mean) > 1) {
        report[[i]]$t0_lb <- HPDI(x$t0_mean)[1]
        report[[i]]$t0_ub <- HPDI(x$t0_mean)[2]
      }
      expect_true(length(x$t0_mean) == length(true_t0))
      report[[i]]$t0_mean_bias <- mean(x$t0_mean - true_t0)
    }
    if ("t0_male_effect" %in% names(x)) {
      report[[i]]$t0_male_effect_mean <- mean(x$t0_male_effect)
      report[[i]]$t0_male_effect_lb <- HPDI(x$t0_male_effect)[1]
      report[[i]]$t0_male_effect_ub <- HPDI(x$t0_male_effect)[2]
      report[[i]]$t0_male_effect_bias <- mean(x$t0_male_effect - true_t0_male_effect)
    }

    if ("k_mean" %in% names(x)) {
      report[[i]]$k_mean <- mean(x$k_mean)
      if (length(x$k_mean) > 1) {
        report[[i]]$k_lb <- HPDI(x$k_mean)[1]
        report[[i]]$k_ub <- HPDI(x$k_mean)[2]
      }
      expect_true(length(x$k_mean) == length(true_t0))
      report[[i]]$k_mean_bias <- mean(x$k_mean - true_k)
    }
    if ("k_male_effect" %in% names(x)) {
      report[[i]]$k_male_effect_mean <- mean(x$k_male_effect)
      report[[i]]$k_male_effect_lb <- HPDI(x$k_male_effect)[1]
      report[[i]]$k_male_effect_ub <- HPDI(x$k_male_effect)[2]
      report[[i]]$k_male_effect_bias <- mean(x$k_male_effect - true_k_male_effect)
    }

    if ("d_mean" %in% names(x)) {
      report[[i]]$d_mean <- mean(x$d_mean)
      if (length(x$d_mean) > 1) {
        report[[i]]$d_lb <- HPDI(x$d_mean)[1]
        report[[i]]$d_ub <- HPDI(x$d_mean)[2]
      }
      expect_true(length(x$d_mean) == length(true_t0))
      report[[i]]$d_mean_bias <- mean(x$d_mean - true_d)
    }
    if ("d_male_effect" %in% names(x)) {
      report[[i]]$d_male_effect_mean <- mean(x$d_male_effect)
      report[[i]]$d_male_effect_lb <- HPDI(x$d_male_effect)[1]
      report[[i]]$d_male_effect_ub <- HPDI(x$d_male_effect)[2]
      report[[i]]$d_male_effect_bias <- mean(x$d_male_effect - true_d_male_effect)
    }

  }

  # calculate two more report entries: glm estimates

  i <- length(report) - 2
  report[[i]] <- calc_glm_hurdle(res)

  i <- length(report) - 1
  report[[i]] <- calc_lm_lognormal(res)

  i <- length(report)
  report[[i]] <- calc_lm_lognormal_plus1(res)

  out <- as.data.frame(bind_rows(report))

  out$N_ind <- res$sim$N_ind
  out$long <- res$sim$longitudinal

  # pull out stan settings

  out$n_chains <- res$stan$n_chains
  out$n_cores <- res$stan$n_cores
  out$n_iter <- res$stan$n_iter
  out$adapt_delta <- res$stan$adapt_delta
  out$max_treedepth <- res$stan$max_treedepth

  out$true_b <- true_b
  out$true_t0_mean <- true_t0_mean
  out$true_t0_male_effect <- true_t0_male_effect
  out$true_k_mean <- true_k_mean
  out$true_k_male_effect <- true_k_male_effect
  out$true_d_mean <- true_d_mean
  out$true_d_male_effect <- true_d_male_effect

  return(out)

}
