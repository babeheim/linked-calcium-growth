
trnorm <- function(n, mean, sd, lower = 0, upper = NA){
   lower_quantile <- 0
   if (!is.na(lower)) lower_quantile <- pnorm(lower, mean = mean, sd = sd)
   upper_quantile <- 1
   if (!is.na(upper)) upper_quantile <- pnorm(upper, mean = mean, sd = sd)
   out <- qnorm(runif(n, min = lower_quantile, max = upper_quantile), mean = mean, sd = sd)
   return(out)
}

sim_cac <- function(pars) {

  if (!hasName(pars, "seed")) stop("seed must be NA or specified")
  if (!hasName(pars, "longitudinal")) stop("longitudinal flag must be specified")
  if (is.na(pars$seed)) pars$seed <- sample(1:1000, 1)

  set.seed(pars$seed)

  pars$l_pop <- rnorm(1, pars$l_pop_mu, pars$l_pop_sigma)

  pars$k_pop <- rnorm(1, pars$k_pop_mu, pars$k_pop_sigma)

  if (pars$b_mu <= 0) stop("b_mu must be positive")
  pars$b <- trnorm(1, pars$b_mu, pars$b_sigma)
  pars$v <- 1/pars$b

  pars$l_sex <- rnorm(2, 0, pars$l_sex_sigma)
  pars$k_sex <- rnorm(2, 0, pars$k_sex_sigma)

  pars$l_ind <- rep(NA, pars$N_ind)
  pars$l <- rep(NA, pars$N_ind)
  pars$t0 <- rep(NA, pars$N_ind)
  pars$k_ind <- rep(NA, pars$N_ind)
  pars$k <- rep(NA, pars$N_ind)

  sim <- list()

  for (i in 1:pars$N_ind) {
    # initialize person
    dat <- list(
      pid = i,
      sex = sample(1:2, 1), # where 1 = male, 2 = female
      year_of_birth = sample(1910:1970, 1),
      age = 1:100
    )
    dat$age_su = (dat$age - 50)/10
    pars$l_ind[i] <- rnorm(1, 0, pars$l_ind_sigma)
    pars$l[i] <- pars$l_pop +
      pars$l_sex[dat$sex] +
      pars$l_ind[i]
    # assumption: no one can have an l before being born
    if (pars$l[i] < (-5)) pars$l[i] <- (-5)
    pars$t0[i] <- rlogis(1, pars$l[i], pars$v)
    # assumption: no one can have a t0 before being born
    if (pars$t0[i] < (-5)) pars$t0[i] <- (-5)
    pars$k_ind[i] <- rnorm(1, 0, pars$k_ind_sigma)
    pars$k[i] <- exp(pars$k_pop +
      pars$k_sex[dat$sex] +
      pars$k_ind[i])
    # assumption: no one can have negative growth (log link)
    # calculate cac over lifespan
    dat$cac <- rep(0, length(dat$age))
    dat$growth_rate <- rep(NA, length(dat$age))
    growing <- which(dat$age_su >= pars$t0[i])
    if (length(growing) > 0) {
      dat$growth_rate[growing] <- rep(pars$k[i], length(growing))
      dat$cac[growing] <- exp(dat$growth_rate[growing] * (dat$age_su[growing] - pars$t0[i]) +
        rnorm(length(growing), 0, pars$s))
    }
    # record results
    dat$cac_start_age <- pars$t0[i] * 10 + 50
    dat$obs_cohort <- dat$year_of_birth + dat$age
    sim[[i]] <- dat
  }

  people <- list()
  for (i in 1:pars$N_ind) people[[i]] <- sim[[i]][c("pid", "sex", "year_of_birth", "cac_start_age")]
  people %>% bind_rows() %>% as.data.frame() -> people

  obs <- list()
  for (i in 1:pars$N_ind) {
    obs[[i]] <- sim[[i]][c("age", "obs_cohort", "cac", "growth_rate")]
    obs[[i]]$pid <- rep(sim[[i]]$pid, length(sim[[i]]$age))
  }
  obs %>% bind_rows() %>% as.data.frame() %>% select(pid, obs_cohort, age, growth_rate, cac) -> obs

  if (!pars$longitudinal) {
    # define an observation year or years
    people$obs_cohort <- 2000
    obs$pid_obs_cohort <- people$obs_cohort[match(obs$pid, people$pid)]
    # subset only to scans in that person's observation year
    obs <- obs[which(obs$obs_cohort == obs$pid_obs_cohort),]
  }

  out <- list(
    pars = pars,
    people = people,
    obs = obs
  )
  return(out)

}

