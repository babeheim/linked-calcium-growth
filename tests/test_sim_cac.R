

test_that("sim_cac won't accept invalid input", {

  pars <- list(
    longitudinal =  TRUE,
    seed         =    NA,
    N_ind        =    -3, # N_ind must be a positive integer
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =    NA,
    N_ind        =   1.2, # N_ind must be a positive integer
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =    NA,
    N_ind        =     0, # N_ind must be a positive integer
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
                          # seed cannot be missing
    N_ind        =     5,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   937,
    N_ind        =     5,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =  -1.5, # b cannot be negative
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   937,
    N_ind        =     5,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5
    # bunch of pars missing
  )

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   937,
    N_ind        =     5,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
#    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )
# one par missing

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   937,
    N_ind        =     5,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
#    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )
# one par missing

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   937,
    N_ind        =     5,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
#    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )
# one par missing

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   937,
    N_ind        =     5,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
#    d_sex_sigma  =  0.01,
    s            =   0.1
  )
# one par missing

  expect_error(sim_cac(pars))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   818,
    N_ind        =   100,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =    -5, # hypergrowth, doubles every 2 days
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_error(sim <- sim_cac(pars)) # does not allow infinite cac growth!

})

test_that("sim_cac accepts valid input", {

  pars <- list(
    longitudinal =  TRUE,
    seed         =     2,
    N_ind        =     1,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.01,
    b_mu         =   1.5,
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(abs(mean(sim$obs$cac) - 173.91) < 0.01)
  expect_true(abs(sim$pars$t0 - (-0.51)) < 0.01)
  expect_true(abs(sim$people$cac_start_age - (44.86)) < 0.01)

# everyone is identical

  pars <- list(
    longitudinal =  TRUE,
    seed         =     2,
    N_ind        =   100,
    l_pop_mu     =   0.0,
    l_pop_sigma  =   0.0,
    l_ind_sigma  =   0.0,
    l_sex_sigma  =   0.0,
    b_mu         =   100,
    b_sigma      =   0.1,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.00,
    d_ind_sigma  =  0.00,
    d_sex_sigma  =  0.00,
    s            =  0.01
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(abs(mean(sim$obs$cac) - 82.80) < 0.01)
  expect_true(all(abs(sim$pars$t0 - (-0.003)) < 0.1))
  expect_true(all(abs(sim$people$cac_start_age - (50)) < 1))

  pars <- list(
    longitudinal =  TRUE,
    seed         =   265,
    N_ind        =  1000,
    l_pop_mu     =   0.0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =  0.16,
    l_sex_sigma  =  0.32,
    b_mu         =   1.5,
    b_sigma      =   0.5, # no growth, doubles every 1480 years
    d_pop_mu     =     5,
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(all(sim$obs$cac < 2))

  pars <- list(
    longitudinal = FALSE,
    seed         =   442,
    N_ind        =  1000,
    l_pop_mu     =   1.0, # everyone starts growing at age 60
    l_pop_sigma  =   0.0,
    l_ind_sigma  =   0.0,
    l_sex_sigma  =   0.0,
    b_mu         =   100,
    b_sigma      =  0.00,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.00,
    d_ind_sigma  =  0.00,
    d_sex_sigma  =  0.00,
    s            =   0.1
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(all(abs(sim$pars$t0 - 1) < 0.1))

  pars <- list(
    longitudinal = FALSE,
    seed         =   442,
    N_ind        =  1000,
    l_pop_mu     =  20.0, # no one starts growing cac at all
    l_pop_sigma  =   0.0,
    l_ind_sigma  =   0.0,
    l_sex_sigma  =   0.0,
    b_mu         =   100,
    b_sigma      =  0.00,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.00,
    d_ind_sigma  =  0.00,
    d_sex_sigma  =  0.00,
    s            =   0.1
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(all(sim$obs$cac == 0))

  pars <- list(
    longitudinal = FALSE,
    seed         =   442,
    N_ind        =  1000,
    l_pop_mu     =  -8.0, # everyone starts growing at birth
    l_pop_sigma  =   0.0,
    l_ind_sigma  =   0.0,
    l_sex_sigma  =   0.0,
    b_mu         =   100,
    b_sigma      =  0.00,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.00,
    d_ind_sigma  =  0.00,
    d_sex_sigma  =  0.00,
    s            =   0.1
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(all(sim$pars$t0 - (-5) < 0.1))

  # turn up individual variation via l_ind_sigma

  pars <- list(
    longitudinal = FALSE,
    seed         =   442,
    N_ind        =  1000,
    l_pop_mu     =   0.0,
    l_pop_sigma  = 0.000,
    l_ind_sigma  =     2, # +/- 2 decades of SD
    l_sex_sigma  = 0.000,
    b_mu         =   100,
    b_sigma      =   0.0,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.00,
    d_ind_sigma  =  0.00,
    d_sex_sigma  =  0.00,
    s            =   0.1
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(abs(mean(sim$pars$t0) - (-0.084)) < 0.01)
  expect_true(abs(sd(sim$pars$t0) - 2.085) < 0.01)
  expect_true(abs(diff(range(sim$pars$t0)) - 10.860) < 0.01) # range is ~109 years

  # turn up individual variation via b

  pars <- list(
    longitudinal = FALSE,
    seed         =   442,
    N_ind        =  1000,
    l_pop_mu     =   0.0,
    l_pop_sigma  = 0.000,
    l_ind_sigma  =     0, # +/- 2 decades of SD
    l_sex_sigma  = 0.000,
    b_mu         =   0.5,
    b_sigma      =   0.0,
    d_pop_mu     =  -0.7,
    d_pop_sigma  =  0.00,
    d_ind_sigma  =  0.00,
    d_sex_sigma  =  0.00,
    s            =   0.1
  )

  expect_silent(sim <- sim_cac(pars))
  expect_true(abs(mean(sim$pars$t0) - (0.100)) < 0.01)
  expect_true(abs(sd(sim$pars$t0) - 3.315) < 0.01)
  expect_true(abs(diff(range(sim$pars$t0)) - 25.652) < 0.01) # range is ~250 years

})
