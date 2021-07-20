
rm(list = ls())

source("project_support.R")

# define my set of populations and study designs

calcs <- list()

sex_effect <- c("l off, k off", "l on, k off", "l on, k on")
ind_variation <- c("none", "lots")
priors <- c("perfect", "regularizing")

N_ind <- c(1, 5, 20)

long_exp <- expand.grid(
  longitudinal = TRUE,
  priors = priors,
  ind_variation = ind_variation,
  sex_effect = sex_effect,
  sim_N_ind = N_ind
)

N_ind <- c(20, 50, 100, 200, 500, 1000)

cross_exp <- expand.grid(
  longitudinal = FALSE,
  priors = priors,
  ind_variation = ind_variation,
  sex_effect = sex_effect,
  sim_N_ind = N_ind
)

all_exp <- bind_rows(long_exp, cross_exp) # 114 designs x 5 is 570 studies!

# initialize remaining variables

all_exp$seed <- "NA"
all_exp$n_runs <- 5

# set the expected t0 values

all_exp$sim_l_pop_mu <- 0
all_exp$sim_l_pop_sigma <- 0.2

all_exp$sim_l_sex_sigma <- case_when(
  all_exp$sex_effect == "l off, k off" ~ 0.0001,
  all_exp$sex_effect == "l on, k off" ~ 0.4,
  all_exp$sex_effect == "l on, k on" ~ 0.4
)
all_exp$sim_l_ind_sigma <- case_when(
  all_exp$ind_variation == "none" ~ 0.0001,
  all_exp$ind_variation == "lots" ~ 0.4
)

# the two sexes are IID normal with mean 0 and sigma sex_effect, but
# their difference is normal with sigma sex_effect * sqrt(2)
# since the variance is s1^2 + s2^2

calcs$studySexDiffOnset <- 0.4 * sqrt(2) # the RMS difference between sexes
calcs$studyIndDiffOnset <- 0.4 # the SD between individuals

# set background variation

all_exp$sim_b_mu <- 2
all_exp$sim_b_sigma <- 0.2
all_exp$sim_s <- 0.1

# set the doubling times

all_exp$sim_d_pop_mu <- (-0.7)
all_exp$sim_d_pop_sigma <- 0.2

all_exp$sim_d_sex_sigma <- case_when(
  all_exp$sex_effect == "l off, k off" ~ 0.0001,
  all_exp$sex_effect == "l on, k off" ~ 0.0001,
  all_exp$sex_effect == "l on, k on" ~ 0.1
)
all_exp$sim_d_ind_sigma <- case_when(
  all_exp$ind_variation == "none" ~ 0.0001,
  all_exp$ind_variation == "lots" ~ 0.1
)

calcs$studySexDiffDoubling <- 0.1 * sqrt(2) # the RMS difference between sexes
calcs$studyIndDiffDoubling <- 0.1 # the SD between individuals

# set the samping parameters

all_exp$stan_n_chains <- 4
all_exp$stan_n_cores <- 4
all_exp$stan_adapt_delta <- 0.9
all_exp$stan_max_treedepth <- 15
all_exp$stan_n_iter <- NA

tar <- which(all_exp$longitudinal)
all_exp$stan_n_iter[tar] <- 3000
all_exp$stan_n_iter[-tar] <- 5000

# set the priors

all_exp$stan_l_pop_mu <- NA
all_exp$stan_b_mu <- NA
all_exp$stan_b_sigma <- NA
all_exp$stan_l_pop_sigma <- NA
all_exp$stan_l_ind_sigma_rate <- NA
all_exp$stan_l_sex_sigma_rate <- NA
all_exp$stan_d_pop_mu <- NA
all_exp$stan_d_pop_sigma <- NA
all_exp$stan_d_ind_sigma_rate <- NA
all_exp$stan_d_sex_sigma_rate <- NA

# these priors are always perfect
all_exp$stan_l_pop_mu <- all_exp$sim_l_pop_mu
all_exp$stan_d_pop_mu <- all_exp$sim_d_pop_mu
all_exp$stan_b_mu <- all_exp$sim_b_mu
all_exp$stan_b_sigma <- all_exp$sim_b_sigma
all_exp$stan_s_rate <- 1/all_exp$sim_s

tar <- which(all_exp$priors == "perfect")
all_exp$stan_l_pop_sigma[tar] <- all_exp$sim_l_pop_sigma[tar]
all_exp$stan_l_ind_sigma_rate[tar] <- 1/all_exp$sim_l_ind_sigma[tar]
all_exp$stan_l_sex_sigma_rate[tar] <- 1/all_exp$sim_l_sex_sigma[tar]
all_exp$stan_d_pop_sigma[tar] <- all_exp$sim_d_pop_sigma[tar]
all_exp$stan_d_ind_sigma_rate[tar] <- 1/all_exp$sim_d_ind_sigma[tar]
all_exp$stan_d_sex_sigma_rate[tar] <- 1/all_exp$sim_d_sex_sigma[tar]

tar <- which(all_exp$priors == "regularizing")
all_exp$stan_l_pop_sigma[tar] <- 1.0 # 1 decade, vs true sigma of +/- 2 years (0.2)
all_exp$stan_l_ind_sigma_rate[tar] <- 1/2 # 2 decades; vs true sigma of 0.0001 or 0.4 (or +/- 4 years)
all_exp$stan_l_sex_sigma_rate[tar] <- 1/2 # 2 decades; vs true sigma of 0.0001 or 0.4 (or +/- 4 years)
all_exp$stan_d_pop_sigma[tar] <- 0.3 # vs true value of 0.2
all_exp$stan_d_ind_sigma_rate[tar] <- 1/0.2 # vs true sigma of 0.0001 or 0.1 in sim
all_exp$stan_d_sex_sigma_rate[tar] <- 1/0.2 # vs true sigma of 0.0001 or 0.1 in sim

expect_true(all(colSums(is.na(all_exp[, -1])) == 0))

# export calculations
calcs$studySexDiffOnset <- sprintf("%1.1f", 10 * calcs$studySexDiffOnset)
calcs$studyIndDiffOnset <- sprintf("%1.1f", 10 * calcs$studyIndDiffOnset)
calcs$studySexDiffDoubling <- sprintf("%1.1f", 10 * calcs$studySexDiffDoubling)
calcs$studyIndDiffDoubling <- sprintf("%1.1f", 10 * calcs$studyIndDiffDoubling)
for (i in 1:length(calcs)) calcs[[i]] <- paste0(calcs[[i]], " years")
writeLines(prep_latex_variables(calcs), "calcsExperimentSetup.tex")


# write each parameter list to file

dir_init("./experiments")

template <- "list(
  sex_effect = $SEX_EFFECT,
  priors = $PRIORS,
  ind_variation = $IND_VARIATION,
  sampling = \"normal\",
  seed = $SEED,
  n_runs = $N_RUNS,
  sim = list(
    longitudinal = $LONGITUDINAL,
    N_ind        =   $SIM_N_IND,
    l_pop_mu     =   $SIM_L_POP_MU,
    l_pop_sigma  =   $SIM_L_POP_SIGMA,
    l_ind_sigma  =   $SIM_L_IND_SIGMA,
    l_sex_sigma  =   $SIM_L_SEX_SIGMA,
    b_mu         =   $SIM_B_MU,
    b_sigma      =   $SIM_B_SIGMA,
    d_pop_mu     =   $SIM_D_POP_MU,
    d_pop_sigma  =   $SIM_D_POP_SIGMA,
    d_ind_sigma  =   $SIM_D_IND_SIGMA,
    d_sex_sigma  =   $SIM_D_SEX_SIGMA,
    s            =   $SIM_S
  ),
  stan = list(
    # prior model parameters
    l_pop_mu          = $STAN_L_POP_MU,
    l_pop_sigma       = $STAN_L_POP_SIGMA,
    l_ind_sigma_rate  = $STAN_L_IND_SIGMA_RATE,
    l_sex_sigma_rate  = $STAN_L_SEX_SIGMA_RATE,
    b_mu              = $STAN_B_MU,
    b_sigma           = $STAN_B_SIGMA,
    d_pop_mu          = $STAN_D_POP_MU,
    d_pop_sigma       = $STAN_D_POP_SIGMA,
    d_ind_sigma_rate  = $STAN_D_IND_SIGMA_RATE,
    d_sex_sigma_rate  = $STAN_D_SEX_SIGMA_RATE,
    s_rate            = $STAN_S_RATE,
    # sampling parameters
    n_chains = $STAN_N_CHAINS,
    n_cores = $STAN_N_CORES,
    n_iter = $STAN_N_ITER,
    adapt_delta = $STAN_ADAPT_DELTA,
    max_treedepth = $STAN_MAX_TREEDEPTH
  )
)
"

all_exp$priors <- as.character(all_exp$priors)
all_exp$ind_variation <- as.character(all_exp$ind_variation)
all_exp$sex_effect <- as.character(all_exp$sex_effect)

for (i in 1:nrow(all_exp)) {
  my_filename <- paste0("experiments/script_", sprintf("%04d", i), ".R")
  exp_script <- template
  for (j in 1:ncol(all_exp)) {
    my_variable <- paste0("\\$", toupper(colnames(all_exp)[j]))
    my_value <- all_exp[i, j]
    if (is.character(my_value)) my_value <- paste0("\\\"", my_value, "\\\"")
    if (grepl(my_variable, exp_script)) {
      exp_script <- gsub(my_variable, my_value, exp_script)
      if (is.na(exp_script)) print(my_variable)
    }
  }
  writeLines(exp_script, my_filename)
}

rename_experiments(path = "experiments")
