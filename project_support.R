
##### project libraries #####

library(dplyr)
library(parallel)
library(posterior)
library(bayesplot)
library(tictoc)
library(rethinking) # https://github.com/rmcelreath/rethinking
library(jsonlite)
library(testthat)
library(knitr)
library(loo)
library(brms)
library(cmdstanr)
library(digest)
library(rlist)
library(tikzDevice)

##### global parameters #####

this_system <- "local-imac"
n_threads <- 2
my_seed <- 1
results_folder <- "./results"

set.seed(my_seed)
rstan_options(javascript = FALSE)
options(warnPartialMatchDollar=TRUE)

##### compile model programs #####

models <- list()
stan_files <- list.files("./stan", pattern = "*.stan$", full.names = TRUE)
for (i in 1:length(stan_files)) {
  model_name <- gsub(".stan", "", basename(stan_files[i]))
  models[[model_name]] <- cmdstan_model(stan_files[i])
}

##### project functions #####

source("R/misc_functions.R")
source("R/simulation_functions.R")
source("R/single_experiment_functions.R")
source("R/summarize_experiment_functions.R")
source("R/multi_experiment_functions.R")
