
rm(list = ls())

source("project_support.R")

dir_init(results_folder)
run_experiments(path = results_folder, n_threads = n_threads,
  this_system = this_system, seed = my_seed)
write_results_reports(results_folder)
