
rename_experiments <- function(path = "experiments") {

  files <- list.files(path, pattern="*.R$", full.names = TRUE, recursive = TRUE)

  rename_counter <- 0
  unchanged_counter <- 0

  for (i in 1:length(files)) {
    hash <- digest(files[i], algo = "sha1", file = TRUE)
    hash_name <- paste(substr(hash, 1, 11) ,".R", sep = "") 
    hash_name <- file.path(dirname(files[i]), hash_name)
    if(files[i] != hash_name){ 
      file.rename(files[i], hash_name)
      rename_counter <- rename_counter + 1
    } else {
      unchanged_counter <- unchanged_counter + 1
    }
  }

  print(paste0(rename_counter, " files have been renamed"))
  print(paste0(unchanged_counter, " files unchanged"))

}

inventory_experiments <- function(path = "experiments") {

  experiment_files <- list.files(path, pattern="*.R$", full.names = TRUE, recursive = TRUE)

  experiments <- list()
  sim <- list()
  stan <- list()
  for (i in 1:length(experiment_files)) {
    experiments[[i]] <- source(experiment_files[i])$value
    experiments[[i]]$file <- experiment_files[i]
    # save stan pars
    stan[[i]] <- experiments[[i]]$stan
    names(stan[[i]]) <- paste0("stan_", names(stan[[i]]))
    stan[[i]]$file <- experiment_files[i]
    # save sim pars
    sim[[i]] <- experiments[[i]]$sim
    names(sim[[i]]) <- paste0("sim_", names(sim[[i]]))
    sim[[i]]$file <- experiment_files[i]
    experiments[[i]]$stan <- NULL
    experiments[[i]]$sim <- NULL
  }

  experiments <- as.data.frame(bind_rows(experiments))
  stan <- as.data.frame(bind_rows(stan))
  sim <- as.data.frame(bind_rows(sim))

  experiments <- left_join(experiments, sim, by = "file")
  experiments <- left_join(experiments, stan, by = "file")

  return(experiments)

}

prep_experiments <- function(files, run_system, exp_seed) {
  # for each experiment, load experiment parameters from source
  set.seed(exp_seed)
  out <- list()
  for (i in seq_along(files)) {
    this_experiment <- source(files[i])$value
    this_experiment$system <- run_system
    if (!"n_runs" %in% names(this_experiment)) this_experiment$n_runs <- 1 
    my_seeds <- sample(100:999, this_experiment$n_runs)
    for (j in 1:this_experiment$n_runs) {
      this_run <- this_experiment
      if (this_experiment$n_runs > 1 | is.na(this_experiment$seed)) {
        this_run$seed <- my_seeds[j]
      }
      this_run$name <- paste(gsub(".R", "", basename(files[i])), this_run$seed, sep = "-")
      out <- c(out, list(this_run))
    }
  }
  return(out)
}

run_experiments <- function(n_threads = 1, this_system = NA, fit_models = TRUE, seed = 1) {
  # wrapper for batches of experiment files to run at once
  if (is.na(this_system)) stop("please give the name of the computer you are using")
  if (!dir.exists("results")) stop("results folder must exist first!")

  # load the parameters of all pending experiments and prep them
  experiment_files <- list.files("experiments", full.names = TRUE, recursive = TRUE, pattern = "*.R$")
  experiments <- prep_experiments(experiment_files, this_system, exp_seed = seed)
  experiment_names <- unlist(lapply(experiments, function(z) z[["name"]]))
  experiment_hashes <- substr(experiment_names, 1, 11)

  # compare experiment list to existing results
  completed_files <- list.files("results", full.names = FALSE)
  completed_names <- gsub("\\.rds", "", basename(completed_files))
  completed_hashes <- substr(completed_names, 1, 11)
  drop <- which(experiment_hashes %in% completed_hashes)

  if (length(drop) > 0) experiments <- experiments[-drop] 
  if (length(experiments) == 0) stop("no experiments to run!")

  # run all pending experiments and save output
  result <- mclapply(
    1:length(experiments),
    function(i) {
      out <- run_experiment(experiments[[i]], fit_models = fit_models)
      outname <- paste0(experiments[[i]]$name, ".rds")
      saveRDS(out, file.path("results", outname))
      print(outname)
    },
    mc.cores = n_threads
  )

  return("experiments run")

}


write_results_reports <- function(path = "results") {

  files <- list.files(path, pattern="*.rds$", full.names = TRUE, recursive = TRUE)

  tabular <- vector("list", length(files))
  ppl <- vector("list", length(files))
  obs <- vector("list", length(files))
  pars <- vector("list", length(files))

  for (i in seq_along(files)) {
    res <- readRDS(files[i])
    dir_name <- file.path(".", path, res$name)
    dir_init(dir_name)

    # summarize fits and export
    tabular[[i]] <- summarize_experiment(res)
    tabular[[i]]$experiment <- res$name
    write.csv(tabular[[i]], file.path(dir_name, paste0(res$name, ".csv")), row.names = FALSE)

    # save ppl and obs from experiment
    ppl[[i]] <- res$data$people
    obs[[i]] <- res$data$obs

    ppl[[i]]$experiment <- res$name
    obs[[i]]$experiment <- res$name

    res$N_ind <- res$sim$N_ind
    res$longitudinal <- res$sim$longitudinal
    obs[[i]]$longitudinal <- res$longitudinal
    ppl[[i]]$longitudinal <- res$longitudinal
    pars[[i]] <- list.remove(res, c("data", "estimates", "diagnostics", "summaries", "sim", "stan"))

    write.csv(ppl[[i]], file.path(dir_name, "people.csv"), row.names = FALSE)
    write.csv(obs[[i]], file.path(dir_name, "observations.csv"), row.names = FALSE)
    write.csv(pars[[i]], file.path(dir_name, "parameters.csv"), row.names = FALSE)
    print(res$name)
 
  }

  out <- as.data.frame(bind_rows(tabular))
  write.csv(out, file.path(path, "model_summaries.csv"), row.names = FALSE)

  out <- as.data.frame(bind_rows(ppl))
  write.csv(out, file.path(path, "people.csv"), row.names = FALSE)

  out <- as.data.frame(bind_rows(obs))
  write.csv(out, file.path(path, "observations.csv"), row.names = FALSE)

  out <- as.data.frame(bind_rows(pars))
  write.csv(out, file.path(path, "parameters.csv"), row.names = FALSE)

}

