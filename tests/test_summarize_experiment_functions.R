
test_that("linear models and summarize_experiment also work", {

    pars <- list(
      seed = 1,
      sim = list(
        longitudinal =  TRUE,
        N_ind        =     5,
        l_pop_mu     =     0,
        l_pop_sigma  =  0.18,
        l_ind_sigma  =   0.5,
        l_sex_sigma  =   1.0,
        b_mu         =   1.5, # b = 1 / v, v is a scale term in decades
        b_sigma      =   0.5,
        d_pop_mu     =  -0.7, # on the log scale
        d_pop_sigma  =   0.01,
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

    expect_error(capture.output(x <- run_experiment(pars, fit_models = TRUE)), NA)

    # check run_experiment output
    expect_true(setequal(names(x), c("seed", "sim", "stan", "start_time", "stop_time", "run_time", "data", "estimates", "summaries", "diagnostics")))
    expect_true(length(x$estimates) == length(models))
    expect_true(length(x$summaries) == length(models))
    expect_true(length(x$diagnostics) == length(models))

    # check post-experiment glm and lm model fitting
    expect_silent(model_results <- calc_glm_hurdle(x))
    expect_true(abs(model_results$l_mean_bias - (-0.015)) < 0.01)
    expect_true(abs(model_results$l_male_effect_bias - (-0.078)) < 0.01)
    expect_silent(model_results <- calc_lm_lognormal(x))
    expect_true(abs(model_results$t0_mean_bias - (-1.014)) < 0.01)
    expect_true(abs(model_results$t0_male_effect_bias - (-1.266)) < 0.01)
    expect_true(abs(model_results$k_mean_bias - (-0.019)) < 0.01)
    expect_true(abs(model_results$k_male_effect_bias - (-0.023)) < 0.01)
    expect_silent(model_results <- calc_lm_lognormal_plus1(x))
    expect_true(abs(model_results$t0_mean_bias - (-3.915)) < 0.01)
    expect_true(abs(model_results$t0_male_effect_bias - (-2.581)) < 0.01)
    expect_true(abs(model_results$k_mean_bias - (-0.905)) < 0.01)
    expect_true(abs(model_results$k_male_effect_bias - (-0.555)) < 0.01)

    # check summarize_experiment
    expect_silent(res <- summarize_experiment(x))
    expect_true(!any(is.null(res)))
    expect_true(!any(is.na(res$true_b)) & length(unique(res$true_b)) == 1)
    expect_true(!any(is.na(res$true_t0_mean)) & length(unique(res$true_t0_mean)) == 1)
    expect_true(!any(is.na(res$true_t0_male_effect)) & length(unique(res$true_t0_male_effect)) == 1)
    expect_true(!any(is.na(res$true_k_mean)) & length(unique(res$true_k_mean)) == 1)
    expect_true(!any(is.na(res$true_k_male_effect)) & length(unique(res$true_k_male_effect)) == 1)
    expect_true(abs(mean(res$t0_mean_bias, na.rm = TRUE) - (-1.338)) < 0.05)
    expect_true(abs(mean(res$t0_male_effect_bias, na.rm = TRUE) - (-1.240)) < 0.01)
    expect_true(abs(mean(res$k_mean_bias, na.rm = TRUE) - (-0.187)) < 0.05)
    expect_true(abs(mean(res$k_male_effect_bias, na.rm = TRUE) - (-0.118)) < 0.01)

})
