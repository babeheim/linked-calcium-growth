
rm(list = ls())

source("project_support.R")

calcs <- list()

# create tableExperimentSummary

para <- read.csv(file.path(results_folder, "parameters.csv"), stringsAsFactors = FALSE)
expect_true(nrow(para) == 540)

calcs$numStudies <- nrow(para)
calcs$numStudiesLong <- sum(para$longitudinal)
calcs$numStudiesCross <- sum(!para$longitudinal)
calcs$numStudiesPerPop <- para$n_runs[1]
calcs$numPops <- calcs$numStudies / calcs$numStudiesPerPop

para %>%
  filter(longitudinal) %>%
  summarize(
    N_exp = n(),
    ind_var_off = sum(ind_variation == "none"),
    ind_var = sum(ind_variation == "lots"),
    sex_off = sum(sex_effect == "l off, k off"),
    sex_t0 = sum(sex_effect != "l off, k off"),
    sex_k = sum(sex_effect == "l on, k on"),
    priors_normal = sum(priors == "regularizing"),
    priors_perfect = sum(priors == "perfect"),
    N_1 = sum(N_ind == 1),
    N_5 = sum(N_ind == 5),
    N_20 = sum(N_ind == 20),
    N_50 = sum(N_ind == 50),
    N_100 = sum(N_ind == 100),
    N_200 = sum(N_ind == 200),
    N_500 = sum(N_ind == 500),
    N_1000 = sum(N_ind == 1000)
  ) -> exp_long

para %>%
  filter(!longitudinal) %>%
  summarize(
    N_exp = n(),
    ind_var_off = sum(ind_variation == "off"),
    ind_var = sum(ind_variation == "lots"),
    sex_off = sum(sex_effect == "l off, k off"),
    sex_t0 = sum(sex_effect != "l off, k off"),
    sex_k = sum(sex_effect == "l on, k on"),
    priors_normal = sum(priors == "regularizing"),
    priors_perfect = sum(priors == "perfect"),
    N_1 = sum(N_ind == 1),
    N_5 = sum(N_ind == 5),
    N_20 = sum(N_ind == 20),
    N_50 = sum(N_ind == 50),
    N_100 = sum(N_ind == 100),
    N_200 = sum(N_ind == 200),
    N_500 = sum(N_ind == 500),
    N_1000 = sum(N_ind == 1000)
  ) -> exp_cross

out <- cbind(t(as.matrix(exp_long)), t(as.matrix(exp_cross)))

out <- cbind(rownames(out), out)

colnames(out) <- c("", "Longitudinal", "Cross-sectional")

out <- insert_row(out, 1, "Inter-Individual Variation")
out <- insert_row(out, 4, "Sex Effect")
out <- insert_row(out, 8, "Prior Accuracy")
out <- insert_row(out, 11, "Study Size")

out[,1] <- gsub("N_exp", "Total Studies", out[,1])
out[,1] <- gsub("ind_var_off", "\\\\hspace{3 mm} No Variation", out[,1])
out[,1] <- gsub("ind_var", "\\\\hspace{3 mm} $\\\\pm$ 3 years on $t_0$", out[,1])
out[,1] <- gsub("sex_off", "\\\\hspace{3 mm} No sex differences", out[,1])
out[,1] <- gsub("sex_t0", "\\\\hspace{3 mm} Sex difference in $t_0$", out[,1])
out[,1] <- gsub("sex_k", "\\\\hspace{3 mm} Sex difference in $t_0$, $d$", out[,1])
out[,1] <- gsub("priors_normal", "\\\\hspace{3 mm} Regularizing Priors", out[,1])
out[,1] <- gsub("priors_perfect", "\\\\hspace{3 mm} Simulation Priors", out[,1])
out[,1] <- gsub("N_1000", "\\\\hspace{3 mm} $N$ = 1000 Patients", out[,1])
out[,1] <- gsub("N_500", "\\\\hspace{3 mm} $N$ = 500 Patients", out[,1])
out[,1] <- gsub("N_200", "\\\\hspace{3 mm} $N$ = 200 Patients", out[,1])
out[,1] <- gsub("N_100", "\\\\hspace{3 mm} $N$ = 100 Patients", out[,1])
out[,1] <- gsub("N_50", "\\\\hspace{3 mm} $N$ = 50 Patients", out[,1])
out[,1] <- gsub("N_20", "\\\\hspace{3 mm} $N$ = 20 Patients", out[,1])
out[,1] <- gsub("N_5", "\\\\hspace{3 mm} $N$ = 5 Patients", out[,1])
out[,1] <- gsub("N_1", "\\\\hspace{3 mm} $N$ = 1 Patient", out[,1])

tex_table <- rbind(colnames(out), out)
writeLines(texttab(tex_table, hlines = 1), "tableExperimentSummary.tex")

print("created tableExperimentSummary")



# create tablePatientSummary

ppl <- read.csv(file.path(results_folder, "people.csv"), stringsAsFactors = FALSE)
expect_true(nrow(ppl) == 113760)

calcs$numPatients <- nrow(ppl)

obs <- read.csv(file.path(results_folder, "observations.csv"), stringsAsFactors = FALSE)
expect_true(nrow(obs) == 268200)

calcs$numObservations <- nrow(obs)

ppl$cac_start_age[ppl$cac_start_age > 100] <- 101

ppl$key <- paste(ppl$experiment, ppl$pid)
obs$key <- paste(obs$experiment, obs$pid)

expect_true(all(obs$key %in% ppl$key))
obs$sex <- ppl$sex[match(obs$key, ppl$key)]

ppl$n_obs <- table(obs$key)[ppl$key]

# males

ppl %>%
  filter(sex == 1) %>%
  summarize(
    n_ppl = n(),
    n_ppl_long = sum(longitudinal),
    n_ppl_cross = sum(!longitudinal),
    t0_avg = mean(cac_start_age),
    t0_sd = sd(cac_start_age)
  ) -> male_avg_by_ppl

obs %>%
  filter(sex == 1) %>%
  summarize(
    n_obs = n(),
    n_obs_long = sum(longitudinal),
    n_obs_cross = sum(!longitudinal),
    k_avg = mean(growth_rate, na.rm = TRUE),
    k_sd = sd(growth_rate, na.rm = TRUE),
    d_avg = mean(10 * log(2)/growth_rate, na.rm = TRUE),
    d_sd = sd(10 * log(2)/growth_rate, na.rm = TRUE),
    prop_0 = mean(cac == 0),
    prop_1_10 = mean(cac >= 1 & cac < 10),
    prop_10_100 = mean(cac >= 10 & cac < 100),
    prop_100_400 = mean(cac >= 100 & cac < 400),
    prop_400_1000 = mean(cac >= 400 & cac < 1000),
    prop_1000_plus = mean(cac >= 1000)
  ) -> male_avg_by_obs

male_avg <- bind_cols(male_avg_by_ppl, male_avg_by_obs)

# females

ppl %>%
  filter(sex == 2) %>%
  summarize(
    n_ppl = n(),
    n_ppl_long = sum(longitudinal),
    n_ppl_cross = sum(!longitudinal),
    t0_avg = mean(cac_start_age),
    t0_sd = sd(cac_start_age)
  ) -> female_avg_by_ppl

obs %>%
  filter(sex == 2) %>%
  summarize(
    n_obs = n(),
    n_obs_long = sum(longitudinal),
    n_obs_cross = sum(!longitudinal),
    k_avg = mean(growth_rate, na.rm = TRUE),
    k_sd = sd(growth_rate, na.rm = TRUE),
    d_avg = mean(10 * log(2)/growth_rate, na.rm = TRUE),
    d_sd = sd(10 * log(2)/growth_rate, na.rm = TRUE),
    prop_0 = mean(cac == 0),
    prop_1_10 = mean(cac >= 1 & cac < 10),
    prop_10_100 = mean(cac >= 10 & cac < 100),
    prop_100_400 = mean(cac >= 100 & cac < 400),
    prop_400_1000 = mean(cac >= 400 & cac < 1000),
    prop_1000_plus = mean(cac >= 1000)
  ) -> female_avg_by_obs

female_avg <- bind_cols(female_avg_by_ppl, female_avg_by_obs)

# total

ppl %>%
  summarize(
    n_ppl = n(),
    n_ppl_long = sum(longitudinal),
    n_ppl_cross = sum(!longitudinal),
    t0_avg = mean(cac_start_age),
    t0_sd = sd(cac_start_age)
  ) -> total_avg_by_ppl

obs %>%
  summarize(
    n_obs = n(),
    n_obs_long = sum(longitudinal),
    n_obs_cross = sum(!longitudinal),
    k_avg = mean(growth_rate, na.rm = TRUE),
    k_sd = sd(growth_rate, na.rm = TRUE),
    d_avg = mean(10 * log(2)/growth_rate, na.rm = TRUE),
    d_sd = sd(10 * log(2)/growth_rate, na.rm = TRUE),
    prop_0 = mean(cac == 0),
    prop_1_10 = mean(cac >= 1 & cac < 10),
    prop_10_100 = mean(cac >= 10 & cac < 100),
    prop_100_400 = mean(cac >= 100 & cac < 400),
    prop_400_1000 = mean(cac >= 400 & cac < 1000),
    prop_1000_plus = mean(cac >= 1000)
  ) -> total_avg_by_obs

total_avg <- bind_cols(total_avg_by_ppl, total_avg_by_obs)

out <- bind_rows(female_avg, male_avg, total_avg)
rownames(out) <- c("women", "men", "total")

out <- select(out,
  n_ppl, n_ppl_long, n_ppl_cross, n_obs, n_obs_long, n_obs_cross,
  t0_avg, t0_sd, k_avg, k_sd, d_avg, d_sd, prop_0, prop_1_10, prop_10_100,
  prop_100_400, prop_400_1000, prop_1000_plus
)

out$t0_avg <- sprintf("%1.1f", out$t0_avg)
out$t0_sd <- sprintf("%1.1f", out$t0_sd)
out$d_avg <- sprintf("%1.1f", out$d_avg)
out$d_sd <- sprintf("%1.1f", out$d_sd)
out$k_avg <- sprintf("%1.1f", out$k_avg)
out$k_sd <- sprintf("%1.1f", out$k_sd)
out$prop_0 <- sprintf("%1.1f", 100 * out$prop_0)
out$prop_1_10 <- sprintf("%1.1f", 100 * out$prop_1_10)
out$prop_10_100 <- sprintf("%1.1f", 100 * out$prop_10_100)
out$prop_100_400 <- sprintf("%1.1f", 100 * out$prop_100_400)
out$prop_400_1000 <- sprintf("%1.1f", 100 * out$prop_400_1000)
out$prop_1000_plus <- sprintf("%1.1f", 100 * out$prop_1000_plus)

out$n_ppl <- format(out$n_ppl, big.mark = ",")
out$n_ppl_long <- format(out$n_ppl_long, big.mark = ",")
out$n_ppl_cross <- format(out$n_ppl_cross, big.mark = ",")
out$n_obs <- format(out$n_obs, big.mark = ",")
out$n_obs_long <- format(out$n_obs_long, big.mark = ",")
out$n_obs_cross <- format(out$n_obs_cross, big.mark = ",")

out$t0_avg <- paste0(out$t0_avg, " (", out$t0_sd, ")")
out$k_avg <- paste0(out$k_avg, " (", out$k_sd, ")")
out$d_avg <- paste0(out$d_avg, " (", out$d_sd, ")")
out <- select(out, -t0_sd, -d_sd, -k_sd)

out <- t(out)

out <- cbind(rownames(out), out)
out[,1] <- gsub("n_ppl_cross", "\\\\hspace{3 mm} Cross-sectional", out[,1])
out[,1] <- gsub("n_ppl_long", "\\\\hspace{3 mm} Longitudinal", out[,1])
out[,1] <- gsub("n_ppl", "Number of Patients", out[,1])
out[,1] <- gsub("n_obs_cross", "\\\\hspace{3 mm} Cross-sectional", out[,1])
out[,1] <- gsub("n_obs_long", "\\\\hspace{3 mm} Longitudinal", out[,1])
out[,1] <- gsub("n_obs", "Number of Observations", out[,1])
out[,1] <- gsub("t0_avg", "Onset $t_0$ (years)", out[,1])
out[,1] <- gsub("k_avg", "Progression $k$", out[,1])
out[,1] <- gsub("d_avg", "Doubling time $d$ (years)", out[,1])
out[,1] <- gsub("prop_0", "\\\\hspace{3 mm} 0", out[,1])
out[,1] <- gsub("prop_1_10", "\\\\hspace{3 mm} 1 to 9", out[,1])
out[,1] <- gsub("prop_10_100", "\\\\hspace{3 mm} 10 to 99", out[,1])
out[,1] <- gsub("prop_100_400", "\\\\hspace{3 mm} 100 to 399", out[,1])
out[,1] <- gsub("prop_400_1000", "\\\\hspace{3 mm} 400 to 999", out[,1])
out[,1] <- gsub("prop_1000_plus", "\\\\hspace{3 mm} $\\\\geq$ 1000", out[,1])

out <- insert_row(out, 9, "CAC category (\\%)")

colnames(out) <- c("", "Women", "Men", "Total")
tex_table <- rbind(colnames(out), out)
writeLines(texttab(tex_table, hlines = 1), "tablePatientSummary.tex")

print("created tablePatientSummary")


# calculate tab:modelBias

d <- read.csv(file.path(results_folder, "model_summaries.csv"), stringsAsFactors = FALSE)
expect_true(nrow(d) == 3780)

calcs$numModels <- length(unique(d$model))

# check missingness patterns

expect_true(!any(is.na(d$experiment)))
expect_true(!any(is.na(d$N_ind)))

tar <- which(d$model %in% c("glm_hurdle", "hurdle_full", "hurdle_lognormal_full", "linked_hurdle_lognormal_full"))
expect_true(!any(is.na(d$l_mean[tar])))
expect_true(!all(!is.na(d$l_male_effect_bias[tar])))
# might be NA in studies with only men or only women -> na.rm = TRUE is okay

tar <- which(d$model %in% c("lm_lognormal") & d$N_ind == 20 & d$long == FALSE)
expect_true(all(is.na(d$t0_mean[tar])))
# all fail at N = 20, b/c we have to drop some obs in the lm_lognormal model
# -> drop these entries
d <- d[-tar,]

tar <- which(d$model %in% c("lm_lognormal", "lm_lognormal_plus1", "lognormal_full", "hurdle_lognormal_full", "linked_hurdle_lognormal_full"))
expect_true(!any(is.na(d$t0_mean[tar])))
expect_true(!all(!is.na(d$t0_male_effect_bias[tar])))
# might be NA in studies with only men or only women -> na.rm = TRUE is okay
expect_true(!any(is.na(d$k_mean[tar])))
expect_true(!all(!is.na(d$k_male_effect_bias[tar])))
# might be NA in studies with only men or only women -> na.rm = TRUE is okay
expect_true(!any(is.na(d$d_mean[tar])))
expect_true(!all(!is.na(d$d_male_effect_bias[tar])))
# might be NA in studies with only men or only women -> na.rm = TRUE is okay

tar <- which(d$model %in% c("glm_hurdle", "lm_lognormal", "lm_lognormal_plus1"))
expect_true(all(is.na(d$sampling_time_min[tar])))
# non-Stan models have no sampling time stats

# drop models that did not converge
drop <- which(d$rhat_good == FALSE)
expect_true(length(drop) == 701)
d <- d[-drop,]

# pathological fits in glm_hurdle and lm_lognormal; drop in any analysis
drop <- which(
  abs(d$l_mean_bias) > 1e1 |
  abs(d$l_male_effect_bias) > 2e1 |
  abs(d$t0_male_effect_bias) > 2e1
)
expect_true(length(drop) == 8)
d <- d[-drop,]

d %>%
  group_by(long, model) %>%
  summarize(
    rhat_good = round(mean(rhat_good), 2) * 100,
    ess_good = round(mean(ess_good), 2) * 100,
    energy_good = round(mean(energy_good), 2) * 100,
    divergence_good = round(mean(divergence_good), 2) * 100,
    sampling_time = round(mean(sampling_time_min), 1),
    l_mean_bias_avg = mean(l_mean_bias),
    l_mean_bias_se = se(l_mean_bias),
    l_male_effect_bias_avg = mean(l_male_effect_bias, na.rm = TRUE),
    l_male_effect_bias_se = se(l_male_effect_bias),
    t0_mean_bias_avg = mean(t0_mean_bias),
    t0_mean_bias_se = se(t0_mean_bias),
    t0_male_effect_bias_avg = mean(t0_male_effect_bias, na.rm = TRUE),
    t0_male_effect_bias_se = se(t0_male_effect_bias),
    k_mean_bias_avg = mean(k_mean_bias),
    k_mean_bias_se = se(k_mean_bias),
    k_male_effect_bias_avg = mean(k_male_effect_bias, na.rm = TRUE),
    k_male_effect_bias_se = se(k_male_effect_bias),
    d_mean_bias_avg = mean(d_mean_bias),
    d_mean_bias_se = se(d_mean_bias),
    d_male_effect_bias_avg = mean(d_male_effect_bias, na.rm = TRUE),
    d_male_effect_bias_se = se(d_male_effect_bias)
  ) %>% as.data.frame() -> model_design_res

model_design_res$l_mean_bias_p <- calc_twotail_p(
  model_design_res$l_mean_bias_avg /
  model_design_res$l_mean_bias_se
)

model_design_res$l_male_effect_bias_p <- calc_twotail_p(
  model_design_res$l_male_effect_bias_avg /
  model_design_res$l_male_effect_bias_se
)

model_design_res$t0_mean_bias_p <- calc_twotail_p(
  model_design_res$t0_mean_bias_avg /
  model_design_res$t0_mean_bias_se
)
model_design_res$t0_male_effect_bias_p <- calc_twotail_p(
  model_design_res$t0_male_effect_bias_avg /
  model_design_res$t0_male_effect_bias_se
)

model_design_res$k_mean_bias_p <- calc_twotail_p(
  model_design_res$k_mean_bias_avg /
  model_design_res$k_mean_bias_se
)

model_design_res$k_male_effect_bias_p <- calc_twotail_p(
  model_design_res$k_male_effect_bias_avg /
  model_design_res$k_male_effect_bias_se
)

model_design_res$d_mean_bias_p <- calc_twotail_p(
  model_design_res$d_mean_bias_avg /
  model_design_res$d_mean_bias_se
)

model_design_res$d_male_effect_bias_p <- calc_twotail_p(
  model_design_res$d_male_effect_bias_avg /
  model_design_res$d_male_effect_bias_se
)

keep <- which(model_design_res$model %in% c(
  "hurdle_full", "lognormal_full", "hurdle_lognormal_full",
  "linked_hurdle_lognormal_full", "lm_lognormal",
  "lm_lognormal_plus1", "glm_hurdle")
)
model_bias_table <- model_design_res[keep,]

# rewrite l as t_0 estimates
# this is an important assumption, since our assertion is E(t0) = l.
tar <- model_bias_table$model %in% c("glm_hurdle", "hurdle_full")
model_bias_table$t0_mean_bias_avg[tar] <- model_bias_table$l_mean_bias_avg[tar]
model_bias_table$t0_mean_bias_se[tar] <- model_bias_table$l_mean_bias_se[tar]
model_bias_table$t0_mean_bias_p[tar] <- model_bias_table$l_mean_bias_p[tar]
model_bias_table$t0_male_effect_bias_avg[tar] <- model_bias_table$l_male_effect_bias_avg[tar]
model_bias_table$t0_male_effect_bias_se[tar] <- model_bias_table$l_male_effect_bias_se[tar]
model_bias_table$t0_male_effect_bias_p[tar] <- model_bias_table$l_male_effect_bias_p[tar]

# add study design designator
model_bias_table$design <- case_when(
  model_bias_table$long == FALSE ~ "cross-sectional",
  model_bias_table$long == TRUE ~ "longitudinal"
)

# order table based on design, then t0 bias
o <- order(model_bias_table$long, abs(model_bias_table$t0_mean_bias_avg))
model_bias_table <- model_bias_table[o,]

# add tags for whether to italicize the estimates or not

threshold <- 0.001

model_bias_table$t0_mean_bias_p_tag <- model_bias_table$t0_mean_bias_p < threshold
model_bias_table$t0_male_effect_bias_p_tag <- model_bias_table$t0_male_effect_bias_p < threshold
model_bias_table$k_mean_bias_p_tag <- model_bias_table$k_mean_bias_p < threshold
model_bias_table$k_male_effect_bias_p_tag <- model_bias_table$k_male_effect_bias_p < threshold
model_bias_table$d_mean_bias_p_tag <- model_bias_table$d_mean_bias_p < threshold
model_bias_table$d_male_effect_bias_p_tag <- model_bias_table$d_male_effect_bias_p < threshold

# format output numbers (convert from decades back to years)

model_bias_table$t0_mean_bias_avg <- sprintf("%1.1f", model_bias_table$t0_mean_bias_avg * 10)
model_bias_table$t0_mean_bias_se <- sprintf("%1.1f", model_bias_table$t0_mean_bias_se * 10)
model_bias_table$t0_mean_bias_p <- sprintf("%1.3f", model_bias_table$t0_mean_bias_p)

model_bias_table$t0_male_effect_bias_avg <- sprintf("%1.1f", model_bias_table$t0_male_effect_bias_avg * 10)
model_bias_table$t0_male_effect_bias_se <- sprintf("%1.1f", model_bias_table$t0_male_effect_bias_se * 10)
model_bias_table$t0_male_effect_bias_p <- sprintf("%1.3f", model_bias_table$t0_male_effect_bias_p)

model_bias_table$d_mean_bias_avg <- sprintf("%1.1f", model_bias_table$d_mean_bias_avg * 10)
model_bias_table$d_mean_bias_se <- sprintf("%1.1f", model_bias_table$d_mean_bias_se * 10)
model_bias_table$d_mean_bias_p <- sprintf("%1.3f", model_bias_table$d_mean_bias_p)

model_bias_table$d_male_effect_bias_avg <- sprintf("%1.1f", model_bias_table$d_male_effect_bias_avg * 10)
model_bias_table$d_male_effect_bias_se <- sprintf("%1.1f", model_bias_table$d_male_effect_bias_se * 10)
model_bias_table$d_male_effect_bias_p <- sprintf("%1.3f", model_bias_table$d_male_effect_bias_p)

model_bias_table$k_mean_bias_avg <- sprintf("%1.1f", model_bias_table$k_mean_bias_avg)
model_bias_table$k_mean_bias_se <- sprintf("%1.1f", model_bias_table$k_mean_bias_se)
model_bias_table$k_mean_bias_p <- sprintf("%1.3f", model_bias_table$k_mean_bias_p)

model_bias_table$k_male_effect_bias_avg <- sprintf("%1.1f", model_bias_table$k_male_effect_bias_avg)
model_bias_table$k_male_effect_bias_se <- sprintf("%1.1f", model_bias_table$k_male_effect_bias_se)
model_bias_table$k_male_effect_bias_p <- sprintf("%1.3f", model_bias_table$k_male_effect_bias_p)

model_bias_table$t0_mean_bias_p <- paste0("P=", model_bias_table$t0_mean_bias_p)
model_bias_table$t0_male_effect_bias_p <- paste0("P=", model_bias_table$t0_male_effect_bias_p)
model_bias_table$d_mean_bias_p <- paste0("P=", model_bias_table$d_mean_bias_p)
model_bias_table$d_male_effect_bias_p <- paste0("P=", model_bias_table$d_male_effect_bias_p)

tar <- which(model_bias_table$t0_mean_bias_p == "P=0.000")
model_bias_table$t0_mean_bias_p[tar] <- "P<0.001"
tar <- which(model_bias_table$t0_male_effect_bias_p == "P=0.000")
model_bias_table$t0_male_effect_bias_p[tar] <- "P<0.001"

tar <- which(model_bias_table$d_mean_bias_p == "P=0.000")
model_bias_table$d_mean_bias_p[tar] <- "P<0.001"
tar <- which(model_bias_table$d_male_effect_bias_p == "P=0.000")
model_bias_table$d_male_effect_bias_p[tar] <- "P<0.001"

tar <- which(model_bias_table$k_mean_bias_p == "P=0.000")
model_bias_table$k_mean_bias_p[tar] <- "P<0.001"
tar <- which(model_bias_table$k_male_effect_bias_p == "P=0.000")
model_bias_table$k_male_effect_bias_p[tar] <- "P<0.001"

# clean up the NaN values
model_bias_table$t0_mean_bias_se <- gsub("NaN|NA", "-", model_bias_table$t0_mean_bias_se)
model_bias_table$t0_mean_bias_p <- gsub("NaN|NA", "-", model_bias_table$t0_mean_bias_p)

model_bias_table$d_mean_bias_avg <- gsub("NaN|NA", "-", model_bias_table$d_mean_bias_avg)
model_bias_table$d_mean_bias_se <- gsub("NaN|NA", "-", model_bias_table$d_mean_bias_se)
model_bias_table$d_mean_bias_p <- gsub("NaN|NA", "-", model_bias_table$d_mean_bias_p)

model_bias_table$k_mean_bias_avg <- gsub("NaN|NA", "-", model_bias_table$k_mean_bias_avg)
model_bias_table$k_mean_bias_se <- gsub("NaN|NA", "-", model_bias_table$k_mean_bias_se)
model_bias_table$k_mean_bias_p <- gsub("NaN|NA", "-", model_bias_table$k_mean_bias_p)

model_bias_table$t0_male_effect_bias_se <- gsub("NaN|NA", "-", model_bias_table$t0_male_effect_bias_se)
model_bias_table$t0_male_effect_bias_p <- gsub("NaN|NA", "-", model_bias_table$t0_male_effect_bias_p)

model_bias_table$d_male_effect_bias_avg <- gsub("NaN|NA", "-", model_bias_table$d_male_effect_bias_avg)
model_bias_table$d_male_effect_bias_se <- gsub("NaN|NA", "-", model_bias_table$d_male_effect_bias_se)
model_bias_table$d_male_effect_bias_p <- gsub("NaN|NA", "-", model_bias_table$d_male_effect_bias_p)

model_bias_table$k_male_effect_bias_avg <- gsub("NaN|NA", "-", model_bias_table$k_male_effect_bias_avg)
model_bias_table$k_male_effect_bias_se <- gsub("NaN|NA", "-", model_bias_table$k_male_effect_bias_se)
model_bias_table$k_male_effect_bias_p <- gsub("NaN|NA", "-", model_bias_table$k_male_effect_bias_p)

# store key calculations for in-text reference

tar <- which(model_bias_table$model == "glm_hurdle" & model_bias_table$design == "cross-sectional")
calcs$biasLogitCrossOnsetEst <- model_bias_table[tar,]$t0_mean_bias_avg
calcs$biasLogitCrossOnsetSE <- model_bias_table[tar,]$t0_mean_bias_se
calcs$biasLogitCrossOnsetP <- model_bias_table[tar,]$t0_mean_bias_p

tar <- which(model_bias_table$model == "hurdle_full" & model_bias_table$design == "cross-sectional")
calcs$biasHurdleCrossOnsetEst <- model_bias_table[tar,]$t0_mean_bias_avg
calcs$biasHurdleCrossOnsetSE <- model_bias_table[tar,]$t0_mean_bias_se
calcs$biasHurdleCrossOnsetP <- model_bias_table[tar,]$t0_mean_bias_p

tar <- which(model_bias_table$model == "glm_hurdle" & model_bias_table$design == "longitudinal")
calcs$biasLogitLongOnsetEst <- model_bias_table[tar,]$t0_mean_bias_avg
calcs$biasLogitLongOnsetSE <- model_bias_table[tar,]$t0_mean_bias_se
calcs$biasLogitLongOnsetP <- model_bias_table[tar,]$t0_mean_bias_p

tar <- which(model_bias_table$model == "linked_hurdle_lognormal_full" & model_bias_table$design == "cross-sectional")
calcs$biasLHLNCrossOnsetEst <- model_bias_table[tar,]$t0_mean_bias_avg
calcs$biasLHLNCrossOnsetSE <- model_bias_table[tar,]$t0_mean_bias_se
calcs$biasLHLNCrossOnsetP <- model_bias_table[tar,]$t0_mean_bias_p

tar <- which(model_bias_table$model == "linked_hurdle_lognormal_full" & model_bias_table$design == "longitudinal")
calcs$biasLHLNLongOnsetEst <- model_bias_table[tar,]$t0_mean_bias_avg
calcs$biasLHLNLongOnsetSE <- model_bias_table[tar,]$t0_mean_bias_se
calcs$biasLHLNLongOnsetP <- model_bias_table[tar,]$t0_mean_bias_p

tar <- which(model_bias_table$model == "lm_lognormal" & model_bias_table$design == "longitudinal")
calcs$biasLNLongSexOnsetEst <- model_bias_table[tar,]$t0_male_effect_bias_avg
calcs$biasLNLongSexOnsetSE <- model_bias_table[tar,]$t0_male_effect_bias_se
calcs$biasLNLongSexOnsetP <- model_bias_table[tar,]$t0_male_effect_bias_p

tar <- which(model_bias_table$model == "lm_lognormal" & model_bias_table$design == "cross-sectional")
calcs$biasLNCrossOnsetEst <- model_bias_table[tar,]$t0_mean_bias_avg
calcs$biasLNCrossOnsetSE <- model_bias_table[tar,]$t0_mean_bias_se
calcs$biasLNCrossOnsetP <- model_bias_table[tar,]$t0_mean_bias_p
calcs$biasLNCrossSexOnsetEst <- model_bias_table[tar,]$t0_male_effect_bias_avg
calcs$biasLNCrossSexOnsetSE <- model_bias_table[tar,]$t0_male_effect_bias_se
calcs$biasLNCrossSexOnsetP <- model_bias_table[tar,]$t0_male_effect_bias_p
calcs$biasLNCrossProgEst <- model_bias_table[tar,]$d_mean_bias_avg
calcs$biasLNCrossProgSE <- model_bias_table[tar,]$d_mean_bias_se
calcs$biasLNCrossProgP <- model_bias_table[tar,]$d_mean_bias_p

tar <- which(model_bias_table$model == "lm_lognormal_plus1" & model_bias_table$design == "cross-sectional")
calcs$biasLNPCrossOnsetEst <- model_bias_table[tar,]$t0_mean_bias_avg
calcs$biasLNPCrossOnsetSE <- model_bias_table[tar,]$t0_mean_bias_se
calcs$biasLNPCrossOnsetP <- model_bias_table[tar,]$t0_mean_bias_p
calcs$biasLNPCrossSexOnsetEst <- model_bias_table[tar,]$t0_male_effect_bias_avg
calcs$biasLNPCrossSexOnsetSE <- model_bias_table[tar,]$t0_male_effect_bias_se
calcs$biasLNPCrossSexOnsetP <- model_bias_table[tar,]$t0_male_effect_bias_p
calcs$biasLNPCrossProgEst <- model_bias_table[tar,]$d_mean_bias_avg
calcs$biasLNPCrossProgSE <- model_bias_table[tar,]$d_mean_bias_se
calcs$biasLNPCrossProgP <- model_bias_table[tar,]$d_mean_bias_p

tar <- which(model_bias_table$model == "linked_hurdle_lognormal_full" & model_bias_table$design == "cross-sectional")
calcs$biasLHLNCrossSexOnsetEst <- model_bias_table[tar,]$t0_male_effect_bias_avg
calcs$biasLHLNCrossSexOnsetSE <- model_bias_table[tar,]$t0_male_effect_bias_se
calcs$biasLHLNCrossSexOnsetP <- model_bias_table[tar,]$t0_male_effect_bias_p

# rename models
model_bias_table$model <- case_when(
  model_bias_table$model == "glm_hurdle" ~ "logit(Pr(CAC\\,>\\,0))",
  model_bias_table$model == "hurdle_full" ~ "logit(Pr(CAC\\,>\\,0))$_i$ ",
  model_bias_table$model == "hurdle_lognormal_full" ~ "hurdle-lognormal$_i$",
  model_bias_table$model == "linked_hurdle_lognormal_full" ~ "linked hurdle-lognormal$_i$",
  model_bias_table$model == "lognormal_full" ~ "ln(CAC\\,|\\,CAC\\,>\\,0)$_i$",
  model_bias_table$model == "lm_lognormal" ~ "ln(CAC\\,|\\,CAC\\,>\\,0)",
  model_bias_table$model == "lm_lognormal_plus1" ~ "ln(CAC\\,+\\,1)",
)

# paste the SD's to the means with parentheses
model_bias_table$t0_mean_bias_avg <- paste0(model_bias_table$t0_mean_bias_avg, " (", model_bias_table$t0_mean_bias_se, ")")
model_bias_table$t0_male_effect_bias_avg <- paste0(model_bias_table$t0_male_effect_bias_avg, " (", model_bias_table$t0_male_effect_bias_se, ")")
model_bias_table$d_mean_bias_avg <- paste0(model_bias_table$d_mean_bias_avg, " (", model_bias_table$d_mean_bias_se, ")")
model_bias_table$d_male_effect_bias_avg <- paste0(model_bias_table$d_male_effect_bias_avg, " (", model_bias_table$d_male_effect_bias_se, ")")
model_bias_table$k_mean_bias_avg <- paste0(model_bias_table$k_mean_bias_avg, " (", model_bias_table$k_mean_bias_se, ")")
model_bias_table$k_male_effect_bias_avg <- paste0(model_bias_table$k_male_effect_bias_avg, " (", model_bias_table$k_male_effect_bias_se, ")")

# simplify missing coefficients
model_bias_table$d_mean_bias_avg <- gsub("- \\(-\\)", "-", model_bias_table$d_mean_bias_avg)
model_bias_table$d_male_effect_bias_avg <- gsub("- \\(-\\)", "-", model_bias_table$d_male_effect_bias_avg)
model_bias_table$k_mean_bias_avg <- gsub("- \\(-\\)", "-", model_bias_table$k_mean_bias_avg)
model_bias_table$k_male_effect_bias_avg <- gsub("- \\(-\\)", "-", model_bias_table$k_male_effect_bias_avg)

# use the significance tags to determine what values to italicize
tar <- which(model_bias_table$t0_mean_bias_p_tag)
model_bias_table$t0_mean_bias_avg[tar] <- paste0("\\textbf{", model_bias_table$t0_mean_bias_avg[tar], "}")

tar <- which(model_bias_table$t0_male_effect_bias_p_tag)
model_bias_table$t0_male_effect_bias_avg[tar] <- paste0("\\textbf{", model_bias_table$t0_male_effect_bias_avg[tar], "}")

tar <- which(model_bias_table$d_mean_bias_p_tag)
model_bias_table$d_mean_bias_avg[tar] <- paste0("\\textbf{", model_bias_table$d_mean_bias_avg[tar], "}")

tar <- which(model_bias_table$d_male_effect_bias_p_tag)
model_bias_table$d_male_effect_bias_avg[tar] <- paste0("\\textbf{", model_bias_table$d_male_effect_bias_avg[tar], "}")

# subset to exactly the variables to be shown
model_bias_table <- select(model_bias_table,
  model,
  design,
  t0_mean_bias_avg,
  t0_male_effect_bias_avg,
  d_mean_bias_avg,
  d_male_effect_bias_avg
)

# rename columns
model_bias_table <- rename(model_bias_table,
  `$t_0$ (years)` = t0_mean_bias_avg,
  `$t_0$ sex (years)` = t0_male_effect_bias_avg,
  `$d$ (years)` = d_mean_bias_avg,
  `$d$ sex (years)` = d_male_effect_bias_avg,
)

# output as tex format
tex_table <- model_bias_table
tex_table <- rbind(colnames(tex_table), tex_table)
writeLines(texttab(tex_table, hlines = c(1, 8)), "tableModelBiasMeans.tex")

print("created tableModelBiasMeans")



# export manuscript calculations

calcs$numPatients <- format(calcs$numPatients, big.mark = ",")
calcs$numObservations <- format(calcs$numObservations, big.mark = ",")
writeLines(prep_latex_variables(calcs), "calcsExperimentOutput.tex")
