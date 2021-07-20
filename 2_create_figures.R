
rm(list = ls())

source("project_support.R")

# initialize person

t0 <- 40 # age of cac beginning, in years

growth_period <- 60
delta_t <- 0.01 # in years
age <- rep(NA, growth_period/delta_t)
y <- rep(NA, growth_period/delta_t)

age[1] <- t0
for (i in 2:length(y)) {
  age[i] <- age[i-1] + delta_t
}

# growth process
cac <- rep(NA, length(age))
growing <- which(age >= t0)
if (length(growing) < length(age)) cac[-growing] <- 0

# doubling time schedule (manual):
simplex <- function(x) x/sum(x)
d_vec <- rep(c(2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10), times = floor(simplex(c(0.2, 0.2, 3, 0.2, 0.2, 0.1, 2, 0.2, 0.2, 5)) * length(y)))
d_vec <- c(d_vec, rep(d_vec[length(d_vec)], length(age) - length(d_vec)))

k_avg <- rep(NA, length(y))
k_vec <- log(2)/d_vec # growth rate, in log-y per year
y[1] <- 1
for (i in 2:length(y)) {
  k_avg[i] <- mean(k_vec[1:i])
  y[i] <- y[i-1] + (k_vec[i-1] * y[i-1]) * delta_t
}

# png("figureSingleTrajectory.png", res = 300, units = "in", height = 5, width = 5)

tikz("figureSingleTrajectory.tex", height = 4, width = 5.5)

age_off <- 3

par(mfrow = c(1, 1))
plot(age, y, log = "y", type = "n", las = 1, xlim = c(t0 - age_off, 80), yaxt = "n",
  frame.plot = FALSE, ylim = c(1, 1000), xlab = "patient age",
  ylab = "CAC score (Agatston units, log scale)", xaxt = "n")
axis(2, at = c(1, 10, 100, 1000), labels = c("1", "", "", ""), las = 1)
axis(2, at = seq(1, 9, by = 1), labels = NA, tcl = (-0.2))
axis(2, at = seq(10, 100, by = 10), labels = NA, tcl = (-0.2))
axis(2, at = seq(100, 1000, by = 100), labels = NA, tcl = (-0.2))
i <- 3000

d_avg <- log(2)/k_avg[i]

axis(1, tcl = 0, at = c(t0 - age_off, 80), labels = FALSE)
axis(1, at = c(t0, age[i], age[i] + d_avg),
  labels = c("$t_0$", "$t$", "$t + d$"), xpd = TRUE)
axis(2, at = c(y[i], 2*y[i]), labels = c("$y$", "$2y$"), las = 1)

points(age, y, type = "l", lwd = 2.5, col = "dodgerblue")

polygon(
  c(age[i], age[i], age[i] + d_avg),
  c(y[i], 2 * y[i], 2 * y[i]),
  border = "gray", col = NULL
)
boxright <- d_avg * 0.15
boxdown <- y[i] * 0.3
polygon(
  c(age[i], age[i] + boxright, age[i] + boxright, age[i]),
  c(2 * y[i] - boxdown, 2 * y[i] - boxdown, 2 * y[i], 2 * y[i]),
  border = "gray", col = NULL
)
text(age[i] + d_avg/2, 2 * y[i], labels = "$d$", pos = 3)
text(age[i], 3/2 * y[i], labels = "$y$", pos = 2)

curve(exp(k_avg[i] * (x - t0)), lty = 1, add = TRUE, from = t0, to = age[i] + d_avg)
lines(c(age[i], age[i]), c(1, y[i]), lty = 2)
lines(c(t0 - age_off, age[i]), c(y[i], y[i]), lty = 2)
text(56, 6, labels = "$k = \\displaystyle \\frac{\\mathrm{ln}\\,y}{t-t_0}$")

points(age[i], y[i], pch = 16)

dev.off()

print("figureSingleTrajectory created")




# plot fig:varianceBias

calcs <- list()

n_samples <- 100
N_ind <- 5000

calcs$biasNumPatients <- format(N_ind, big.mark = ",")
calcs$biasNumSamples <- as.character(n_samples)

t0_pop <- 0
k <- 1
t0_ind <- rnorm(N_ind, 0, 1)

t0_ind_sigma <- seq(0, 3.2, by = 0.1)

out <- data.frame(
  t0_ind_sigma = t0_ind_sigma,
  t0_bias_mean = rep(NA, length(t0_ind_sigma)),
  t0_bias_lb = rep(NA, length(t0_ind_sigma)),
  t0_bias_ub = rep(NA, length(t0_ind_sigma)),
  k_bias_mean = rep(NA, length(t0_ind_sigma)),
  k_bias_lb = rep(NA, length(t0_ind_sigma)),
  k_bias_ub = rep(NA, length(t0_ind_sigma)),
  d_bias_mean = rep(NA, length(t0_ind_sigma)),
  d_bias_lb = rep(NA, length(t0_ind_sigma)),
  d_bias_ub = rep(NA, length(t0_ind_sigma)),
  t0_plus1_bias_mean = rep(NA, length(t0_ind_sigma)),
  t0_plus1_bias_lb = rep(NA, length(t0_ind_sigma)),
  t0_plus1_bias_ub = rep(NA, length(t0_ind_sigma)),
  k_plus1_bias_mean = rep(NA, length(t0_ind_sigma)),
  k_plus1_bias_lb = rep(NA, length(t0_ind_sigma)),
  k_plus1_bias_ub = rep(NA, length(t0_ind_sigma)),
  d_plus1_bias_mean = rep(NA, length(t0_ind_sigma)),
  d_plus1_bias_lb = rep(NA, length(t0_ind_sigma)),
  d_plus1_bias_ub = rep(NA, length(t0_ind_sigma))
)

for (i in seq_along(t0_ind_sigma)) {

  t0 <- t0_pop + t0_ind * t0_ind_sigma[i]

  t0_bias_est <- rep(NA, n_samples)
  k_bias_est <- rep(NA, n_samples)
  d_bias_est <- rep(NA, n_samples)
  t0_plus1_bias_est <- rep(NA, n_samples)
  k_plus1_bias_est <- rep(NA, n_samples)
  d_plus1_bias_est <- rep(NA, n_samples)

  for (j in 1:n_samples) {

    ages <- rnorm(N_ind, 2, 1)
    y <- exp(k * (ages - t0))
    y[ages < t0] <- 0

    d <- data.frame(age = ages, y = y)
    d$log_y <- log(d$y)
    d$log_yplus1 <- log(d$y + 1)

    m_lm <- lm(log_y ~ age, data = d[d$y > 0,])
    t0_bias_est[j] <- (-coef(m_lm)[1]/coef(m_lm)[2]) - t0_pop
    k_bias_est[j] <- coef(m_lm)[2] - k
    d_bias_est[j] <- log(2)/coef(m_lm)[2] - log(2)/k

    m_lm <- lm(log_yplus1 ~ age, data = d[d$y > 0,])
    t0_plus1_bias_est[j] <- (-coef(m_lm)[1]/coef(m_lm)[2]) - t0_pop
    k_plus1_bias_est[j] <- coef(m_lm)[2] - k
    d_plus1_bias_est[j] <- log(2)/coef(m_lm)[2] - log(2)/k

  }

  out$t0_bias_mean[i] <- mean(t0_bias_est)
  out$t0_bias_lb[i] <- HPDI(t0_bias_est)[1]
  out$t0_bias_ub[i] <- HPDI(t0_bias_est)[2]
  out$k_bias_mean[i] <- mean(k_bias_est)
  out$k_bias_lb[i] <- HPDI(k_bias_est)[1]
  out$k_bias_ub[i] <- HPDI(k_bias_est)[2]
  out$d_bias_mean[i] <- mean(d_bias_est)
  out$d_bias_lb[i] <- HPDI(d_bias_est)[1]
  out$d_bias_ub[i] <- HPDI(d_bias_est)[2]
  
  out$t0_plus1_bias_mean[i] <- mean(t0_plus1_bias_est)
  out$t0_plus1_bias_lb[i] <- HPDI(t0_plus1_bias_est)[1]
  out$t0_plus1_bias_ub[i] <- HPDI(t0_plus1_bias_est)[2]
  out$k_plus1_bias_mean[i] <- mean(k_plus1_bias_est)
  out$k_plus1_bias_lb[i] <- HPDI(k_plus1_bias_est)[1]
  out$k_plus1_bias_ub[i] <- HPDI(k_plus1_bias_est)[2]
  out$d_plus1_bias_mean[i] <- mean(d_plus1_bias_est)
  out$d_plus1_bias_lb[i] <- HPDI(d_plus1_bias_est)[1]
  out$d_plus1_bias_ub[i] <- HPDI(d_plus1_bias_est)[2]
  
  if (i %% 5 == 0) print(i)

}

biasDiffOnset <- mean(out$t0_plus1_bias_mean - out$t0_bias_mean) # 0.589, six years
biasDiffDoubling <- mean(out$d_plus1_bias_mean - out$d_bias_mean) # 0.096, one year

# export manuscript calculations

calcs$biasDiffOnset <- paste(round(10 * biasDiffOnset, 1), "years")
calcs$biasDiffDoubling <- paste(round(10 * biasDiffDoubling, 1), "years")
calcs$biasNumPatients <- format(N_ind, big.mark = ",")
calcs$biasNumSamples <- as.character(n_samples)
writeLines(prep_latex_variables(calcs), "calcsVarianceBias.tex")


# png("figureVarianceBias.png", res = 300, height = 5, width = 10, units = "in")

tikz("figureVarianceBias.tex", height = 3.5, width = 7)

par(mfrow = c(1, 2))
plot(out$t0_ind_sigma, out$t0_bias_mean, type = "l", xlim = c(0, 3), ylim = c(-5, 0),
  xlab = "$\\sigma_{t0}$ onset variability (years)", ylab = "$t_0$ onset bias (years)", yaxt = "n")
axis(2, at = seq(-5, 0, by = 1), labels = seq(-50, 0, by = 10))
polygon(c(out$t0_ind_sigma, rev(out$t0_ind_sigma)),
  c(out$t0_bias_lb, rev(out$t0_bias_ub)),
  border = NA, col = col_alpha("gray", 0.4))

points(out$t0_ind_sigma, out$t0_plus1_bias_mean, type = "l", col = "red")
polygon(c(out$t0_ind_sigma, rev(out$t0_ind_sigma)),
  c(out$t0_plus1_bias_lb, rev(out$t0_plus1_bias_ub)),
  border = NA, col = col_alpha("red", 0.4))

abline(h = 0, lty = 2)

text(2.5, -1.8, "ln(CAC|CAC>0)")
text(1.25, -2.5, "ln(CAC+1)", col = "red")

plot(out$t0_ind_sigma, out$d_bias_mean, type = "l",
  xlab = "$\\sigma_{t0}$ onset variability (years)", ylab = "$d$ doubling-time bias (years)",
  xlim = c(0, 3), ylim = c(0, 1), yaxt = "n")
axis(2, at = seq(0, 1, by = 0.2), labels = seq(0, 10, by = 2))
polygon(c(out$t0_ind_sigma, rev(out$t0_ind_sigma)),
  c(out$d_bias_lb, rev(out$d_bias_ub)),
  border = NA, col = col_alpha("gray", 0.4))

points(out$t0_ind_sigma, out$d_plus1_bias_mean, type = "l", col = "red")
polygon(c(out$t0_ind_sigma, rev(out$t0_ind_sigma)),
  c(out$d_plus1_bias_lb, rev(out$d_plus1_bias_ub)),
  border = NA, col = col_alpha("red", 0.4))

abline(h = 0, lty = 2)

text(2.5, 0.3, "ln(CAC|CAC>0)")
text(1.25, 0.45, "ln(CAC+1)", col = "red")

dev.off()

print("figureVarianceBias created")





# only fit the full lhln model
models <- models["linked_hurdle_lognormal_full"]

pars <- list(
  seed = 5,
  sim = list(
    longitudinal =  TRUE,
    N_ind        =     5,
    l_pop_mu     =     0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =   0.5,
    l_sex_sigma  =   1.0,
    b_mu         =   1.5, # b = 1/v, v is a scale term in decades
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7, # on the log scale
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  ),
  stan = list(
    # prior model parameters
    l_pop_mu          = 0,
    l_pop_sigma       = 0.18,
    l_ind_sigma_rate  = 1/0.5,
    l_sex_sigma_rate  = 1/1.0,
    b_mu              = 1.0,
    b_sigma           = 0.5,
    d_pop_mu          = -0.7,
    d_pop_sigma       = 0.01,
    d_ind_sigma_rate  = 1/0.04,
    d_sex_sigma_rate  = 1/0.01,
    s_rate            = 1/0.1,
    # sampling parameters
    n_chains = 3,
    n_cores = 3,
    n_iter = 2000,
    adapt_delta = 0.8,
    max_treedepth = 15
  )
)

res_long <- run_experiment(pars, fit_models = TRUE)

pars <- list(
  seed = 5,
  sim = list(
    longitudinal = FALSE,
    N_ind        =    40,
    l_pop_mu     =     0,
    l_pop_sigma  =  0.18,
    l_ind_sigma  =   0.5,
    l_sex_sigma  =  0.01,
    b_mu         =   1.5, # b = 1/v, v is a scale term in decades
    b_sigma      =   0.5,
    d_pop_mu     =  -0.7, # on the log scale
    d_pop_sigma  =  0.01,
    d_ind_sigma  =  0.04,
    d_sex_sigma  =  0.01,
    s            =   0.1
  ),
  stan = list(
    # prior model parameters
    l_pop_mu          = 0,
    l_pop_sigma       = 0.18,
    l_ind_sigma_rate  = 1/0.5,
    l_sex_sigma_rate  = 1/0.01,
    b_mu              = 1.5,
    b_sigma           = 0.5,
    d_pop_mu          = -0.7,
    d_pop_sigma       = 0.01,
    d_ind_sigma_rate  = 1/0.04,
    d_sex_sigma_rate  = 1/0.01,
    s_rate            = 1/0.1,
    # sampling parameters
    n_chains = 4,
    n_cores = 4,
    n_iter = 14000,
    adapt_delta = 0.90,
    max_treedepth = 10
  )
)

res_wide <- run_experiment(pars, fit_models = TRUE)

# png("figureLongCrossExamples.png", res = 300, units = "in", height = 5, width = 9.5)

tikz("figureLongCrossExamples.tex", height = 3.5, width = 6.8)

par(mfrow = c(1, 2))

res <- res_long
d <- res$data$obs
d$sampled <- FALSE
d$ln_cac <- log(d$cac)
d$ln_cac_plus1 <- log(d$cac + 1)
m_ln_plus1 <- lm(ln_cac_plus1 ~ age_su, data = d)
m_ln <- lm(ln_cac ~ age_su, data = d[d$cac > 0,])
plot(d$age_su, d$cac, log = "y", las = 1, xaxt = "n", xlim = c(-2, 5), ylim = c(1, 1000), ylab = "CAC score (Agatston units, log scale)", xlab = "patient age (years)", col = col_alpha("black", 0.4))
axis(1, at = seq(-5, 5, by = 1), labels = seq(0, 100, by = 10))
axis(2, at = seq(1, 10, by = 1), labels = NA, tcl = (-0.3))
axis(2, at = seq(10, 100, by = 10), labels = NA, tcl = (-0.3))
axis(2, at = seq(100, 1000, by = 100), labels = NA, tcl = (-0.3))
curve(exp(coef(m_ln)[1] + coef(m_ln)[2] * x), add = TRUE, from = (-coef(m_ln)[1]/coef(m_ln)[2]), to = 100)
curve(exp(coef(m_ln_plus1)[1] + coef(m_ln_plus1)[2] * x), add = TRUE, col = "red", from = (-coef(m_ln_plus1)[1]/coef(m_ln_plus1)[2]), to = 100)
points(mean(res$sim$t0), 1, pch = 17)
for (i in 1:res$sim$N_ind) {
  curve(exp(res$estimates$linked_hurdle_lognormal_full$k_mean[i] * (x - res$estimates$linked_hurdle_lognormal_full$t0_mean[i])), add = TRUE, col = col_alpha("dodgerblue", 0.7), from = res$estimates$linked_hurdle_lognormal_full$t0_mean[i], to = 100)
}

res <- res_wide
d <- res$data$obs
d$sampled <- FALSE
d$ln_cac <- log(d$cac)
d$ln_cac_plus1 <- log(d$cac + 1)
m_ln_plus1 <- lm(ln_cac_plus1 ~ age_su, data = d)
m_ln <- lm(ln_cac ~ age_su, data = d[d$cac > 0,])
plot(d$age_su, d$cac, log = "y", las = 1, xaxt = "n", xlim = c(-2, 5), ylim = c(1, 1000),
  ylab = "CAC score (Agatston units, log scale)", xlab = "patient age (years)", col = col_alpha("black", 0.4))
axis(1, at = seq(-5, 5, by = 1), labels = seq(0, 100, by = 10))
axis(2, at = seq(1, 10, by = 1), labels = NA, tcl = (-0.3))
axis(2, at = seq(10, 100, by = 10), labels = NA, tcl = (-0.3))
axis(2, at = seq(100, 1000, by = 100), labels = NA, tcl = (-0.3))
curve(exp(coef(m_ln)[1] + coef(m_ln)[2] * x), add = TRUE)
curve(exp(coef(m_ln_plus1)[1] + coef(m_ln_plus1)[2] * x), add = TRUE, col = "red")
points(mean(res$sim$t0), 1, pch = 17)
for (i in 1:res$sim$N_ind) {
  curve(exp(res$estimates$linked_hurdle_lognormal_full$k_mean[i] * (x - res$estimates$linked_hurdle_lognormal_full$t0_mean[i])), add = TRUE, col = col_alpha("dodgerblue", 0.3), from = res$estimates$linked_hurdle_lognormal_full$t0_mean[i], to = 100)
}

dev.off()

print("figureLongCrossExamples created")
