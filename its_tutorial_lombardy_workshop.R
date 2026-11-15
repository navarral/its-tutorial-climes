# ==================================================================
# TUTORIAL: INTERRUPTED TIME SERIES (ITS) WITH brms
# Question: did a policy lower the preterm rate in Lombardy?
# Data: simulated. Period: 4 weeks. Unit: province.
# Run time: about 10 minutes in class (plus first-time model compiling).
# ==================================================================

# ------------------------------------------------------------------
# PART 0. WORDS USED IN THIS TUTORIAL
# ------------------------------------------------------------------
# Period:            one 4-week block. One year has 13 periods.
# Births:            all live births in one province in one period.
# Preterm births:    births before 37 weeks of pregnancy.
# Preterm rate:      preterm births divided by births.
# Policy start:      the first period after the policy begins.
# Pre-policy:        all periods before the policy start.
# Post-policy:       all periods from the policy start onward.
# Level change:      the immediate jump in the preterm rate at the policy start.
# Slope change:      the extra change per year in the rate after the policy start.
# Counterfactual:    the preterm rate we expect without the policy.
# Prior:             what we believe about a value before we see the data.
# Posterior:         what we believe about a value after we see the data.
# Credible interval: the range that holds the true value with 95% probability.
# ROPE:              the range of changes too small to matter in practice.
# pd:                the posterior probability that an effect has a given sign.

# ------------------------------------------------------------------
# PART 1. THE IDEA OF INTERRUPTED TIME SERIES
# ------------------------------------------------------------------
# We measure the preterm rate in many periods before the policy start.
# We use the pre-policy periods to project the counterfactual.
# We compare the counterfactual with the observed post-policy rate.
# The gap between the two lines is the policy effect.
# The model has three parts:
#   1. A pre-policy slope (time_c).
#   2. A level change at the policy start (policy).
#   3. A slope change after the policy start (time_after).
# A policy can change the level, the slope, or both.
# We add two seasonal terms (sin1, cos1) because preterm rates vary by season.
# We add province intercepts because provinces differ in baseline rate.

# ------------------------------------------------------------------
# PART 2. THE BAYESIAN IDEA, IN BRIEF
# ------------------------------------------------------------------
# A Bayesian model starts with a prior for each value.
# The data update each prior into a posterior.
# We read the posterior directly as probabilities.
# Example: "The probability that the policy lowered the rate is 98%."
# MCMC is the standard way to compute a posterior. MCMC is slow.
# ADVI is a fast approximation of the posterior. ADVI takes seconds.
# We use ADVI here ("fullrank") so the class can run all the code.
# ADVI can understate uncertainty. Use MCMC for a real analysis.
# To switch, replace algorithm = "fullrank" with algorithm = "sampling".

# ------------------------------------------------------------------
# PART 3. SET UP
# ------------------------------------------------------------------
# install.packages(c("brms", "tidyverse", "bayestestR"))   # Run once.
# brms needs the rstan package.
library(brms)
library(tidyverse)
library(bayestestR)
set.seed(2026)

# Setting I like for ggplot: #648fff
theme_set(
  theme_bw(base_size = 16) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(colour = "black")
    )
)

# Approximate mean births per province per period (invented for teaching).
provinces <- data.frame(
  province    = c("Milano", "Brescia", "Bergamo", "Monza e Brianza", "Varese",
                  "Como", "Pavia", "Mantova", "Cremona", "Lecco", "Lodi",
                  "Sondrio"),
  births_annual = c(
    26984, 10558, 9425, 7052,
    6971, 4884, 3936, 3165,
    2625, 2586, 1949, 1453
  )
)

# Add mean level every 4 weeks
provinces$births_mean <- as.integer(provinces$births_annual / 13)

# ------------------------------------------------------------------
# PART 4. SIMULATE DATA
# ------------------------------------------------------------------
# One function simulates one study. We reuse it for the power calculation.
# Arguments:
#   n_pre         number of pre-policy periods
#   n_post        number of post-policy periods
#   rr            rate ratio of the level change (0.80 means a 20% decrease)
#   slope_rr      rate ratio of the slope change per year
#                 (0.90 means an extra 10% decrease each year)
#   baseline      preterm rate before the policy (0.07 means 7%)
#   pre_slope     yearly slope of the pre-policy log-odds
#   births_scale  multiplier for births (0.5 means half the births)
#   rho           correlation between neighbouring periods (autocorrelation)
#   noise_sd      size of random variation beyond binomial chance
#
# Example: rr = 0.80 and slope_rr = 0.90 give a 20% drop at the policy start.
# After one year, the rate sits 28% below the counterfactual (0.80 x 0.90).
# After two years, the rate sits 35% below it (0.80 x 0.90 x 0.90).

simulate_its <- function(n_pre = 26, n_post = 26, rr = 0.80, slope_rr = 0.90,
                         baseline = 0.07, pre_slope = 0,
                         births_scale = 1, rho = 0.3, noise_sd = 0.05) {
  n <- n_pre + n_post
  df <- expand.grid(period = 1:n, province = provinces$province)
  idx <- as.integer(df$province)

  # Time is in years. Time 0 is the policy start.
  df$time_c     <- (df$period - (n_pre + 1)) / 13
  df$policy     <- as.numeric(df$period > n_pre)
  df$time_after <- pmax(df$time_c, 0)

  # Seasonal terms repeat every 13 periods.
  df$sin1 <- sin(2 * pi * df$period / 13)
  df$cos1 <- cos(2 * pi * df$period / 13)

  # The model works on the log-odds scale. We convert both rate ratios.
  # The slope conversion is approximate because preterm births are rare.
  level_change <- qlogis(baseline * rr) - qlogis(baseline)
  slope_change <- log(slope_rr)

  # Each province has its own baseline. Each province has its own noise.
  province_effect <- rnorm(nrow(provinces), mean = 0, sd = 0.10)
  noise <- unlist(lapply(seq_len(nrow(provinces)), function(i) {
    as.numeric(stats::filter(rnorm(n, 0, noise_sd), rho, method = "recursive"))
  }))

  log_odds <- qlogis(baseline) + province_effect[idx] +
    pre_slope * df$time_c + slope_change * df$time_after +
    level_change * df$policy +
    0.05 * df$sin1 + 0.03 * df$cos1 + noise

  df$births  <- pmax(rpois(nrow(df), provinces$births_mean[idx] * births_scale), 10)
  df$preterm <- rbinom(nrow(df), size = df$births, prob = plogis(log_odds))
  df
}

n_pre  <- 26   # 26 periods = 2 years before the policy
n_post <- 26   # 26 periods = 2 years after the policy
df <- simulate_its(n_pre, n_post, rr = 0.80, slope_rr = 0.90)
head(df)

# ------------------------------------------------------------------
# PART 5. LOOK AT THE DATA
# ------------------------------------------------------------------
# We always plot the data before we fit any model.
df %>%
  mutate(
    rate = 100 * preterm / births,
    province = factor(province, levels = provinces$province)
  ) %>%
  ggplot(aes(x = period, y = rate)) +
  geom_line(colour = "black", linewidth = 0.5, alpha = 0.5) +
  geom_point(colour = "black", size = 1, alpha = 0.8) +
  geom_vline(
    xintercept = n_pre + 0.5,
    linetype = "dashed",
    colour = "#fe6100"
  ) +
  facet_wrap(~ province) +
  labs(
    title = "Simulated preterm rates by province",
    subtitle = "The dashed line marks policy introduction",
    x = "Time [4 weeks]",
    y = "Preterm rate [%]"
  )

# We also add up all provinces to get one Lombardy series.
lombardy <- aggregate(cbind(preterm, births) ~ period, data = df, FUN = sum)
lombardy$rate <- lombardy$preterm / lombardy$births

plot(lombardy$period, 100 * lombardy$rate, type = "b", pch = 16,
     xlab = "Time [4 weeks]", ylab = "Preterm rate [%]",
     main = "Simulated preterm rate in Lombardy")
abline(v = n_pre + 0.5, lty = 2)   # The dashed line marks the policy start.

# ------------------------------------------------------------------
# PART 6. DESIGN CHOICES
# ------------------------------------------------------------------
# Choice 1. Outcome scale.
#   We model preterm births out of births with a binomial family.
#   The binomial family respects the counts in small provinces.
# Choice 2. Time scale.
#   We measure time in years and centre it at the policy start.
#   The policy coefficient then equals the level change at the policy start.
# Choice 3. Level change and slope change.
#   The level change measures an immediate effect.
#   The slope change measures an effect that grows or fades over time.
#   The slope change needs many post-policy periods to detect (see Part 9).
# Choice 4. Seasonality.
#   We use sine and cosine terms with a 13-period cycle.
#   Seasonal terms need at least one full year before the policy start.
# Choice 5. Provinces.
#   A random intercept lets each province keep its own baseline rate.
#   Small provinces borrow strength from large provinces.
# Choice 6. Extra variation.
#   Real data vary more than chance alone predicts.
#   The binomial family ignores this extra variation, so intervals run narrow.
#   The beta-binomial family absorbs it but needs starting values in ADVI.
# Choice 7. Autocorrelation.
#   Neighbouring periods often resemble each other.
#   We check autocorrelation after fitting (Part 8).
# Choice 8. Phase-in period.
#   Policies often take months to work.
#   Consider dropping periods right after the policy start.
# Choice 9. Other events.
#   Other events can mimic a policy effect. COVID-19 in 2020 is one example.
#   A control series (a region without the policy) protects against this.
# Choice 10. Practical relevance.
#   A tiny effect can look certain in a large dataset.
#   We define a ROPE: changes under 5% count as no practical change.

# ------------------------------------------------------------------
# PART 7. FIT THE BAYESIAN MODEL
# ------------------------------------------------------------------
# "0 + Intercept" makes the intercept the rate at the policy start.
# The first run compiles the model. Compiling takes about one minute.

formula_its <- bf(
  preterm | trials(births) ~ 0 + Intercept + time_c + policy + time_after +
    sin1 + cos1 + (1 | province)
)

# Weakly informative priors. Each prior states a plausible range.
priors_its <- c(
  prior(normal(-2.6, 1),   class = "b", coef = "Intercept"),  # near 7% preterm
  prior(normal(0, 0.5),    class = "b"),                      # effects on log-odds
  prior(exponential(10),   class = "sd")                      # province variation
)

fit <- brm(
  formula_its, data = df, family = binomial(), prior = priors_its,
  # Uncomment this line and comment 4 lines below to run MCMC
  algorithm = "sampling", iter = 2000, warmup = 1000, chains = 2, cores = 2, threads = 2,
  # algorithm = "fullrank",        # fast ADVI instead of MCMC 
  # iter = 20000,                # maximum number of ADVI steps
  # output_samples = 2000,       # posterior draws kept
  # tol_rel_obj = 0.001,         # stopping rule
  seed = 2026, refresh = 0
)

# Stan may print ADVI messages and warnings
summary(fit)

# Check fit
pp_check(fit, ndraws = 100)

# Check overdispersion: if bionomial is well-adjusted we should see the 95% of 
# real points within the 95% bands
pp_check(fit, type = "intervals")

# ------------------------------------------------------------------
# PART 8. READ THE RESULTS
# ------------------------------------------------------------------
dr <- as.data.frame(as_draws_df(fit))

# Rate ratio of the level change at the policy start. Below 1 means a decrease.
rate_ratio <- plogis(dr$b_Intercept + dr$b_policy) / plogis(dr$b_Intercept)

# ---- Summarise ITS model fits ---------------------------------------
# source 1: https://easystats.github.io/bayestestR/articles/
# source 2: https://easystats.github.io/bayestestR/articles/guidelines.html
#
# We report four numbers for each effect:
#   Median:         the central estimate of the rate ratio.
#   95% HDI:        the credible interval (highest density interval).
#   pd:             the probability that the effect is negative.
#   ROPE_Percentage: the share of the 95% HDI inside the ROPE.
#
# Decision rules for the ROPE:
#   The HDI lies fully outside the ROPE: the effect matters in practice.
#   The HDI lies fully inside the ROPE: the effect does not matter in practice.
#   The HDI overlaps the ROPE edge: the data cannot decide.
#
# We set the ROPE to rate ratios between 0.95 and 1.05.
# We work on the log scale, because the model works on that scale.
rope_ratio <- c(0.95, 1.05)

draws_log <- data.frame(
  level_change = log(rate_ratio),   # log rate ratio at the policy start
  slope_change = dr$b_time_after    # log rate ratio per year
)

summary_its <- describe_posterior(
  draws_log,
  centrality = "median",
  ci = 0.95, ci_method = "hdi",
  test = c("p_direction", "rope"),
  rope_range = log(rope_ratio), rope_ci = 0.95
)

# We convert the median and the HDI back to rate ratios.
summary_its %>%
  mutate(across(c(Median, CI_low, CI_high), exp)) %>%
  select(Parameter, Median, CI_low, CI_high, pd, ROPE_Percentage) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print()

# Question: How well does the model pick up the true effect? 

cat("\nTrue level change in the simulation:  rate ratio 0.80\n")
cat("True slope change in the simulation:  rate ratio 0.90 per year\n")

# Total effect at the end of the post-policy window.
# We use the fitted and counterfactual rates from the next block.
n <- n_pre + n_post
X <- cbind(Intercept  = 1,
           time_c     = (1:n - (n_pre + 1)) / 13,
           policy     = as.numeric(1:n > n_pre),
           time_after = pmax((1:n - (n_pre + 1)) / 13, 0),
           sin1       = sin(2 * pi * (1:n) / 13),
           cos1       = cos(2 * pi * (1:n) / 13))
B <- as.matrix(dr[, paste0("b_", colnames(X))])

# The counterfactual switches the policy off.
X0 <- X                          
X0[, c("policy", "time_after")] <- 0

rate_fit <- plogis(B %*% t(X))
rate_cf  <- plogis(B %*% t(X0))

rr_end <- rate_fit[, n] / rate_cf[, n]
cat("\nRate ratio at the end of the post-policy window (median, 95% interval):\n")
print(round(quantile(rr_end, c(0.5, 0.025, 0.975)), 3))

# Plot the fitted line and the counterfactual.
# colors plot: #fe6100, #648fff
tibble(
  period    = 1:n,
  observed  = 100 * lombardy$rate,
  fit_median    = 100 * apply(rate_fit, 2, median),
  fit_low    = 100 * apply(rate_fit, 2, quantile, 0.025),
  fit_high   = 100 * apply(rate_fit, 2, quantile, 0.975),
  cf_median = 100 * apply(rate_cf, 2, median),
  cf_low    = 100 * apply(rate_cf, 2, quantile, 0.025),
  cf_high   = 100 * apply(rate_cf, 2, quantile, 0.975)
) %>%
  ggplot(aes(x = period)) +
  geom_ribbon(aes(ymin = cf_low, ymax = cf_high, alpha = period > n_pre), fill = "#648fff") +
  geom_line(aes(y = cf_median), colour = "#648fff", linetype = "dotted", linewidth = 1) +
  geom_line(aes(y = fit_median), linewidth = 1) +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_high), fill = "black", alpha = 0.2) +
  geom_point(aes(y = observed), colour = "black", size = 2) +
  geom_vline(xintercept = n_pre + 0.5, linetype = "dashed") +
  labs(
    title = "Observed rate, fitted rate, and counterfactual",
    x = "Time [4 weeks]",
    y = "Preterm rate [%]"
  ) +
  scale_alpha_manual(
    values = c("FALSE" = 0, "TRUE" = 0.2),
    guide = "none"
  )

# Check autocorrelation. Bars outside the dashed lines signal a problem.
# A problem means the credible interval is too narrow.
residual_rate <- lombardy$rate - apply(rate_fit, 2, median)
acf(residual_rate, main = "Autocorrelation of residuals")

# Provincial plot
# Keep the same row order for predictions and observed rates.
plot_data <- df %>%
  arrange(province, period)

stopifnot(all(plot_data$births > 0))

# Switch off the intervention for the counterfactual.
counterfactual_data <- plot_data %>%
  mutate(
    policy = 0,
    time_after = 0
  )

# Predict expected counts, including province effects.
count_fit <- posterior_epred(
  fit,
  newdata = plot_data,
  re_formula = NULL
)

count_cf <- posterior_epred(
  fit,
  newdata = counterfactual_data,
  re_formula = NULL
)

# Convert expected counts into percentages.
rate_fit <- 100 * sweep(count_fit, 2, plot_data$births, "/")
rate_cf  <- 100 * sweep(count_cf, 2, plot_data$births, "/")

# Summarise posterior draws for each province and period.
plot_data <- plot_data %>%
  mutate(
    observed = 100 * preterm / births,
    
    fit_median = apply(rate_fit, 2, median),
    fit_low = apply(rate_fit, 2, quantile, probs = 0.025),
    fit_high = apply(rate_fit, 2, quantile, probs = 0.975),
    
    cf_median = apply(rate_cf, 2, median),
    cf_low = apply(rate_cf, 2, quantile, probs = 0.025),
    cf_high = apply(rate_cf, 2, quantile, probs = 0.975)
  )

ggplot(plot_data, aes(x = period)) +
  geom_ribbon( aes(ymin = fit_low, ymax = fit_high), fill = "black", alpha = 0.5  ) +
  geom_ribbon( data = filter(plot_data, period > n_pre), aes(ymin = cf_low, ymax = cf_high), fill = "#648fff", alpha = 0.5) +
  geom_line( aes(y = cf_median), colour = "#648fff", linetype = "dotted", linewidth = 0.8) +
  geom_line( aes(y = fit_median), colour = "black", linewidth = 0.8) +
  geom_point( aes(y = observed), colour = "black", size = 1) +
  geom_vline( xintercept = n_pre + 0.5, linetype = "dashed" ) +
  facet_wrap(~ province, ncol = 3) +
  labs(
    title = "Preterm rates by province",
    subtitle = "Black: fitted. Blue: counterfactual. Shading: 95% credible intervals.",
    x = "Time [4 weeks]",
    y = "Preterm rate [%]"
  )

# ------------------------------------------------------------------
# PART 9. POWER CALCULATION BY SIMULATION
# ------------------------------------------------------------------
# Power is the probability that a study detects a true effect.
# We estimate power in four steps:
#   1. Simulate many studies with a known effect.
#   2. Analyse each study.
#   3. Count the studies that detect the effect.
#   4. Divide that count by the number of studies.
# We aim for power of 80% or more.
#
# We estimate power for two effects in every study:
#   The level change (the jump at the policy start).
#   The slope change (the extra change per year after the policy start).
#
# Speed: brms is too slow for thousands of fits in class.
# We use a fast frequentist model for the main power grids.
# The frequentist model uses the same terms as the Bayesian model.
# It declares an effect "detected" when p < 0.05 and the rate falls.
# Part 10 checks one design with brms.
#
# Four design features drive power:
#   Events:  the number of preterm births in each segment.
#   Periods: the number of periods before and after the policy start.
#   Effect:  the size of the level change and of the slope change.
#   Noise:   the random variation beyond binomial chance.

power_glm <- function(n_pre, n_post, rr, slope_rr, births_scale = 1, n_sim = 100) {
  detected <- replicate(n_sim, {
    s <- simulate_its(n_pre, n_post, rr, slope_rr, births_scale = births_scale)
    m <- glm(cbind(preterm, births - preterm) ~ province + time_c + policy +
               time_after + sin1 + cos1, family = quasibinomial(), data = s)
    est <- summary(m)$coefficients
    c(level = est["policy", "Estimate"] < 0 & est["policy", "Pr(>|t|)"] < 0.05,
      slope = est["time_after", "Estimate"] < 0 & est["time_after", "Pr(>|t|)"] < 0.05)
  })
  rowMeans(detected)   # one power value for the level and one for the slope
}

# ---- Helpers ----------------------------------------------------------
# Three effect scenarios. Each scenario has a level change and a slope change.
scenarios <- data.frame(size     = c("Small", "Medium", "Large"),
                        rr       = c(0.90, 0.80, 0.70),
                        slope_rr = c(0.95, 0.90, 0.85))
scenarios$scenario <- sprintf("%s (level -%d%%, slope -%d%% per year)",
                              scenarios$size,
                              round((1 - scenarios$rr) * 100),
                              round((1 - scenarios$slope_rr) * 100))
scenarios$scenario <- factor(scenarios$scenario, levels = scenarios$scenario)

# This function runs power_glm() on every row of a grid.
run_grid <- function(grid) {
  res <- t(mapply(power_glm, n_pre = grid$n_pre, n_post = grid$n_post,
                  rr = grid$rr, slope_rr = grid$slope_rr,
                  births_scale = grid$births_scale))
  cbind(grid, res)
}

# This function stacks the level power and the slope power for plotting.
to_long <- function(grid) {
  grid %>%
    pivot_longer(c(level, slope), names_to = "type", values_to = "power") %>%
    mutate(type = recode(type, level = "Level change", slope = "Slope change"))
}

# ---- Grid 1: window length and effect size --------------------------
# Question: how many periods do we need before and after the policy start?
grid1 <- merge(expand.grid(periods_each_side = c(6, 13, 26, 52, 78),
                           size = scenarios$size, stringsAsFactors = FALSE),
               scenarios)
grid1$n_pre <- grid1$periods_each_side
grid1$n_post <- grid1$periods_each_side
grid1$births_scale <- 1
grid1_long <- to_long(run_grid(grid1))

cat("\nPOWER by window length (periods on each side) and scenario\n")
grid1_long %>%
  select(periods_each_side, scenario, type, power) %>%
  pivot_wider(names_from = scenario, values_from = power) %>%
  arrange(type, periods_each_side) %>%
  print()

# The solid line shows the simulation.
ggplot(grid1_long, aes(x = periods_each_side, y = power, colour = scenario)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ type) +
  scale_x_continuous(breaks = c(6, 13, 26, 52, 78)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Power by window length and effect size",
       x = "Periods before and after the policy start (13 periods = 1 year)",
       y = "Power", colour = "True effect") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

# ---- Grid 2: births per period and effect size ----------------------
# Question: how many births per period do we need?
# Total births per period across all provinces = sum(provinces$births_mean) x births_scale.
grid2 <- merge(expand.grid(births_scale = c(0.02, 0.05, 0.10, 0.25, 1, 2),
                           size = scenarios$size, stringsAsFactors = FALSE),
               scenarios)
grid2$n_pre <- 26
grid2$n_post <- 26
grid2$births_per_period <- round(sum(provinces$births_mean) * grid2$births_scale)
grid2_long <- to_long(run_grid(grid2))

cat("\nPOWER by births per period and scenario, 26 periods on each side\n")
grid2_long %>%
  select(births_per_period, scenario, type, power) %>%
  pivot_wider(names_from = scenario, values_from = power) %>%
  arrange(type, births_per_period) %>%
  print()

ggplot(grid2_long, aes(x = births_per_period, y = power, colour = scenario)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ type) +
  scale_x_log10(breaks = unique(grid2$births_per_period)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Power by births per period",
       x = "Births per period, all provinces together (log scale)",
       y = "Power", colour = "True effect") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

# ---- Grid 3: unequal windows ----------------------------------------
# Question: does it matter if we have more periods before than after?
# We use the medium scenario: level -20%, slope -10% per year.
grid3 <- expand.grid(n_pre = c(6, 13, 26, 52, 78, 104), n_post = c(6, 13, 26, 52, 78, 104))
grid3$rr <- 0.80
grid3$slope_rr <- 0.90
grid3$births_scale <- 1
grid3_long <- to_long(run_grid(grid3))

cat("\nPOWER for the medium scenario with unequal windows\n")
grid3_long %>%
  pivot_wider(names_from = n_post, values_from = power,
              id_cols = c(type, n_pre)) %>%
  print()

# The heat map shows power in each cell
ggplot(grid3_long, aes(x = factor(n_post), y = factor(n_pre), fill = power)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", power),
                colour = ifelse(power > 0.6, "black", "white")), size = 5) +
  facet_wrap(~ type) +
  scale_colour_identity() +
  scale_fill_viridis_c(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Power with unequal windows (medium scenario)",
       x = "Post-policy periods", y = "Pre-policy periods", fill = "Power")

# ------------------------------------------------------------------
# PART 10. POWER CHECK WITH brms (OPTIONAL, ABOUT 2-3 MINUTES)
# ------------------------------------------------------------------
# We reuse the compiled model. We refit it on 20 new simulated studies.
# A study detects an effect when its 95% credible interval lies fully below the ROPE.
# This rule is stricter than p < 0.05. It demands a practically relevant decrease.
# Bayesian power can therefore sit below the frequentist power in Part 9.
# Twenty studies give a rough estimate only. Use more studies at home.

run_bayes_power <- TRUE
if (run_bayes_power) {
  # The model reports odds ratios. We convert the true level change to match.
  truth <- data.frame(type  = c("Level change", "Slope change"),
                      ratio = c(exp(qlogis(0.07 * 0.80) - qlogis(0.07)), 0.90))

  studies <- lapply(1:20, function(i) {
    df_new <- simulate_its(n_pre, n_post, rr = 0.80, slope_rr = 0.90)
    f <- tryCatch(update(fit, newdata = df_new, recompile = FALSE,
                         seed = i, refresh = 0),
                  error = function(e) NULL)
    if (is.null(f)) return(NULL)
    est <- fixef(f)[c("policy", "time_after"), ]
    data.frame(study = i,
               type  = c("Level change", "Slope change"),
               ratio = exp(est[, "Estimate"]),
               lower = exp(est[, "Q2.5"]),
               upper = exp(est[, "Q97.5"]),
               row.names = NULL)
  })
  studies <- bind_rows(studies) %>%
    mutate(detected = upper < rope_ratio[1])

  cat("\nBayesian power (interval fully below the ROPE):\n")
  studies %>%
    group_by(type) %>%
    summarise(studies = n(), power = mean(detected)) %>%
    print()

  # Each line shows the 95% credible interval of one simulated study.
  # Blue lines detect the effect. Grey lines miss it.
  # The grey band marks the ROPE. The red line marks the true value.
  ggplot(studies, aes(x = ratio, y = factor(study), colour = detected)) +
    annotate("rect", xmin = rope_ratio[1], xmax = rope_ratio[2],
             ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5) +
    geom_vline(data = truth, aes(xintercept = ratio), colour = "red") +
    geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.3) +
    facet_wrap(~ type, scales = "free_x") +
    scale_colour_manual(values = c("TRUE" = "#648fff", "FALSE" = "grey60")) +
    labs(title = "Bayesian power check: one credible interval per study",
         subtitle = "Red: true ratio. Grey band: ROPE.",
         x = "Ratio (odds ratio of the level change; ratio per year of the slope change)",
         y = "Simulated study", colour = "Effect detected")
}

# ------------------------------------------------------------------
# PART 11. HOW TO CHOOSE A DESIGN
# ------------------------------------------------------------------
# 1. Fix the smallest level change and slope change worth detecting.
# 2. Use at least 13 periods before the policy start, to estimate seasons.
# 3. Use the grids to find the periods and births giving 80% power.
# 4. Expect the slope change to need far more post-policy periods than the level change.
# 5. Rerun the grids with a larger noise_sd and a larger rho.
#    Real data are usually noisier than simulated data.
# 6. Plan for a control series if other events might mimic the policy.
# 7. Fit the final model with MCMC (algorithm = "sampling").

# ------------------------------------------------------------------
# PART 12. EXERCISES
# ------------------------------------------------------------------
# 1. Change rr to 0.90 in simulate_its(). Refit the model. What changes?
# 2. Set rho = 0.7 and noise_sd = 0.10. Check the autocorrelation plot.
# 3. Set slope_rr = 1 (no slope change). Does the model still report one?
# 4. Set n_post = 13. Does the slope change stay detectable?
# 5. Widen the ROPE to c(0.90, 1.10). Which conclusions change?
# 6. Compare ADVI with MCMC: set algorithm = "sampling" and compare intervals.
