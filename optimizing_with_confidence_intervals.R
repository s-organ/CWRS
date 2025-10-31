##############################################################
# Example Workflow: GAM Optimization and Bootstrap Confidence
#
# Description:
#   This example demonstrates how to:
#     1. Generate example process data
#     2. Fit Generalized Additive Models (GAMs)
#     3. Use bootstrapping for confidence intervals
#     4. Visualize response curves for different predictors
#     5. choose optimal bound
#
##############################################################

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------

library(tidyverse)
library(mgcv)
library(boot)
library(ggpubr)
library(lmerTest)
library(caret)
library(rsample)
library(ggdist)
library(VGAM)
library(MASS)
library(glmnet)
library(plotrix)
library(rpart)

# ------------------------------------------------------------
# 1. Helper Functions
# ------------------------------------------------------------

# Normalize numeric vectors between 0 and 1
normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# Consistent ggplot theme
theme_custom <- theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 11)
  )

# Generic bootstrap GAM fit
boot_gam_fit <- function(data, xvar, yvar = "response", nboot = 500, length_out = 100) {
  # Create sequence of x values
  xseq <- seq(min(data[[xvar]], na.rm = TRUE),
              max(data[[xvar]], na.rm = TRUE),
              length.out = length_out)
  
  # Create prediction grid
  grid <- setNames(data.frame(xseq), xvar)
  
  # Bootstrap function
  boot_fn <- function(data, indices) {
    d <- data[indices, ]
    mod <- gam(reformulate(xvar, response = yvar), data = d)
    predict(mod, newdata = grid)
  }
  
  # Run bootstrap
  set.seed(123)
  boot_results <- boot(data, statistic = boot_fn, R = nboot)
  
  # Summarize bootstrap results (mean + 95% CI)
  grid$fit <- colMeans(boot_results$t)
  grid$lwr <- apply(boot_results$t, 2, quantile, probs = 0.025)
  grid$upr <- apply(boot_results$t, 2, quantile, probs = 0.975)
  
  # Return summarized data frame with predictor and intervals
  return(grid)
}

# Generic plotting function
plot_boot <- function(grid, xvar, title_text, ylab_text = "Predicted Response") {
  # Calculate percentiles of fitted values
  p50 <- quantile(grid$fit, 0.5)
  p75 <- quantile(grid$fit, 0.75)
  p90 <- quantile(grid$fit, 0.9)
  
  ggplot(grid, aes(x = .data[[xvar]], y = fit)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), color = "deepskyblue", alpha = 1) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_hline(yintercept = c(p50, p75, p90),
               linetype = "dashed",
               color = c("darkgreen", "goldenrod", "firebrick4")) +
    annotate("text", x = min(grid[[xvar]], na.rm = TRUE),
             y = p50, label = "50th percentile", vjust = -0.5, hjust = 0,
             color = "darkgreen", size = 4, fontface = "bold") +
    annotate("text", x = min(grid[[xvar]], na.rm = TRUE),
             y = p75, label = "75th percentile", vjust = -0.5, hjust = 0,
             color = "goldenrod", size = 4, fontface = "bold") +
    annotate("text", x = min(grid[[xvar]], na.rm = TRUE),
             y = p90, label = "90th percentile", vjust = -0.5, hjust = 0,
             color = "firebrick4", size = 4, fontface = "bold") +
    labs(title = title_text, x = xvar, y = ylab_text) +
    theme_custom
}
# ------------------------------------------------------------
# 2. Generate Example Dataset
# replace with actual data
# ------------------------------------------------------------

set.seed(42)
n <- 200
data_example <- tibble(
  flow_rate    = runif(n, 1000, 5000),
  alum_dose    = runif(n, 5, 25),
  polymer_dose = runif(n, 0.5, 10),
  mix_ph       = rnorm(n, 6, 0.3),
  raw_uv254    = runif(n, 0.05, 0.25)
) %>%
  mutate(
    runtime = 5 +
      0.03 * flow_rate -
      0.2 * alum_dose^2 / 100 +
      0.5 * mix_ph +
      rnorm(n, 0, 3),
    response = normalize(runtime)
  )

# Split dataset into “low” and “high” UV groups
data_low  <- data_example %>% filter(raw_uv254 < 0.15)
data_high <- data_example %>% filter(raw_uv254 >= 0.15)

# ------------------------------------------------------------
# 3. Example Analyses
# ------------------------------------------------------------

## 3.1 Alum dose effect (low UV)
alum_boot <- boot_gam_fit(data_low, "alum_dose", "response")
p1 <- plot_boot(alum_boot, "alum_dose", "Effect of Alum Dose (Low UV)")

## 3.2 pH effect (low UV)
ph_boot <- boot_gam_fit(data_low, "mix_ph", "response")
p2 <- plot_boot(ph_boot, "mix_ph", "Effect of pH (Low UV)")

## 3.3 Polymer effect (low UV)
poly_boot <- boot_gam_fit(data_low, "polymer_dose", "response")
p3 <- plot_boot(poly_boot, "polymer_dose", "Effect of Polymer Dose (Low UV)")

# Combine plots
combined <- ggarrange(p1, p2, p3, ncol = 3)
print(combined)

alum_boot
ph_boot
poly_boot

# -----------------------------------------------------------
# INTERPRETATION: want confidence intervals that are above/below desired percentile
# the x axis value they are associaed with statistically have best response
# ----------------------------------------------------------

# ------------------------------------------------------------
# 4. Model Comparison Example
# REPLACE WITH OWN DATA/OPTIMAL VALUES
# ------------------------------------------------------------

# Compare two GAMs (e.g., “current” vs. “optimized” parameterization)
data_example <- data_example %>%
  mutate(
    mix_ph_adj   = ifelse(raw_uv254 < 0.15, mix_ph - 0.1, mix_ph + 0.1),
    alum_dose_adj = ifelse(raw_uv254 < 0.15, alum_dose - 2, alum_dose + 1)
  )

model_current <- gam(response ~ s(mix_ph) + s(alum_dose) + s(raw_uv254), data = data_example)
model_new     <- gam(response ~ s(mix_ph_adj) + s(alum_dose_adj) + s(raw_uv254), data = data_example)

data_example <- data_example %>%
  mutate(
    fit_current = predict(model_current, type = "response"),
    fit_new     = predict(model_new, type = "response")
  )

ggplot(data_example, aes(x = flow_rate)) +
  geom_line(aes(y = fit_current, color = "Current")) +
  geom_line(aes(y = fit_new, color = "Optimized")) +
  labs(
    title = "Comparison of Model Fits",
    x = "Flow Rate",
    y = "Predicted Response",
    color = "Model"
  ) +
  scale_color_manual(values = c("skyblue", "blue4")) +
  theme_custom
