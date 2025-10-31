##############################################################
# Example Workflow: Mixed-Effects Model to Test Operational Treatment Reccomendations
#
# Description:
#   Demonstrates how to:
#     1. Fit a mixed-effects model using lme4::lmer()
#     2. Simulate bootstrapped confidence intervals via bootMer()
#     3. Generate predicted effluent TSS values under different scenarios
#     4. Visualize predicted vs observed data with uncertainty bounds
#
#   All data here are synthetic and reproducible.
##############################################################

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------

library(tidyverse)
library(lme4)
library(ggplot2)
library(ggpubr)

# ------------------------------------------------------------
# 1. Generate Example Dataset
# REplace with actual data
# ------------------------------------------------------------

set.seed(42)

n <- 150
example_data <- tibble(
  flow = runif(n, 3000, 4000),
  ph_meter = runif(n, 6, 7.5),
  alum_dose = runif(n, 2, 7),
  poly_dose = runif(n, 1, 2.5),
  total_solids_clarifer = runif(n, 200, 1000)
) %>%
  mutate(
    eff_tss = 20 + 0.004 * flow + 
      -3 * (ph_meter - 6.8) +
      0.8 * alum_dose +
      0.5 * poly_dose +
      rnorm(n, 0, 5),
    flow_binned = ifelse(flow < 3515, "<3515", ">3515"),
    ph = case_when(
      ph_meter < 6.5 ~ "<6.5",
      ph_meter >= 7 ~ ">7",
      TRUE ~ "6.5-7"
    )
  )

# ------------------------------------------------------------
# 2. Fit Mixed-Effects Model
# ------------------------------------------------------------

model <- lmer(eff_tss ~ alum_dose + poly_dose + total_solids_clarifer + ph_meter +
                flow + (1 | flow_binned) + (1 | ph),
              data = example_data)

summary(model)

# ------------------------------------------------------------
# 3. Create New Scenario Data
# again replace with actual data and actual recomendations
# ------------------------------------------------------------

new_data <- expand_grid(
  flow = seq(3000, 4000, length.out = 50),
  ph_meter = seq(6, 7.5, length.out = 10)
) %>%
  mutate(
    alum_dose = ifelse(flow < 3515, 4, 5),
    poly_dose = 1.5,
    total_solids_clarifer = ifelse(flow < 3515, 500, 1000),
    flow_binned = ifelse(flow < 3515, "<3515", ">3515"),
    ph = case_when(
      ph_meter < 6.5 ~ "<6.5",
      ph_meter >= 7 ~ ">7",
      TRUE ~ "6.5-7"
    )
  )

# ------------------------------------------------------------
# 4. Define Bootstrapping Functions
# ------------------------------------------------------------

predfn <- function(.) predict(., newdata = new_data, re.form = NULL)

sumBoot <- function(merBoot) {
  data.frame(
    fit = apply(merBoot$t, 2, function(x) quantile(x, 0.5, na.rm = TRUE)),
    lwr = apply(merBoot$t, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
    upr = apply(merBoot$t, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  )
}

# ------------------------------------------------------------
# 5. Bootstrapping Model Predictions
# ------------------------------------------------------------

boot_res <- bootMer(model, predfn, nsim = 100, use.u = TRUE, type = "parametric")

boot_summary <- sumBoot(boot_res)

# Combine bootstrap summary with scenario data
results <- bind_cols(new_data, boot_summary)

# ------------------------------------------------------------
# 6. Plot Bootstrapped Predictions
# ------------------------------------------------------------

ggplot(results, aes(x = flow, y = fit, color = ph)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = ph), alpha = 0.2, color = NA) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Bootstrapped Mixed-Effects Model Predictions",
    x = "Influent Flow (m³/h)",
    y = "Predicted Effluent TSS (mg/L)",
    color = "pH Category",
    fill = "pH Category"
  )

# ------------------------------------------------------------
# 7. Example Comparison: Observed vs Predicted
# ------------------------------------------------------------

# Simulate observed data (for example purposes, use actual data)
results$Actual_Effluent_TSS <- results$fit + rnorm(nrow(results), 0, 3)
results$Target <- ifelse(results$fit < 40, "Predicted <40", "Predicted ≥40")

library(ggplot2)
library(dplyr)

# Assume 'results' dataframe contains:
# Date (as.Date), fit, lwr, upr, Target ("Predicted <25"/"<40"), 
# and `Actual Effluent TSS` (numeric)

# Example color scheme
Value <- c("Observed" = "black", "Predicted <25" = "firebrick4", "Predicted <40" = "lightskyblue")

# Format dates for clean x-axis labels

# simuate it here for ex, use actual data

results <- results %>%
  dplyr::mutate(
    Date = seq.Date(from = as.Date("2024-05-01"), by = "day", length.out = n()),
    DateLabel = format(Date, "%m-%d"),
    DateLabel = factor(DateLabel, levels = unique(DateLabel))
  )

# results <- results %>%
#   dplyr::mutate(DateLabel = format(Date, "%m-%d")) %>%
#   dplyr::mutate(DateLabel = factor(DateLabel, levels = format(sort(unique(Date)), "%m-%d")))

# Position adjustment to dodge overlapping bars
dodge <- position_dodge(width = 0.8)

# --- Plot ---
ggplot(results, aes(x = DateLabel, y = fit, color = Target)) +
  geom_errorbar(
    aes(ymin = lwr, ymax = upr, color = Target),
    width = 0.6, size = 0.8, position = dodge
  ) +
  geom_point(aes(y = `Actual_Effluent_TSS`, color = "Observed"), size = 1.5) +
  scale_color_manual(values = Value) +
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0)) +
  labs(
    title = NULL,
    x = "Sampling Date",
    y = "Effluent TSS (mg/L)",
    color = "Effluent TSS"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6)
  )

# ------------------------------------------------------------
# 8. Export CI Data
# ------------------------------------------------------------

# write.csv(results, "output/mixed_model_bootstrap_results.csv", row.names = FALSE)
