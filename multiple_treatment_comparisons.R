##############################################################
# Example Workflow: Confidence Intervals and Group Comparisons
#
# Description:
#   Demonstrates how to:
#     1. Create synthetic grouped data
#     2. Fit linear models to estimate group means
#     3. Extract and visualize confidence intervals
#     4. Perform ANOVA and pairwise comparisons
#
#   Used for comparing treatments/groups 
##############################################################

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------

library(tidyverse)
library(readxl)    # included for generality, not needed in synthetic example
library(emmeans)
library(sjPlot)
library(ggpubr)

# ------------------------------------------------------------
# 1. Generate Example Dataset
# Replace with actual data
# ------------------------------------------------------------

set.seed(123)

# Example: simulated effluent TSS data under three pH bins
example_data <- tibble(
  ph_meter = runif(150, 6, 7.5),
  eff_tss  = 50 + 
    ifelse(ph_meter < 6.5, rnorm(150, -5, 4),
           ifelse(ph_meter < 7, rnorm(150, 0, 4), rnorm(150, 5, 4)))
) %>%
  mutate(
    ph_bin = case_when(
      ph_meter < 6.5 ~ "<6.5",
      ph_meter >= 6.5 & ph_meter < 7 ~ "6.5-7",
      TRUE ~ ">7"
    )
  )

# ------------------------------------------------------------
# 2. Model and Confidence Intervals
# ------------------------------------------------------------

# Fit model
model <- lm(eff_tss ~ ph_bin, data = example_data)

# Estimated marginal means (with 95% CI)
emm <- emmeans(model, ~ ph_bin)
emm_df <- as.data.frame(emm)

# Pairwise comparisons with Bonferroni correction
comparison_results <- pairs(emm, adjust = "bonferroni")

# Print CI data frame and comparisons
print(emm_df)
print(comparison_results)

# ------------------------------------------------------------
# 3. Plot Estimated Means + CI
# ------------------------------------------------------------

ggplot(emm_df, aes(x = ph_bin, y = emmean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "skyblue4", size = 1) +
  geom_text(aes(label = sprintf("Mean = %.1f", emmean)),
            vjust = -1.0, color = "steelblue4", size = 3.5) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Estimated Group Means with 95% Confidence Intervals",
    x = "pH Group",
    y = "Estimated Mean Effluent TSS"
  )

# Alternative automatic plot
plot_model(model, type = "pred", terms = "ph_bin")

# ------------------------------------------------------------
# 4. ANOVA and Pairwise Tests
# ------------------------------------------------------------

anova_res <- aov(eff_tss ~ ph_bin, data = example_data)
summary(anova_res)

# Boxplot with statistical comparisons
ggplot(example_data, aes(x = ph_bin, y = eff_tss, color = ph_bin)) +
  geom_boxplot() +
  theme_minimal(base_size = 13) +
  stat_compare_means(
    comparisons = list(c("<6.5", "6.5-7"), c("<6.5", ">7"), c("6.5-7", ">7")),
    p.adjust.method = "bonferroni"
  ) +
  labs(title = "Effluent TSS Across pH Groups",
       x = "pH Bin", y = "Effluent TSS (mg/L)")

# Nonparametric pairwise test (Wilcoxon)
pairwise_results <- pairwise.wilcox.test(
  x = example_data$eff_tss,
  g = example_data$ph_bin,
  p.adjust.method = "bonferroni"
)

print(pairwise_results)

# ------------------------------------------------------------
# 5. Save CI Results (Optional)
# ------------------------------------------------------------

# write.csv(emm_df, "output/emmeans_confidence_intervals.csv", row.names = FALSE)