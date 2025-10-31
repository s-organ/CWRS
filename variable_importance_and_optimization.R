##############################################################
# Example Analysis Pipeline: Variable Importance and Logistic Model
# Description: 
#   Example workflow for identifying important predictors, 
#   fitting a logistic regression, and visualizing results.
##############################################################

# ------------------------------------------------------------
# 0.0 Load Libraries
# ------------------------------------------------------------

# Note: Not all of these are needed for this example; 
# adjust to your own workflow.
library(tidyverse)
library(dplyr)
library(stringr)
library(readxl)
library(plyr)
library(glmnet)
library(MASS)
library(ggplot2)
library(tidyr)
library(corrplot)
library(sjPlot)
library(ggpubr)
library(car)
library(forecast)
library(tseries) 
library(astsa)
library(ggdist)
library(vars)
library(orcutt)
library(prais)
library(reshape2)
library(mlbench)
library(rsample)
library(GGally)
library(ggeffects)
library(glmm)
library(caret)
library(lme4)
library(rsm)
library(VGAM)
library(plotrix)
library(rpart)

# ------------------------------------------------------------
# 1.0 Import Data
# ------------------------------------------------------------

# Replace with your own dataset
# Example: data <- read_excel("path/to/your/data.xlsx")

# For demonstration purposes, create a mock dataset:
data <- data.frame(
  flow = runif(100, 500, 5000),
  total_solids_clarifer = runif(100, 50, 200),
  ph_meter = rnorm(100, 7, 0.5),
  temp_meter = rnorm(100, 20, 2),
  alum_dose = runif(100, 10, 50),
  poly_dose = runif(100, 1, 10),
  eff_tss = rnorm(100, 30, 10)
)

# Create binary response variables
data <- data %>%
  mutate(low_tss_40 = ifelse(eff_tss < 40, 1, 0),
         low_tss_25 = ifelse(eff_tss < 25, 1, 0)) %>%
  drop_na(flow, low_tss_40)

# ------------------------------------------------------------
# 2.0 Variable Importance (Decision Tree)
# ------------------------------------------------------------

# Fit a simple tree model
fit <- rpart(
  eff_tss ~ flow + total_solids_clarifer + ph_meter + temp_meter + alum_dose + poly_dose,
  data = data
)

# Extract and visualize variable importance
df <- data.frame(imp = fit$variable.importance) %>%
  tibble::rownames_to_column("variable") %>%
  arrange(imp) %>%
  mutate(variable = forcats::fct_inorder(variable))

g1 <- ggplot(df, aes(x = variable, y = imp)) +
  geom_col(fill = "lightskyblue") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Variable Importance (Example)",
    x = "Variable",
    y = "Importance Measure"
  )

print(g1)

# ------------------------------------------------------------
# 3.0 Logistic Regression on the Most Important Variable
# ------------------------------------------------------------

# Example using 'flow' as the key predictor
model <- glm(low_tss_40 ~ flow, family = binomial(link = "logit"), data = data)
summary(model)

# Predict probabilities
data$pred_prob <- predict(model, type = "response")

# ------------------------------------------------------------
# 4.0 Identify Slope Change / Threshold
# ------------------------------------------------------------

# Use LOESS smoothing to identify changes in fitted probabilities
smooth_fit <- loess(data$low_tss_40 ~ data$flow)
dl <- diff(smooth_fit$fitted)

# Find the steepest negative slope (greatest decrease)
largest_slope <- which.min(dl)
bound <- data$flow[largest_slope]

# OR Find the steepest positive slope (greatest increase)
largest_slope <- which.max(dl)
bound <- data$flow[largest_slope]


cat("Estimated threshold (bound):", round(bound, 2), "\n")

# ------------------------------------------------------------
# 5.0 Visualization
# ------------------------------------------------------------

plot_data <- data %>%
  dplyr::select(flow, pred_prob) %>%
  mutate(`TSS Target` = "<40 (mg/L)")

g2 <- ggplot(plot_data, aes(x = flow, y = pred_prob, color = `TSS Target`)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = bound, linetype = 2, color = "black") +
  annotate("text", x = bound + 300, y = 0.25, 
           label = paste("Flow =", round(bound, 2)), color = "black") +
  theme_bw() +
  labs(
    title = "Predicted Probability vs. Flow Rate",
    x = "Flow Rate (m³/h)",
    y = "Predicted Probability of Meeting TSS Target"
  )

print(g2)

# ------------------------------------------------------------
# 6.0 Statistical Test Example
# ------------------------------------------------------------

data <- data %>%
  dplyr::mutate(flow_split = ifelse(flow < bound, "low flow", "high flow"))

test_result <- wilcox.test(
  data$flow[data$flow_split == "low flow"],
  data$flow[data$flow_split == "high flow"]
)

print(test_result)
