############################################################
# Title: Ammonia Removal Steady-State Analysis (Example)
#
# Description: Analyze data to 
#              find steady-state periods using moving window 
#              statistics and stationarity tests.
#              Uses simulated data for demonstration.
############################################################

##### 1.0 Libraries #######
library(dplyr)
library(tidyr)
library(tseries)
library(ggplot2)

##### 2.0 Simulate Data #######
# REPLACE WITH ACTUAL DATA
set.seed(123)

n_samples <- 37
dates <- seq.Date(from = as.Date("2025-01-01"), by = "day", length.out = n_samples)

# Simulated ammonia removal (%) with some random noise
ammonia_rem <- c(
  rnorm(15, mean = 50, sd = 10),   # initial unstable period
  rnorm(10, mean = 80, sd = 2),    # steady-state period
  rnorm(12, mean = 75, sd = 5)     # later variation
)

amm_rem_c <- data.frame(
  date = dates,
  ammonia_rem = ammonia_rem
)

##### 3.1 c_effluent Time Series #######
amm_ts_c <- ts(amm_rem_c$ammonia_rem)
plot.ts(amm_ts_c, main = "c_effluent Ammonia Removal Time Series")

##### 3.2 Moving Window Function #######
lag_apply <- function(x, n, callback){
  k <- length(x)
  result <- rep(NA, k)
  for(i in 1:(k - n + 1)){
    result[i] <- callback(x[i:(i + n - 1)])
  }
  return(result)
}

window_size <- 5
moving_sd <- lag_apply(amm_ts_c, window_size, sd)
moving_avg <- lag_apply(amm_ts_c, window_size, mean)

##### 3.3 Prepare Data for Plotting #######
plot_dat <- data.frame(
  Value = c(moving_sd, moving_avg),
  Date = rep(amm_rem_c$date, 2),
  Measure = rep(c("Moving SD", "Moving Average"), each = length(amm_ts_c)),
  Sample = rep(1:length(amm_ts_c), 2)
)

##### 3.4 Plot Moving Statistics #######
ggplot(plot_dat, aes(x = Sample, y = Value, color = Measure)) +
  geom_point() +
  geom_line(linetype = 2) +
  theme_minimal() +
  geom_hline(yintercept = mean(moving_sd, na.rm = TRUE), color = "black", linetype = 2) +
  scale_x_continuous("Sample", labels = as.character(plot_dat$Sample[1:length(amm_ts_c)]), breaks = 1:length(amm_ts_c)) +
  ggtitle("Moving SD and Moving Average of c_effluent Ammonia Removal (Simulated Data)")

# look for where moving average starts to stablize

##### 3.5 Stationarity Check #######
# Steady-state visually around sample 16-25
steady_start <- 16
steady_end <- 25

# Augmented Dickey-Fuller Test for steady-state period
adf_test_result <- adf.test(amm_ts_c[steady_start:steady_end])
print(adf_test_result)

# Check p-values if extending steady-state
p_values <- sapply(steady_end:length(amm_ts_c), function(i){
  tseries::adf.test(amm_ts_c[steady_start:i])$p.value
})
print(p_values)

##### 3.6 Plot with Steady-State Highlighted #######
ggplot(plot_dat, aes(x = Date, y = Value, color = Measure)) +
  geom_point(size = 2) +
  geom_line(linetype = 2, size = 1) +
  theme_minimal() +
  geom_hline(yintercept = mean(moving_sd, na.rm = TRUE), color = "black", linetype = 2, size = 1) +
  geom_vline(xintercept = as.numeric(plot_dat$Date[c(steady_start, steady_end)]),
             linetype = 4, color = "red", size = 1) +
  scale_color_manual(values = c("lightskyblue", "goldenrod")) +
  ggtitle("Steady-State Range for c_effluent Ammonia Removal (Simulated Data)")

