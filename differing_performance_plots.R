############################################################
# Title: finding where performance differs between two treatments
#
# Description: Analyze data to 
#              find periods using moving window 
#              statistics where performance differs.
############################################################

# 1.0 librarys #############
library(tidyverse)
library(lubridate)
library(ggpubr)


# 2.0 simulated date REPLACE WITH ACTUAL DATA #####################
set.seed(123)
n_days <- 1200
dates <- seq.Date(from = as.Date("2015-01-01"), by = "day", length.out = n_days)

# Simulate filter run times
filt1_runtime <- c(rnorm(200, 300, 20), rnorm(400, 250, 10), rnorm(600, 280, 15))
filt2_runtime <- c(rnorm(200, 310, 25), rnorm(400, 255, 12), rnorm(600, 285, 18))

JDK <- data.frame(
  date = dates,
  filt1_runtime = filt1_runtime,
  filt2_runtime = filt2_runtime
)

# Reshape for plotting
runtime_data <- JDK %>%
  pivot_longer(cols = c(filt1_runtime, filt2_runtime),
               names_to = "Filter", values_to = "RunTime") %>%
  mutate(Filter = ifelse(Filter=="filt1_runtime","F1","F2"))

# 1. Define 6-month windows
start_date <- min(JDK$date)
end_date <- max(JDK$date)
window_starts <- seq(start_date, end_date, by = "6 months")
window_ends <- c(window_starts[-1]-1, end_date)

windows <- data.frame(start=window_starts[-length(window_starts)], end=window_ends[-length(window_ends)])

# 2. Test F1 vs F2 in each window
windows$p_value <- NA

for(i in 1:nrow(windows)){
  f1 <- JDK %>% filter(date >= windows$start[i] & date <= windows$end[i]) %>% pull(filt1_runtime)
  f2 <- JDK %>% filter(date >= windows$start[i] & date <= windows$end[i]) %>% pull(filt2_runtime)
  windows$p_value[i] <- t.test(f1, f2)$p.value
}

# 3. Keep windows where p < 0.05
diff_windows <- windows %>% filter(p_value < 0.05)

# 4. Plot with shaded windows
runtime_plot <- ggplot(runtime_data, aes(x=date, y=RunTime, color=Filter)) +
  geom_line() +
  geom_point(size=1) +
  theme_bw() +
  scale_color_manual(values=c("lightskyblue","deepskyblue4")) +
  ggtitle("Filter Run Time with Difference Windows") +
  geom_rect(data=diff_windows, inherit.aes=FALSE,
            aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
            fill="lightblue", alpha=0.75)

runtime_plot
