library(readxl)
library(dplyr)
library(ggplot2)

data_low <- read_excel("IQR_low latitude.xlsx")
colnames(data_low) <- c("IQR", "Midpoint")

set.seed(123)
n_curves <- 10
random_curves <- list()
smooth_curves <- list()

for (i in 1:n_curves) {
  
  sampled_curve <- data_low %>%
    group_by(Midpoint) %>%
    summarise(
      EV = mean(sample(IQR, size = length(IQR) * 2, replace = TRUE)),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    mutate(Curve = paste0("Curve_", i))
  
  random_curves[[i]] <- sampled_curve
  
  ss <- smooth.spline(sampled_curve$Midpoint, sampled_curve$EV, spar = 0.6)
  smooth_curve <- data.frame(
    Midpoint = seq(min(sampled_curve$Midpoint), max(sampled_curve$Midpoint), length.out = 200),
    EV_smooth = predict(ss, seq(min(sampled_curve$Midpoint), max(sampled_curve$Midpoint), length.out = 200))$y,
    Curve = paste0("Curve_", i)
  )
  smooth_curves[[i]] <- smooth_curve
}

boot_results_low <- bind_rows(random_curves)
smooth_results_low <- bind_rows(smooth_curves)

max_ev <- max(c(boot_results_low$EV, smooth_results_low$EV_smooth), na.rm = TRUE)
y_breaks <- seq(0.1, ceiling(max_ev*10)/10, by = 0.2)

ggplot() +
  geom_line(data = boot_results_low,
            aes(x = Midpoint, y = EV, group = Curve, color = Curve),
            size = 0.8, alpha = 0.5) +
  geom_line(data = smooth_results_low,
            aes(x = Midpoint, y = EV_smooth, group = Curve, color = Curve),
            size = 1.1, alpha = 0.9) +
  scale_x_continuous(limits = c(67200, 64800), trans = "reverse") +
  scale_y_continuous(breaks = y_breaks) +
  labs(
    x = "Time (kyr, 67.2–64.8 Ma)",
    y = "Environmental Variability (EV)",
    title = "Upsampled Low-latitude EV curves (mean of 2x IQR per Midpoint)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
