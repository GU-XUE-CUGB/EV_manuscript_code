#!/usr/bin/env Rscript

# Fig. S8 drawing only. This script reads the completed formal 10,000-
# simulation joint summary and does not run uncertainty propagation, EV
# reconstruction, or any other numerical analysis. The 95% simulation interval
# is plotted from the unsmoothed binwise percentiles. Display smoothing is
# npreg::ss(spar = 0.2) fitted only to the 49 bin-level simulation medians.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(npreg)
  library(grid)
  library(gridExtra)
  library(svglite)
})

options(warn = 1)

args <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)])
work_root <- normalizePath(file.path(dirname(normalizePath(script_file)), ".."))
input_dir <- file.path(
  work_root, "reference_figures", "03_results", "EV_uncertainty_inventory",
  "Joint_conservative_measurement_chronology_sensitivity_formal_50kyr_10000"
)
output_dir <- file.path(
  work_root, "04_figures",
  "Joint_conservative_measurement_chronology_Fig_S8"
)
unified_dir <- file.path(
  work_root, "04_figures", "Supplement_final_unified_figures"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(unified_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(
  input_dir, "Joint_BB_formal_50kyr_Global_EV_summary.csv"
)
stopifnot(file.exists(input_path))

font_family <- "sans"
spar_value <- 0.2
kpg_age <- 66.016
major_age_breaks <- seq(65.0, 67.0, by = 0.5)
minor_age_breaks <- seq(64.8, 67.2, by = 0.1)

col_observed <- "grey45"
col_median <- "#2B6F9F"
col_interval <- "#AFC6D9"
col_kpg <- "#7B3F28"

phase_background <- tibble(
  phase = factor(
    c("Pre-escalation", "Phase 1", "Phase 2", "Recovery"),
    levels = c("Pre-escalation", "Phase 1", "Phase 2", "Recovery")
  ),
  xmin = c(66.925, 66.425, 66.025, 64.800),
  xmax = c(67.200, 66.925, 66.425, 65.975),
  midpoint = c(67.075, 66.675, 66.225, 65.375),
  label_short = c("Pre", "Phase 1", "Phase 2", "Recovery"),
  fill = c("#EDF4F8", "#F3F3F3", "#FDF1E6", "#F2F2F2")
)

ev <- read_csv(input_path, show_col_types = FALSE) %>%
  transmute(
    Age_Ma = as.numeric(midpoint_Ma),
    observed_EV = as.numeric(observed_EV),
    simulation_median = as.numeric(simulation_median),
    simulation_p2_5 = as.numeric(.data[["simulation_2.5th"]]),
    simulation_p97_5 = as.numeric(.data[["simulation_97.5th"]]),
    interval_width = as.numeric(interval_width),
    valid_simulations = as.integer(valid_simulations)
  ) %>%
  arrange(Age_Ma)

stopifnot(
  nrow(ev) == 49L,
  all(ev$valid_simulations == 10000L),
  all(is.finite(as.matrix(ev[, 1:6]))),
  all(ev$simulation_p2_5 <= ev$simulation_median),
  all(ev$simulation_median <= ev$simulation_p97_5),
  max(ev$simulation_p97_5) < 0.95
)

fit_model <- npreg::ss(
  x = ev$Age_Ma, y = ev$simulation_median, spar = spar_value
)
fit_prediction <- predict(
  fit_model,
  x = seq(min(ev$Age_Ma), max(ev$Age_Ma), length.out = 481),
  se.fit = FALSE
)
fit <- tibble(
  Age_Ma = as.numeric(fit_prediction$x),
  npreg_fitted = as.numeric(fit_prediction$y)
)

ev_long <- bind_rows(
  ev %>% transmute(Age_Ma, value = observed_EV, series = "Observed Global EV"),
  ev %>% transmute(Age_Ma, value = simulation_median, series = "Simulation median")
) %>%
  mutate(
    series = factor(
      series, levels = c("Observed Global EV", "Simulation median")
    )
  )

formal_age_scale <- scale_x_reverse(
  limits = c(67.2, 64.8),
  breaks = major_age_breaks,
  labels = function(x) sprintf("%.1f", x),
  minor_breaks = minor_age_breaks,
  expand = expansion(mult = c(0, 0)),
  guide = guide_axis(minor.ticks = TRUE)
)

formal_theme <- theme_classic(base_size = 9.5, base_family = font_family) +
  theme(
    text = element_text(family = font_family, colour = "black"),
    axis.title = element_text(size = 9.3, colour = "black"),
    axis.text = element_text(size = 8.2, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    axis.ticks = element_line(colour = "black", linewidth = 0.42),
    axis.minor.ticks.x.bottom = element_line(colour = "black", linewidth = 0.32),
    axis.ticks.length = unit(4.0, "pt"),
    axis.minor.ticks.length = unit(2.1, "pt"),
    panel.grid.major.y = element_line(colour = "grey89", linewidth = 0.32),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.1),
    legend.key.width = unit(1.05, "lines"),
    plot.margin = margin(4, 8, 5, 6)
  )

plot_object <- ggplot() +
  geom_rect(
    data = phase_background,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
    inherit.aes = FALSE, alpha = 0.23, colour = NA
  ) +
  geom_ribbon(
    data = ev,
    aes(
      x = Age_Ma, ymin = simulation_p2_5, ymax = simulation_p97_5,
      fill = "95% simulation interval"
    ),
    inherit.aes = FALSE, alpha = 0.48, colour = NA
  ) +
  geom_line(
    data = ev_long,
    aes(
      x = Age_Ma, y = value, colour = series,
      linewidth = series, group = series
    )
  ) +
  geom_point(
    data = ev,
    aes(x = Age_Ma, y = observed_EV, colour = "Observed Global EV"),
    size = 1.18, stroke = 0
  ) +
  geom_line(
    data = fit,
    aes(x = Age_Ma, y = npreg_fitted, colour = "npreg, spar = 0.2"),
    linewidth = 0.90
  ) +
  geom_vline(
    xintercept = kpg_age, linetype = "dashed",
    colour = col_kpg, linewidth = 0.45
  ) +
  geom_text(
    data = phase_background,
    aes(x = midpoint, y = 0.938, label = label_short),
    inherit.aes = FALSE, colour = "grey38", size = 2.45, vjust = 1
  ) +
  annotate(
    "text", x = kpg_age - 0.015, y = 0.942, label = "K-Pg",
    colour = col_kpg, size = 2.40, hjust = 0
  ) +
  scale_fill_manual(
    values = c(
      setNames(phase_background$fill, phase_background$phase),
      "95% simulation interval" = col_interval
    ),
    breaks = "95% simulation interval", name = NULL
  ) +
  scale_colour_manual(
    values = c(
      "Observed Global EV" = col_observed,
      "Simulation median" = col_median,
      "npreg, spar = 0.2" = "black"
    ),
    breaks = c(
      "Observed Global EV", "Simulation median", "npreg, spar = 0.2"
    ),
    name = NULL
  ) +
  scale_linewidth_manual(
    values = c("Observed Global EV" = 0.30, "Simulation median" = 0.55),
    guide = "none"
  ) +
  formal_age_scale +
  scale_y_continuous(
    limits = c(0.15, 0.95), breaks = seq(0.2, 0.9, by = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Age (Ma)", y = "Global EV (mean record-level IQR)") +
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        linewidth = c(0.45, 0.65, 0.95), shape = c(16, NA, NA)
      )
    ),
    fill = guide_legend(order = 2)
  ) +
  formal_theme

title_grob <- textGrob(
  "Joint sensitivity of Global environmental variability reconstruction to conservative measurement and chronological uncertainty",
  x = unit(6, "pt"), y = unit(0.48, "npc"),
  just = c("left", "center"),
  gp = gpar(fontfamily = font_family, fontsize = 12.2, fontface = "bold")
)
method_grob <- textGrob(
  paste0(
    "10,000 simulations; 95% simulation interval; npreg, spar = 0.2; ",
    "npreg fitted to 49 bin-level simulation medians; interval boundaries unsmoothed"
  ),
  x = unit(6, "pt"), y = unit(0.52, "npc"),
  just = c("left", "center"),
  gp = gpar(fontfamily = font_family, fontsize = 8.2, col = "grey28")
)

figure_grob <- arrangeGrob(
  title_grob, method_grob, ggplotGrob(plot_object),
  ncol = 1, heights = c(0.48, 0.33, 4.15)
)

export_vector_set <- function(directory, stem, grob_object, width, height) {
  pdf_path <- file.path(directory, paste0(stem, ".pdf"))
  svg_path <- file.path(directory, paste0(stem, ".svg"))
  png_path <- file.path(directory, paste0(stem, ".png"))

  grDevices::pdf(
    file = pdf_path, width = width, height = height,
    family = "Helvetica", useDingbats = FALSE,
    compress = TRUE, onefile = TRUE
  )
  grid.newpage(); grid.draw(grob_object); grDevices::dev.off()

  svglite::svglite(
    file = svg_path, width = width, height = height,
    system_fonts = list(sans = "Arial"), fix_text_size = FALSE
  )
  grid.newpage(); grid.draw(grob_object); grDevices::dev.off()

  ragg::agg_png(
    filename = png_path, width = width, height = height,
    units = "in", res = 600, background = "white", scaling = 1
  )
  grid.newpage(); grid.draw(grob_object); grDevices::dev.off()
}

write_csv(ev, file.path(output_dir, "Fig_S8_joint_formal_plot_data.csv"))
write_csv(fit, file.path(output_dir, "Fig_S8_joint_npreg_predictions.csv"))
write_csv(
  tibble(
    formal_name = "Joint conservative measurement and chronological uncertainty",
    simulations = 10000L,
    interval = "binwise 2.5th-97.5th percentiles; unsmoothed",
    smoothing = "npreg fitted to 49 bin-level simulation medians",
    spar = spar_value,
    kpg_boundary_Ma = kpg_age,
    coverage_panel_included = FALSE,
    source = input_path
  ),
  file.path(output_dir, "Fig_S8_joint_metadata.csv")
)

export_vector_set(output_dir, "Fig_S8_joint_conservative_measurement_chronological_uncertainty", figure_grob, 10.35, 5.45)
export_vector_set(unified_dir, "Fig_S8_final", figure_grob, 10.35, 5.45)

message("Created Fig. S8 from the completed formal 10,000-simulation joint summary.")
message("No simulation, EV reconstruction, Dataset, manuscript, Supplement, or Fig. S4-S7 file was modified.")
