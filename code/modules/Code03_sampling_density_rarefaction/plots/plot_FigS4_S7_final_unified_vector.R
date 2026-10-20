#!/usr/bin/env Rscript

# Final formatting-only rebuild of Supplementary Figs. S4-S7.
# This script reads completed bin-level summaries, intervals, contrasts, and
# coverage tables. It does not run rarefaction, Monte Carlo, age perturbation,
# EV reconstruction, or any other numerical analysis.
#
# Display smoothing is limited to npreg::ss(spar = 0.2) fitted to the already
# completed bin-level median trajectories. Binwise 2.5th and 97.5th percentile
# boundaries are plotted without smoothing.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(npreg)
  library(grid)
  library(gridExtra)
  library(svglite)
})

options(warn = 1)

args <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)])
component_root <- normalizePath(file.path(dirname(normalizePath(script_file)), ".."))
work_root <- file.path(component_root, "plot_snapshot_inputs")
results_root <- file.path(work_root, "03_results")
output_dir <- file.path(
  component_root, "outputs", "plots", "Supplement_final_unified_figures"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Locked display and export parameters
# ---------------------------------------------------------------------------

font_family <- "sans"
spar_value <- 0.2
kpg_age <- 66.016
major_age_breaks <- seq(65.0, 67.0, by = 0.5)
minor_age_breaks <- seq(64.8, 67.2, by = 0.1)

col_observed <- "grey45"
col_median <- "#2B6F9F"
col_interval <- "#AFC6D9"
col_kpg <- "#7B3F28"
col_regression <- "#B44B4B"
col_points <- "#4C78A8"

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

formal_age_scale <- function() {
  scale_x_reverse(
    limits = c(67.2, 64.8),
    breaks = major_age_breaks,
    labels = function(x) sprintf("%.1f", x),
    minor_breaks = minor_age_breaks,
    expand = expansion(mult = c(0, 0)),
    guide = guide_axis(minor.ticks = TRUE)
  )
}

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
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold", colour = "black", hjust = 0, size = 8.8
    ),
    panel.spacing.x = unit(1.05, "lines"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.1),
    legend.key.width = unit(1.05, "lines"),
    plot.margin = margin(4, 8, 4, 6)
  )

title_grob <- function(text) {
  textGrob(
    text, x = unit(6, "pt"), y = unit(0.48, "npc"),
    just = c("left", "center"),
    gp = gpar(fontfamily = font_family, fontsize = 12.2, fontface = "bold")
  )
}

panel_header_grob <- function(letter, subtitle, method_line = NULL) {
  grobTree(
    textGrob(
      letter, x = unit(6, "pt"), y = unit(0.69, "npc"),
      just = c("left", "center"),
      gp = gpar(fontfamily = font_family, fontsize = 11.2, fontface = "bold")
    ),
    textGrob(
      subtitle, x = unit(27, "pt"), y = unit(0.69, "npc"),
      just = c("left", "center"),
      gp = gpar(fontfamily = font_family, fontsize = 9.6)
    ),
    if (!is.null(method_line)) {
      textGrob(
        method_line, x = unit(6, "pt"), y = unit(0.19, "npc"),
        just = c("left", "center"),
        gp = gpar(fontfamily = font_family, fontsize = 8.2, col = "grey28")
      )
    } else {
      nullGrob()
    }
  )
}

align_plot_grobs <- function(...) {
  grobs <- list(...)
  widths <- do.call(unit.pmax, lapply(grobs, function(x) x$widths))
  lapply(grobs, function(x) {
    x$widths <- widths
    x
  })
}

extract_legend <- function(plot_object) {
  gt <- ggplotGrob(plot_object)
  idx <- which(grepl("^guide-box", gt$layout$name))
  if (!length(idx)) return(nullGrob())
  candidates <- gt$grobs[idx]
  nonempty <- which(vapply(
    candidates,
    function(x) !inherits(x, "zeroGrob") && length(x$grobs) > 0,
    logical(1)
  ))
  if (!length(nonempty)) return(nullGrob())
  candidates[[nonempty[1]]]
}

draw_grob <- function(grob_object) {
  grid.newpage()
  grid.draw(grob_object)
}

export_vector_set <- function(stem, grob_object, width, height) {
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(output_dir, paste0(stem, ".svg"))
  png_path <- file.path(output_dir, paste0(stem, ".png"))

  grDevices::pdf(
    file = pdf_path, width = width, height = height,
    family = "Helvetica", useDingbats = FALSE,
    compress = TRUE, onefile = TRUE
  )
  draw_grob(grob_object)
  grDevices::dev.off()

  svglite::svglite(
    file = svg_path, width = width, height = height,
    system_fonts = list(sans = "Arial"), fix_text_size = FALSE
  )
  draw_grob(grob_object)
  grDevices::dev.off()

  ragg::agg_png(
    filename = png_path, width = width, height = height,
    units = "in", res = 600, background = "white",
    scaling = 1
  )
  draw_grob(grob_object)
  grDevices::dev.off()
}

as_logical_eligible <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) == "true"
}

fit_bin_medians <- function(data, x_name, y_name, group_name = NULL) {
  fit_one <- function(d) {
    d <- d %>% arrange(.data[[x_name]])
    model <- npreg::ss(
      x = d[[x_name]], y = d[[y_name]], spar = spar_value
    )
    pred <- predict(
      model,
      x = seq(min(d[[x_name]]), max(d[[x_name]]), length.out = 481),
      se.fit = FALSE
    )
    tibble(
      Age_Ma = as.numeric(pred$x),
      npreg_fitted = as.numeric(pred$y)
    )
  }
  if (is.null(group_name)) return(fit_one(data))
  data %>%
    group_split(.data[[group_name]]) %>%
    lapply(function(d) {
      group_value <- unique(d[[group_name]])
      out <- fit_one(d)
      out[[group_name]] <- group_value
      out
    }) %>%
    bind_rows()
}

# ---------------------------------------------------------------------------
# Fig. S4: EV-sampling intensity associations
# ---------------------------------------------------------------------------

density_dir <- file.path(results_root, "stage3_sampling_density")
inventory <- read_csv(
  file.path(density_dir, "audited_50kyr_sampling_composition_table.csv"),
  show_col_types = FALSE
) %>% arrange(desc(bin_center_Ma))
density_stats <- read_csv(
  file.path(density_dir, "density_associations_moving_block_bootstrap.csv"),
  show_col_types = FALSE
) %>%
  filter(method == "Spearman", block_length_bins == 4)

predictor_map <- tibble(
  predictor = c(
    "contributing_operational_record_count",
    "raw_measurement_count_contributing_record_bins_primary",
    "raw_measurement_count_all_retained_secondary",
    "contributing_canonical_site_count"
  ),
  source_column = c(
    "contributing_record_count",
    "raw_measurement_count_contributing_record_bins",
    "raw_measurement_count_all_retained",
    "contributing_canonical_site_count"
  ),
  raw_label = c(
    "Contributing records", "Contributing measurements",
    "All retained measurements", "Contributing sites"
  ),
  difference_label = c(
    "Change in contributing records", "Change in contributing measurements",
    "Change in all retained measurements", "Change in contributing sites"
  )
)

build_density_plot_data <- function(scale_key) {
  bind_rows(lapply(seq_len(nrow(predictor_map)), function(i) {
    key <- predictor_map$predictor[i]
    source <- predictor_map$source_column[i]
    x <- inventory[[source]]
    y <- inventory$observed_global_EV
    if (scale_key == "first_difference_older_to_younger") {
      x <- diff(x)
      y <- diff(y)
      panel <- predictor_map$difference_label[i]
    } else {
      panel <- predictor_map$raw_label[i]
    }
    tibble(
      predictor = key, panel = panel, x = x, y = y,
      scale = scale_key
    )
  }))
}

build_density_annotation <- function(scale_key) {
  label_column <- if (scale_key == "raw_bins") "raw_label" else "difference_label"
  density_stats %>%
    filter(scale == scale_key) %>%
    inner_join(
      predictor_map %>% select(predictor, panel = all_of(label_column)),
      by = "predictor"
    ) %>%
    transmute(
      panel,
      rho_label = sprintf("rho == %.2f", estimate),
      ci_label = sprintf("95%% CI [%.2f, %.2f]", CI_2_5, CI_97_5)
    )
}

raw_density_data <- build_density_plot_data("raw_bins")
diff_density_data <- build_density_plot_data("first_difference_older_to_younger")
raw_density_annotation <- build_density_annotation("raw_bins")
diff_density_annotation <- build_density_annotation(
  "first_difference_older_to_younger"
)

density_theme <- formal_theme +
  theme(
    legend.position = "none",
    plot.title = element_text(
      size = 8.7, face = "bold", hjust = 0, margin = margin(b = 2)
    ),
    axis.title.x = element_text(size = 7.8, margin = margin(t = 4)),
    axis.title.y = element_text(size = 7.8, margin = margin(r = 4)),
    plot.margin = margin(3, 5, 5, 5)
  )

build_density_small_plot <- function(
  data, annotation, panel_name, x_label, y_label,
  first_difference = FALSE
) {
  d <- data %>% filter(panel == panel_name)
  a <- annotation %>% filter(panel == panel_name)
  p <- ggplot(d, aes(x = x, y = y))
  if (first_difference) {
    p <- p +
      geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.38) +
      geom_vline(xintercept = 0, colour = "grey82", linewidth = 0.38)
  }
  p <- p +
    geom_point(size = 1.15, colour = col_points, alpha = 0.70) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      colour = col_regression, linewidth = 0.58
    ) +
    geom_text(
      data = a,
      aes(x = Inf, y = Inf, label = rho_label),
      inherit.aes = FALSE, hjust = 1.04, vjust = 1.15,
      size = 2.30, parse = TRUE
    ) +
    geom_text(
      data = a,
      aes(x = Inf, y = Inf, label = ci_label),
      inherit.aes = FALSE, hjust = 1.04, vjust = 2.85,
      size = 2.30
    )
  if (first_difference) {
    p <- p + scale_y_continuous(
      limits = c(-0.26, 0.26), breaks = c(-0.2, 0, 0.2),
      expand = expansion(mult = c(0.01, 0.02))
    )
  } else {
    p <- p + scale_y_continuous(
      limits = c(0.20, 0.71), breaks = seq(0.2, 0.7, by = 0.1),
      expand = expansion(mult = c(0.01, 0.02))
    )
  }
  p +
    labs(x = x_label, y = y_label) +
    density_theme
}

raw_x_labels <- c(
  "Contributing records", "Contributing measurements",
  "All retained measurements", "Contributing sites"
)
diff_x_labels <- c(
  "Change in contributing records",
  "Change in contributing\nmeasurements",
  "Change in all retained\nmeasurements",
  "Change in contributing sites"
)

fig_s4_raw_plots <- lapply(seq_along(predictor_map$raw_label), function(i) {
  build_density_small_plot(
    raw_density_data, raw_density_annotation,
    predictor_map$raw_label[i], raw_x_labels[i], "Global EV",
    first_difference = FALSE
  )
})
fig_s4_diff_plots <- lapply(
  seq_along(predictor_map$difference_label),
  function(i) {
    build_density_small_plot(
      diff_density_data, diff_density_annotation,
      predictor_map$difference_label[i], diff_x_labels[i],
      "Change in Global EV",
      first_difference = TRUE
    )
  }
)

s4_individual_grobs <- do.call(
  align_plot_grobs,
  c(
    lapply(fig_s4_raw_plots, ggplotGrob),
    lapply(fig_s4_diff_plots, ggplotGrob)
  )
)
s4_raw_row <- arrangeGrob(grobs = s4_individual_grobs[1:4], nrow = 1)
s4_diff_row <- arrangeGrob(grobs = s4_individual_grobs[5:8], nrow = 1)
fig_s4 <- arrangeGrob(
  title_grob(
    "Associations between Global environmental variability and sampling intensity"
  ),
  panel_header_grob("A", "Raw-bin associations"),
  s4_raw_row,
  panel_header_grob("B", "First-difference associations"),
  s4_diff_row,
  ncol = 1,
  heights = c(0.45, 0.34, 2.50, 0.34, 2.50)
)

# ---------------------------------------------------------------------------
# Fig. S5: sampling-density sensitivity
# ---------------------------------------------------------------------------

a1 <- read_csv(
  file.path(density_dir, "A1_fixed_record_rarefaction_bin_summary.csv"),
  show_col_types = FALSE
) %>% mutate(eligible = as_logical_eligible(eligible))
a2 <- read_csv(
  file.path(
    density_dir, "A2_measurement_depth_sensitivity",
    "A2_measurement_depth_sensitivity_bin_summary.csv"
  ),
  show_col_types = FALSE
) %>% mutate(eligible = as_logical_eligible(eligible))
a2_contrasts <- read_csv(
  file.path(
    density_dir, "A2_measurement_depth_sensitivity",
    "A2_measurement_depth_sensitivity_contrast_summary.csv"
  ),
  show_col_types = FALSE
)

stopifnot(
  setequal(unique(a1$N_records), c(3, 5, 7)),
  setequal(unique(a2$m_measurements), c(3, 4, 5)),
  all(a1$iterations == 10000),
  all(a2$iterations == 10000),
  sum(a2$eligible[a2$m_measurements == 5]) == 47,
  identical(
    sort(a2$bin_center_Ma[a2$m_measurements == 5 & !a2$eligible]),
    c(67.10, 67.15)
  )
)

panel_levels_a1 <- c(
  "3 records per bin", "5 records per bin", "7 records per bin"
)
panel_levels_a2 <- c(
  "3 measurements per record", "4 measurements per record",
  "5 measurements per record"
)
series_colours <- c(
  "N = 3" = "#2B8CBE", "N = 5" = "#8856A7", "N = 7" = "#E6550D",
  "m = 3" = "#2B8CBE", "m = 4" = "#8856A7", "m = 5" = "#E6550D"
)

a1_plot_data <- a1 %>%
  mutate(
    series = paste0("N = ", N_records),
    panel = factor(paste(N_records, "records per bin"), levels = panel_levels_a1)
  ) %>%
  arrange(N_records, bin_center_Ma)
a2_plot_data <- a2 %>%
  mutate(
    series = paste0("m = ", m_measurements),
    panel = factor(
      paste(m_measurements, "measurements per record"),
      levels = panel_levels_a2
    )
  ) %>%
  arrange(m_measurements, bin_center_Ma)

fit_rarefaction_medians <- function(data, group_name) {
  data %>%
    filter(eligible) %>%
    group_split(.data[[group_name]]) %>%
    lapply(function(d) {
      group_value <- unique(d[[group_name]])
      model <- npreg::ss(
        x = d$bin_center_Ma, y = d$median_rarefied_EV, spar = spar_value
      )
      pred <- predict(
        model,
        x = seq(min(d$bin_center_Ma), max(d$bin_center_Ma), length.out = 401),
        se.fit = FALSE
      )
      out <- tibble(
        bin_center_Ma = as.numeric(pred$x),
        npreg_fitted = as.numeric(pred$y)
      )
      out[[group_name]] <- group_value
      out
    }) %>%
    bind_rows()
}

a1_fits <- fit_rarefaction_medians(a1_plot_data, "N_records") %>%
  mutate(
    panel = factor(paste(N_records, "records per bin"), levels = panel_levels_a1)
  )
a2_fits <- fit_rarefaction_medians(a2_plot_data, "m_measurements") %>%
  mutate(
    panel = factor(
      paste(m_measurements, "measurements per record"),
      levels = panel_levels_a2
    ),
    npreg_fitted = if_else(
      m_measurements == 5 &
        bin_center_Ma >= 67.075 & bin_center_Ma <= 67.175,
      NA_real_, npreg_fitted
    )
  )

rarefaction_layers <- function(data, fits) {
  list(
    geom_rect(
      data = phase_background,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
      inherit.aes = FALSE, alpha = 0.23, colour = NA
    ),
    geom_ribbon(
      data = data,
      aes(
        x = bin_center_Ma, ymin = p2_5, ymax = p97_5, group = series
      ),
      inherit.aes = FALSE, fill = "grey70", alpha = 0.34,
      colour = NA, na.rm = TRUE
    ),
    geom_line(
      data = data,
      aes(
        x = bin_center_Ma, y = median_rarefied_EV,
        colour = series, group = series
      ),
      linewidth = 0.38, na.rm = TRUE
    ),
    geom_point(
      data = data %>% filter(eligible),
      aes(x = bin_center_Ma, y = median_rarefied_EV, colour = series),
      size = 1.05
    ),
    geom_line(
      data = fits,
      aes(x = bin_center_Ma, y = npreg_fitted),
      inherit.aes = FALSE, colour = "black", linewidth = 0.90,
      na.rm = TRUE
    ),
    geom_vline(
      xintercept = kpg_age, linetype = "dashed",
      colour = col_kpg, linewidth = 0.45
    ),
    scale_fill_manual(
      values = setNames(phase_background$fill, phase_background$phase),
      guide = "none"
    )
  )
}

build_rarefaction_triptych <- function(data, fits, panel_levels) {
  phase_labels <- phase_background %>%
    mutate(panel = factor(panel_levels[1], levels = panel_levels))
  kpg_label <- tibble(
    bin_center_Ma = kpg_age - 0.015, y = 1.385, label = "K-Pg",
    panel = factor(panel_levels[1], levels = panel_levels)
  )

  ggplot() +
    rarefaction_layers(data, fits) +
    geom_text(
      data = phase_labels,
      aes(x = midpoint, y = Inf, label = label_short),
      inherit.aes = FALSE, vjust = 1.35,
      colour = "grey38", size = 2.10
    ) +
    geom_text(
      data = kpg_label,
      aes(x = bin_center_Ma, y = y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      colour = col_kpg, size = 2.25
    ) +
    facet_wrap(~panel, nrow = 1) +
    scale_colour_manual(values = series_colours, guide = "none") +
    formal_age_scale() +
    scale_y_continuous(
      limits = c(0, 1.4), breaks = seq(0, 1.4, by = 0.2),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Age (Ma)", y = "Rarefied Global EV") +
    formal_theme
}

fig_s5a_time <- build_rarefaction_triptych(
  a1_plot_data, a1_fits, panel_levels_a1
)
fig_s5b_time <- build_rarefaction_triptych(
  a2_plot_data, a2_fits, panel_levels_a2
)

contrast_labels <- c(
  phase1_minus_pre = "Phase 1 - pre-escalation",
  phase2_minus_pre = "Phase 2 - pre-escalation",
  recovery_minus_phase2 = "Recovery - Phase 2"
)
a2_contrast_plot_data <- a2_contrasts %>%
  filter(contrast %in% names(contrast_labels)) %>%
  mutate(
    contrast = factor(contrast, levels = rev(names(contrast_labels))),
    contrast_label = factor(
      unname(contrast_labels[as.character(contrast)]),
      levels = unname(contrast_labels[rev(names(contrast_labels))])
    ),
    series = paste0("m = ", m_measurements)
  )

fig_s5b_contrasts <- ggplot(
  a2_contrast_plot_data,
  aes(x = median, y = contrast_label, colour = series)
) +
  geom_vline(
    xintercept = 0, linetype = "dashed",
    colour = "grey40", linewidth = 0.45
  ) +
  geom_errorbarh(
    aes(xmin = p2_5, xmax = p97_5),
    position = position_dodge(width = 0.55),
    height = 0.16, linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 1.90) +
  scale_colour_manual(values = series_colours, name = NULL) +
  labs(
    x = "Contrast in EV",
    y = NULL,
    subtitle = "Phase contrasts (95% rarefaction interval)"
  ) +
  formal_theme +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8.5, colour = "grey25"),
    plot.margin = margin(4, 8, 5, 6)
  )

s5_time_grobs <- align_plot_grobs(
  ggplotGrob(fig_s5a_time), ggplotGrob(fig_s5b_time)
)
fig_s5 <- arrangeGrob(
  title_grob(
    "Sampling-density sensitivity of Global environmental variability reconstruction"
  ),
  panel_header_grob(
    "A", "Fixed-record rarefaction",
    paste0(
      "10,000 rarefaction iterations; 95% rarefaction interval; ",
      "npreg, spar = 0.2; npreg fitted to bin-level median values; ",
      "interval boundaries unsmoothed"
    )
  ),
  s5_time_grobs[[1]],
  panel_header_grob(
    "B", "Joint record-number and measurement-depth rarefaction",
    paste0(
      "5 records per bin; 10,000 rarefaction iterations; ",
      "95% rarefaction interval; npreg, spar = 0.2; ",
      "npreg fitted to bin-level median values; interval boundaries unsmoothed"
    )
  ),
  s5_time_grobs[[2]],
  ggplotGrob(fig_s5b_contrasts),
  ncol = 1,
  heights = c(0.47, 0.55, 3.65, 0.55, 3.65, 2.35)
)

# ---------------------------------------------------------------------------
# Figs. S6-S7: simulation summaries with observed-bin coverage
# ---------------------------------------------------------------------------

inventory_dir <- file.path(results_root, "EV_uncertainty_inventory")

read_measurement_simulation <- function(which_panel) {
  if (which_panel == "documented") {
    input_dir <- file.path(
      inventory_dir,
      "Layer1A_documented_measurement_uncertainty_formal_50kyr_10000"
    )
    summary_path <- file.path(
      input_dir, "Layer1A_formal_50kyr_Global_EV_summary.csv"
    )
    coverage_path <- file.path(
      input_dir, "Layer1A_formal_perturbation_coverage_50kyr.csv"
    )
  } else {
    input_dir <- file.path(
      inventory_dir,
      "Layer2_conservative_measurement_sensitivity_formal_50kyr_10000"
    )
    summary_path <- file.path(
      input_dir, "Layer2_formal_50kyr_Global_EV_summary.csv"
    )
    coverage_path <- file.path(
      input_dir, "Layer2_formal_perturbation_coverage_50kyr.csv"
    )
  }

  ev <- read_csv(summary_path, show_col_types = FALSE) %>%
    filter(iterations == 10000) %>%
    transmute(
      Age_Ma = as.numeric(midpoint_Ma),
      observed_EV = as.numeric(observed_EV),
      simulation_median = as.numeric(MC_median),
      simulation_p2_5 = as.numeric(.data[["MC_2.5th"]]),
      simulation_p97_5 = as.numeric(.data[["MC_97.5th"]])
    ) %>%
    arrange(Age_Ma)
  coverage <- read_csv(coverage_path, show_col_types = FALSE) %>%
    transmute(
      Age_Ma = as.numeric(midpoint_Ma),
      coverage = as.numeric(perturbation_coverage_proportion)
    ) %>%
    arrange(Age_Ma)

  stopifnot(
    nrow(ev) == 49L, nrow(coverage) == 49L,
    max(abs(ev$Age_Ma - coverage$Age_Ma)) < 1e-10,
    all(ev$simulation_p2_5 <= ev$simulation_median),
    all(ev$simulation_median <= ev$simulation_p97_5)
  )
  list(ev = ev, coverage = coverage)
}

read_chronology_simulation <- function(which_panel) {
  fig_dir <- file.path(
    work_root, "04_figures", "Chronological_sensitivity_Fig_S7"
  )
  if (which_panel == "documented") {
    plot_path <- file.path(
      fig_dir, "Fig_S7_chronological_sensitivity_plot_data.csv"
    )
    coverage_path <- file.path(
      fig_dir, "Fig_S7A_chronology_perturbation_coverage.csv"
    )
  } else {
    plot_path <- file.path(
      fig_dir,
      "Fig_S7B_conservative_method_based_chronological_sensitivity_plot_data.csv"
    )
    coverage_path <- file.path(
      fig_dir, "Fig_S7B_chronology_perturbation_coverage.csv"
    )
  }
  ev <- read_csv(plot_path, show_col_types = FALSE) %>%
    transmute(
      Age_Ma = as.numeric(Age_Ma),
      observed_EV = as.numeric(observed_EV),
      simulation_median = as.numeric(simulation_median),
      simulation_p2_5 = as.numeric(simulation_p2_5),
      simulation_p97_5 = as.numeric(simulation_p97_5)
    ) %>%
    arrange(Age_Ma)
  coverage <- read_csv(coverage_path, show_col_types = FALSE) %>%
    transmute(
      Age_Ma = as.numeric(Age_Ma),
      coverage = as.numeric(coverage)
    ) %>%
    arrange(Age_Ma)
  stopifnot(
    nrow(ev) == 49L, nrow(coverage) == 49L,
    max(abs(ev$Age_Ma - coverage$Age_Ma)) < 1e-10,
    all(ev$simulation_p2_5 <= ev$simulation_median),
    all(ev$simulation_median <= ev$simulation_p97_5)
  )
  list(ev = ev, coverage = coverage)
}

build_simulation_main <- function(ev) {
  fit <- fit_bin_medians(
    ev, "Age_Ma", "simulation_median", group_name = NULL
  )
  ev_long <- bind_rows(
    ev %>%
      transmute(Age_Ma, value = observed_EV, series = "Observed Global EV"),
    ev %>%
      transmute(Age_Ma, value = simulation_median, series = "Simulation median")
  ) %>%
    mutate(
      series = factor(
        series, levels = c("Observed Global EV", "Simulation median")
      )
    )

  ggplot() +
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
      aes(
        x = Age_Ma, y = npreg_fitted,
        colour = "npreg, spar = 0.2"
      ),
      linewidth = 0.90
    ) +
    geom_vline(
      xintercept = kpg_age, linetype = "dashed",
      colour = col_kpg, linewidth = 0.45
    ) +
    geom_text(
      data = phase_background,
      aes(x = midpoint, y = 0.888, label = label_short),
      inherit.aes = FALSE, colour = "grey38", size = 2.45, vjust = 1
    ) +
    annotate(
      "text", x = kpg_age - 0.015, y = 0.892, label = "K-Pg",
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
    formal_age_scale() +
    scale_y_continuous(
      limits = c(0.15, 0.90), breaks = seq(0.2, 0.9, by = 0.1),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = "Global EV (mean record-level IQR)") +
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(
          linewidth = c(0.45, 0.65, 0.95), shape = c(16, NA, NA)
        )
      ),
      fill = guide_legend(order = 2)
    ) +
    formal_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.minor.ticks.x.bottom = element_blank(),
      plot.margin = margin(4, 8, 1, 6)
    )
}

build_coverage_panel <- function(coverage, y_label) {
  ggplot(coverage, aes(x = Age_Ma, y = coverage)) +
    geom_area(fill = col_interval, alpha = 0.72) +
    geom_line(colour = col_median, linewidth = 0.55) +
    geom_vline(
      xintercept = kpg_age, linetype = "dashed",
      colour = col_kpg, linewidth = 0.45
    ) +
    formal_age_scale() +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
      labels = function(x) paste0(round(100 * x), "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Age (Ma)", y = y_label) +
    formal_theme +
    theme(
      legend.position = "none",
      plot.margin = margin(1, 8, 5, 6)
    )
}

build_simulation_figure <- function(
  total_title, panel_a_subtitle, panel_b_subtitle,
  inputs_a, inputs_b, coverage_label
) {
  main_a_with_legend <- build_simulation_main(inputs_a$ev)
  main_b_with_legend <- build_simulation_main(inputs_b$ev)
  legend <- extract_legend(main_a_with_legend)
  main_a <- main_a_with_legend + theme(legend.position = "none")
  main_b <- main_b_with_legend + theme(legend.position = "none")
  coverage_a <- build_coverage_panel(inputs_a$coverage, coverage_label)
  coverage_b <- build_coverage_panel(inputs_b$coverage, coverage_label)

  aligned <- align_plot_grobs(
    ggplotGrob(main_a), ggplotGrob(coverage_a),
    ggplotGrob(main_b), ggplotGrob(coverage_b)
  )

  arrangeGrob(
    title_grob(total_title),
    panel_header_grob(
      "A", panel_a_subtitle,
      paste0(
        "10,000 simulations; 95% simulation interval; npreg, spar = 0.2; ",
        "npreg fitted to bin-level median values; interval boundaries unsmoothed"
      )
    ),
    aligned[[1]], aligned[[2]],
    panel_header_grob(
      "B", panel_b_subtitle,
      paste0(
        "10,000 simulations; 95% simulation interval; npreg, spar = 0.2; ",
        "npreg fitted to bin-level median values; interval boundaries unsmoothed"
      )
    ),
    aligned[[3]], aligned[[4]],
    legend,
    ncol = 1,
    heights = c(0.47, 0.58, 3.20, 0.95, 0.58, 3.20, 0.95, 0.42)
  )
}

measurement_a <- read_measurement_simulation("documented")
measurement_b <- read_measurement_simulation("conservative")
stopifnot(
  max(abs(
    measurement_a$ev$observed_EV - measurement_b$ev$observed_EV
  )) < 1e-10
)
fig_s6 <- build_simulation_figure(
  "Measurement-uncertainty sensitivity of Global environmental variability reconstruction",
  "Documented measurement uncertainty",
  "Conservative method-based measurement uncertainty",
  measurement_a, measurement_b,
  "Measurement perturbation\ncoverage (%)"
)

chronology_a <- read_chronology_simulation("documented")
chronology_b <- read_chronology_simulation("conservative")
stopifnot(
  max(abs(
    chronology_a$ev$observed_EV - chronology_b$ev$observed_EV
  )) < 1e-10
)
fig_s7 <- build_simulation_figure(
  "Chronological-uncertainty sensitivity of Global environmental variability reconstruction",
  "Documented chronological uncertainty",
  "Conservative method-based chronological uncertainty",
  chronology_a, chronology_b,
  "Chronology perturbation\ncoverage (%)"
)

# ---------------------------------------------------------------------------
# Editable vector and high-resolution exports
# ---------------------------------------------------------------------------

figure_arg <- sub(
  "^--figure=", "", args[grepl("^--figure=", args)]
)
requested_figure <- if (length(figure_arg)) toupper(figure_arg[1]) else "ALL"
stopifnot(requested_figure %in% c("ALL", "S4", "S5", "S6", "S7"))

if (requested_figure %in% c("ALL", "S4")) {
  export_vector_set("Fig_S4_final", fig_s4, width = 10.35, height = 6.25)
}
if (requested_figure %in% c("ALL", "S5")) {
  export_vector_set("Fig_S5_final", fig_s5, width = 10.35, height = 12.20)
}
if (requested_figure %in% c("ALL", "S6")) {
  export_vector_set("Fig_S6_final", fig_s6, width = 10.35, height = 11.15)
}
if (requested_figure %in% c("ALL", "S7")) {
  export_vector_set("Fig_S7_final", fig_s7, width = 10.35, height = 11.15)
}

if (requested_figure == "ALL") {
  message("Created final unified Supplementary Figs. S4-S7.")
} else {
  message("Created final unified Supplementary Fig. ", requested_figure, ".")
}
message("SVG and PDF outputs are vector exports with live text; no raster embedding was used.")
message("K-Pg boundary = 66.016 Ma in Figs. S5-S7; Fig. S4 has no K-Pg line.")
message("No rarefaction, Monte Carlo, uncertainty analysis, Dataset, manuscript, or caption was modified.")
