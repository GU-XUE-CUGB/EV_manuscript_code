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
work_root <- normalizePath(file.path(dirname(normalizePath(script_file)), ".."))
results_root <- file.path(work_root, "reference_figures", "03_results")
output_dir <- file.path(
  work_root, "04_figures", "Supplement_final_unified_figures"
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


export_vector_set("Fig_S6_final", fig_s6, width = 10.35, height = 11.15)
