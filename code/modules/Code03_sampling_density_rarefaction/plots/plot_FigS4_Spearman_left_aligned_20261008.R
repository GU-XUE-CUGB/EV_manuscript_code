#!/usr/bin/env Rscript

# Plot Fig. S4 with box frames and bold axis titles.
# Read rho, P, n and block-bootstrap intervals from the input tables.
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
  component_root, "outputs", "plots", "Fig_S4_Spearman_left_aligned_20261008"
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
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    axis.minor.ticks.x.bottom = element_line(colour = "black", linewidth = 0.32),
    axis.ticks.length = unit(4.0, "pt"),
    axis.minor.ticks.length = unit(2.1, "pt"),
    panel.grid.major.y = element_blank(),
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
  heights <- do.call(unit.pmax, lapply(grobs, function(x) x$heights))
  lapply(grobs, function(x) {
    x$widths <- widths
    x$heights <- heights
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


# Reuse the saved rho/P/n table; no correlation or bootstrap is rerun.
nominal_results <- read_csv(file.path(
  work_root, "04_figures", "Fig_S4_S11_style_rectangular_rho_P_20261008",
  "Spearman_rho_nominal_P_existing_CI.csv"), show_col_types = FALSE)
stopifnot(nrow(nominal_results) == 8L, all(is.finite(nominal_results$nominal_P)))

density_stats <- density_stats %>%
  left_join(nominal_results %>% select(predictor, scale, nominal_P),
            by = c("predictor", "scale"))

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
      rho_label = sprintf("\"Spearman's\"~rho == %.2f", estimate),
      p_label = sprintf("P = %.3g", nominal_P)
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
    axis.title.x = element_text(size = 7.8, face = "bold", margin = margin(t = 4)),
    axis.title.y = element_text(size = 7.8, face = "bold", margin = margin(r = 4)),
    plot.margin = margin(3, 5, 5, 5)
  )

build_density_small_plot <- function(
  data, annotation, panel_name, x_label, y_label,
  first_difference = FALSE
) {
  d <- data %>% filter(panel == panel_name)
  a <- annotation %>% filter(panel == panel_name)
  a$label_x <- min(d$x) + 0.005 * diff(range(d$x))
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
      aes(x = label_x, y = Inf, label = rho_label),
      inherit.aes = FALSE, hjust = 0, vjust = 1.15,
      size = 2.30, parse = TRUE
    ) +
    geom_text(
      data = a,
      aes(x = label_x, y = Inf, label = p_label),
      inherit.aes = FALSE, hjust = 0, vjust = 2.85,
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

export_vector_set("Fig_S4_Spearman_left_aligned", fig_s4, width = 12, height = 5)
print(nominal_results %>% select(scale, predictor, n_pairs, rho, nominal_P))
cat("Spearman rho/P left-aligned figure: ", output_dir, "\n", sep = "")


