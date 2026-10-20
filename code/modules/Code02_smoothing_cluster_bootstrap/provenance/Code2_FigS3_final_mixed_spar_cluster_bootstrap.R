#!/usr/bin/env Rscript

# Fig. S3 with bin-width-specific npreg smoothing.
# Read 50- and 100-kyr fits and envelopes from the supplied tables.
# Use spar=0.30 for the 150- and 200-kyr fits and calculate their envelopes
# using 5,000 record-cluster resamples, seed 20260820.
# Exclude zero-IQR record-bins; use 2.5th and 97.5th percentile bounds.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(npreg)
  library(svglite)
})

options(warn = 1)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
working_root <- dirname(script_dir)

previous_fig_dir <- file.path(
  working_root, "04_figures", "Fig_S3_npreg_cluster_bootstrap_revised"
)
comparison_dir <- file.path(
  working_root, "04_figures", "Fig_S3_spar_comparison_150_200"
)
output_dir <- file.path(
  working_root, "04_figures", "Fig_S3_final_mixed_spar_cluster_bootstrap"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

observed_file <- file.path(
  previous_fig_dir, "Fig_S3_revised_observed_Global_EV.csv"
)
previous_envelope_file <- file.path(
  previous_fig_dir, "Fig_S3_revised_npreg_bootstrap_envelopes.csv"
)
previous_metadata_file <- file.path(
  previous_fig_dir, "Fig_S3_revised_bootstrap_metadata.csv"
)
spar_comparison_file <- file.path(
  comparison_dir, "Fig_S3_spar_020_025_030_fitted_values_150_200kyr.csv"
)
record_iqr_file <- file.path(
  working_root, "01_data_revised",
  "record_level_IQR_all_bin_widths_revised_v4_n2.csv"
)

required_files <- c(
  observed_file, previous_envelope_file, previous_metadata_file,
  spar_comparison_file, record_iqr_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required input file(s) missing: ", paste(missing_files, collapse = "; "))
}

bootstrap_iterations <- 5000L
random_seed <- 20260820L
prediction_grid_points <- 481L

schemes <- list(
  list(
    id = "200", panel = "A", spar = 0.30,
    title = "69-63 Ma (time bin=200 kyr)",
    xlim = c(69, 63), breaks = seq(63, 69, by = 1)
  ),
  list(
    id = "150", panel = "B", spar = 0.30,
    title = "68-64 Ma (time bin=150 kyr)",
    xlim = c(68, 64), breaks = seq(64, 68, by = 1)
  ),
  list(
    id = "100", panel = "C", spar = 0.20,
    title = "67.2-64.8 Ma (time bin=100 kyr)",
    xlim = c(67.2, 64.8), breaks = seq(65, 67, by = 0.5)
  ),
  list(
    id = "50", panel = "D", spar = 0.20,
    title = "67.2-64.8 Ma (time bin=50 kyr)",
    xlim = c(67.2, 64.8), breaks = seq(65, 67, by = 0.5)
  )
)

expected_bin_counts <- c(`200` = 29L, `150` = 27L, `100` = 25L, `50` = 49L)

observed_all <- read_csv(observed_file, show_col_types = FALSE)
previous_envelopes <- read_csv(previous_envelope_file, show_col_types = FALSE)
previous_metadata <- read_csv(previous_metadata_file, show_col_types = FALSE)
spar_comparison <- read_csv(spar_comparison_file, show_col_types = FALSE)
record_iqr_all <- read_csv(record_iqr_file, show_col_types = FALSE)

required_record_columns <- c(
  "scheme", "alignment", "record_id", "bin_center_Ma", "IQR"
)
if (!all(required_record_columns %in% names(record_iqr_all))) {
  stop("Record-IQR source is missing required fields")
}

observed_by_scheme <- list()
final_envelopes_by_scheme <- list()
final_metadata_by_scheme <- list()

# First, preserve the completed 50/100-kyr outputs verbatim.
for (scheme_id in c("100", "50")) {
  spec <- schemes[[match(scheme_id, vapply(schemes, `[[`, character(1), "id"))]]
  observed <- observed_all |>
    filter(bin_width_kyr == as.integer(scheme_id)) |>
    arrange(Age_Ma)
  existing <- previous_envelopes |>
    filter(bin_width_kyr == as.integer(scheme_id)) |>
    arrange(Age_Ma)
  metadata <- previous_metadata |>
    filter(bin_width_kyr == as.integer(scheme_id))

  if (nrow(observed) != expected_bin_counts[[scheme_id]]) {
    stop("Unexpected observed bin count for reused ", scheme_id, "-kyr panel")
  }
  if (nrow(existing) != prediction_grid_points || anyDuplicated(existing$Age_Ma)) {
    stop("Unexpected reused envelope grid for ", scheme_id, "-kyr panel")
  }
  if (nrow(metadata) != 1L || any(existing$spar != 0.2)) {
    stop("Reused ", scheme_id, "-kyr source is not the confirmed spar=0.20 result")
  }

  final_envelopes_by_scheme[[scheme_id]] <- existing |>
    transmute(
      panel = spec$panel,
      bin_width_kyr = as.integer(scheme_id),
      Age_Ma,
      observed_npreg_fit,
      bootstrap_median,
      bootstrap_p2_5,
      bootstrap_p97_5,
      spar = 0.20,
      calculation_action = "reused_without_recalculation"
    )
  final_metadata_by_scheme[[scheme_id]] <- metadata |>
    mutate(
      final_spar = 0.20,
      calculation_action = "reused_without_recalculation"
    )
  observed_by_scheme[[scheme_id]] <- observed
}

# Recompute only the 150/200-kyr cluster-bootstrap envelopes at spar=0.30.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

for (scheme_id in c("200", "150")) {
  spec <- schemes[[match(scheme_id, vapply(schemes, `[[`, character(1), "id"))]]
  observed <- observed_all |>
    filter(bin_width_kyr == as.integer(scheme_id)) |>
    arrange(Age_Ma)
  selected_fit <- spar_comparison |>
    filter(bin_width_kyr == as.integer(scheme_id)) |>
    arrange(Age_Ma)

  if (nrow(observed) != expected_bin_counts[[scheme_id]]) {
    stop("Unexpected observed bin count for ", scheme_id, "-kyr panel")
  }
  if (nrow(selected_fit) != prediction_grid_points || anyDuplicated(selected_fit$Age_Ma)) {
    stop("Unexpected selected spar=0.30 fit grid for ", scheme_id, "-kyr panel")
  }

  source_all <- record_iqr_all |>
    filter(as.character(scheme) == scheme_id, alignment == "standard") |>
    transmute(
      record_id = as.character(record_id),
      Age_Ma = as.numeric(bin_center_Ma),
      IQR = as.numeric(IQR)
    ) |>
    arrange(record_id, Age_Ma)

  source <- source_all |>
    filter(IQR != 0)

  if (anyDuplicated(source[c("record_id", "Age_Ma")])) {
    stop("Duplicate operational-record/bin IQR rows found for ", scheme_id, " kyr")
  }

  # Confirm that the record-level source exactly reproduces the already locked
  # observed bin-level Global EV values under the zero-IQR exclusion rule.
  source_check <- source |>
    group_by(Age_Ma) |>
    summarise(Global_EV_check = mean(IQR), .groups = "drop") |>
    arrange(Age_Ma)
  if (nrow(source_check) != nrow(observed) ||
      max(abs(source_check$Age_Ma - observed$Age_Ma)) > 1e-12 ||
      max(abs(source_check$Global_EV_check - observed$Global_EV)) > 1e-12) {
    stop("Record-IQR source does not reproduce locked Global EV for ", scheme_id, " kyr")
  }

  prediction_grid_ma <- selected_fit$Age_Ma
  selected_observed_fit <- selected_fit$spar_0_30
  record_ids <- sort(unique(source$record_id))
  n_records <- length(record_ids)
  source$record_index <- match(source$record_id, record_ids)

  bootstrap_predictions <- matrix(
    NA_real_, nrow = prediction_grid_points, ncol = bootstrap_iterations
  )
  available_bin_counts <- integer(bootstrap_iterations)

  # The confirmed seed is reset independently for each bin-width panel, exactly
  # as in the preceding four-panel cluster-bootstrap workflow.
  set.seed(random_seed)

  for (iteration in seq_len(bootstrap_iterations)) {
    sampled_indices <- sample.int(
      n_records, size = n_records, replace = TRUE
    )
    record_weights <- tabulate(sampled_indices, nbins = n_records)

    bootstrap_bins <- source |>
      mutate(cluster_weight = record_weights[record_index]) |>
      filter(cluster_weight > 0) |>
      group_by(Age_Ma) |>
      summarise(
        Global_EV = weighted.mean(IQR, w = cluster_weight),
        .groups = "drop"
      ) |>
      arrange(Age_Ma)

    available_bin_counts[iteration] <- nrow(bootstrap_bins)

    bootstrap_fit <- npreg::ss(
      x = bootstrap_bins$Age_Ma,
      y = bootstrap_bins$Global_EV,
      spar = 0.30
    )
    bootstrap_prediction <- predict(
      bootstrap_fit, x = prediction_grid_ma, se.fit = FALSE
    )
    bootstrap_predictions[, iteration] <- as.numeric(bootstrap_prediction$y)

    if (iteration %% 500L == 0L) {
      message(
        scheme_id, " kyr spar=0.30: completed ", iteration,
        " / ", bootstrap_iterations, " cluster-bootstrap iterations"
      )
    }
  }

  final_envelopes_by_scheme[[scheme_id]] <- data.frame(
    panel = spec$panel,
    bin_width_kyr = as.integer(scheme_id),
    Age_Ma = prediction_grid_ma,
    observed_npreg_fit = as.numeric(selected_observed_fit),
    bootstrap_median = apply(bootstrap_predictions, 1, median, na.rm = TRUE),
    bootstrap_p2_5 = apply(
      bootstrap_predictions, 1, quantile,
      probs = 0.025, na.rm = TRUE, names = FALSE
    ),
    bootstrap_p97_5 = apply(
      bootstrap_predictions, 1, quantile,
      probs = 0.975, na.rm = TRUE, names = FALSE
    ),
    spar = 0.30,
    calculation_action = "recomputed_cluster_bootstrap",
    stringsAsFactors = FALSE
  )

  final_metadata_by_scheme[[scheme_id]] <- data.frame(
    panel = spec$panel,
    bin_width_kyr = as.integer(scheme_id),
    observed_bin_count = nrow(observed),
    operational_record_clusters = n_records,
    nonzero_record_bin_IQR_rows = nrow(source),
    excluded_zero_IQR_rows = nrow(source_all) - nrow(source),
    bootstrap_iterations = bootstrap_iterations,
    random_seed = random_seed,
    spar = 0.30,
    prediction_grid_points = prediction_grid_points,
    iterations_with_all_bins = sum(available_bin_counts == nrow(observed)),
    iterations_with_one_or_more_empty_bins = sum(available_bin_counts < nrow(observed)),
    minimum_available_bins = min(available_bin_counts),
    final_spar = 0.30,
    calculation_action = "recomputed_cluster_bootstrap",
    stringsAsFactors = FALSE
  )
  observed_by_scheme[[scheme_id]] <- observed
}

# Restore submitted Fig. S3 panel order A-D.
scheme_order <- c("200", "150", "100", "50")
final_envelopes <- bind_rows(final_envelopes_by_scheme[scheme_order])
final_metadata <- bind_rows(final_metadata_by_scheme[scheme_order])

write_csv(
  final_envelopes,
  file.path(output_dir, "Fig_S3_final_npreg_bootstrap_envelopes.csv")
)
write_csv(
  final_metadata,
  file.path(output_dir, "Fig_S3_final_bootstrap_metadata.csv")
)

# Preserve the established common y-axis scale unless a final envelope would be
# clipped. All other layout and styling match the current revised Fig. S3.
envelope_min <- min(final_envelopes$bootstrap_p2_5, na.rm = TRUE)
envelope_max <- max(final_envelopes$bootstrap_p97_5, na.rm = TRUE)
y_lower <- if (envelope_min < 0) floor(envelope_min * 10) / 10 else 0
y_upper <- if (envelope_max > 1.1) ceiling(envelope_max * 10) / 10 else 1.1
y_breaks <- seq(ceiling(y_lower * 5) / 5, floor(y_upper * 5) / 5, by = 0.2)

draw_final_fig_s3 <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(
    mfrow = c(2, 2),
    mar = c(4.3, 4.6, 3.1, 1.0),
    oma = c(0.2, 0.2, 0.2, 0.2),
    mgp = c(2.7, 0.75, 0),
    tcl = -0.3,
    family = "sans",
    font.lab = 2,
    cex.lab = 1.15,
    cex.axis = 0.95,
    cex.main = 1.05,
    las = 1
  )

  for (spec in schemes) {
    observed <- observed_by_scheme[[spec$id]]
    envelope <- final_envelopes_by_scheme[[spec$id]]

    plot(
      NA_real_, NA_real_,
      xlim = spec$xlim,
      ylim = c(y_lower, y_upper),
      xaxs = "i", yaxs = "i",
      axes = FALSE,
      xlab = "Age (Ma)", ylab = "Global EV",
      main = spec$title,
      bty = "o"
    )
    axis(1, at = spec$breaks, labels = format(spec$breaks, trim = TRUE))
    axis(2, at = y_breaks, labels = format(y_breaks, nsmall = 1, trim = TRUE))
    box(lwd = 0.8)

    polygon(
      c(envelope$Age_Ma, rev(envelope$Age_Ma)),
      c(envelope$bootstrap_p2_5, rev(envelope$bootstrap_p97_5)),
      col = "grey88", border = NA
    )

    # Submitted/current Fig. S3 has no K-Pg reference line; no new annotation
    # is introduced in this final update.
    lines(observed$Age_Ma, observed$Global_EV, col = "grey62", lwd = 0.7)
    points(observed$Age_Ma, observed$Global_EV,
           pch = 16, col = "grey38", cex = 0.55)
    lines(envelope$Age_Ma, envelope$observed_npreg_fit,
          col = "black", lwd = 2.0)

    mtext(
      spec$panel, side = 3, line = 0.7, adj = -0.15,
      font = 2, cex = 1.15
    )
  }
}

base_name <- file.path(
  output_dir, "Fig_S3_final_mixed_spar_cluster_bootstrap"
)

png(
  paste0(base_name, ".png"),
  width = 4200, height = 3040, res = 400, bg = "white"
)
draw_final_fig_s3()
dev.off()

tiff(
  paste0(base_name, ".tiff"),
  width = 4200, height = 3040, res = 400,
  compression = "lzw", bg = "white"
)
draw_final_fig_s3()
dev.off()

pdf(
  paste0(base_name, ".pdf"),
  width = 10.5, height = 7.6, useDingbats = FALSE
)
draw_final_fig_s3()
dev.off()

svglite::svglite(
  paste0(base_name, ".svg"),
  width = 10.5, height = 7.6, bg = "white"
)
draw_final_fig_s3()
dev.off()

message("Final Fig. S3 created in: ", output_dir)
message("50/100 kyr reused at spar=0.20; 150/200 kyr recomputed at spar=0.30.")
