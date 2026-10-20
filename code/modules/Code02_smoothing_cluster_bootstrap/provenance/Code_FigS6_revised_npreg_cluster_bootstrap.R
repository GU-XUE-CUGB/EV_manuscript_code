#!/usr/bin/env Rscript

# Proxy-panel npreg fits and record-cluster bootstrap.
#
# Observed curves are read directly from the already-calculated 50-kyr
# EV/mean-IQR columns in the formal workbooks. No EV or IQR is recalculated.
# Confidence envelopes use 5,000 operational-record cluster bootstrap samples:
# complete record-level IQR trajectories are resampled with replacement, a
# bootstrap bin mean is formed, and npreg::ss(spar = 0.2) is fitted to each
# trajectory. The plotted center line is always the fit to the observed,
# workbook-supplied bin-level EV series (not the bootstrap median).

suppressPackageStartupMessages({
  library(readxl)
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

formal_dir <- file.path(working_root, "01_data_revised", "formal_datasets")
output_dir <- file.path(
  working_root, "04_figures", "Fig_S6_npreg_cluster_bootstrap_revised"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  global = file.path(formal_dir, "Dataset_S1_revised_v5.xlsx"),
  temperature = file.path(formal_dir, "Dataset_S2_revised.xlsx"),
  d13C = file.path(formal_dir, "Dataset_S3_revised_v3.xlsx"),
  d18O = file.path(formal_dir, "Dataset_S4_revised_v3.xlsx"),
  weathering = file.path(formal_dir, "Dataset_S5_revised_v4.xlsx"),
  pCO2 = file.path(formal_dir, "Dataset_S6_revised_v4.xlsx")
)
if (any(!file.exists(files))) {
  stop("Missing formal revised workbook(s): ", paste(files[!file.exists(files)], collapse = ", "))
}

bootstrap_iterations <- 5000L
random_seed <- 20260820L
spar_value <- 0.2
prediction_grid <- seq(64.8, 67.2, length.out = 481L)

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

read_numeric_pair <- function(path, sheet, value_col, age_col) {
  dat <- suppressMessages(read_excel(path, sheet = sheet, col_names = TRUE))
  tibble(
    Age_Ma = suppressWarnings(as.numeric(dat[[age_col]])) / 1000,
    EV = suppressWarnings(as.numeric(dat[[value_col]]))
  ) |>
    filter(is.finite(Age_Ma), is.finite(EV)) |>
    arrange(Age_Ma)
}

# Directly read the already-calculated 50-kyr bin-level series.
observed <- list(
  global = read_numeric_pair(files[["global"]], "total IQR", 15, 16),
  temperature = read_numeric_pair(files[["temperature"]], "total IQR", 15, 16),
  d18O = read_numeric_pair(files[["d18O"]], "total IQR", 15, 16),
  d13C = read_numeric_pair(files[["d13C"]], "total IQR", 15, 16),
  weathering = read_numeric_pair(files[["weathering"]], "total IQR", 15, 16)
)

if (nrow(observed$global) != 49L) {
  stop("Dataset S1 does not contain the expected 49 observed 50-kyr Global EV bins")
}

# Extract already-calculated 50-kyr record-level IQRs from every operational
# record sheet. Sheet identity is the operational-record cluster identifier.
read_record_sheets <- function(path, category, excluded_sheets) {
  record_sheets <- setdiff(excel_sheets(path), excluded_sheets)
  out <- vector("list", length(record_sheets))
  for (i in seq_along(record_sheets)) {
    sheet_name <- record_sheets[[i]]
    dat <- suppressMessages(read_excel(path, sheet = sheet_name, col_names = TRUE))
    headers <- names(dat)
    iqr_cols <- which(grepl("^IQR", headers))
    if (!length(iqr_cols)) stop("No IQR field in ", basename(path), " / ", sheet_name)
    iqr_col <- max(iqr_cols)
    age_col <- iqr_col + 1L
    if (age_col > ncol(dat)) stop("Missing 50-kyr midpoint field in ", sheet_name)
    out[[i]] <- tibble(
      category = category,
      record_id = paste(category, sheet_name, sep = "::"),
      Age_Ma = suppressWarnings(as.numeric(dat[[age_col]])) / 1000,
      IQR = suppressWarnings(as.numeric(dat[[iqr_col]]))
    ) |>
      filter(is.finite(Age_Ma), is.finite(IQR))
  }
  bind_rows(out)
}

records <- bind_rows(
  read_record_sheets(files[["temperature"]], "temperature", c("total data", "total IQR", "total z-score")),
  read_record_sheets(files[["d18O"]], "d18O", c("total data", "total IQR", "total z-score")),
  read_record_sheets(files[["d13C"]], "d13C", c("total data", "total IQR", "total z-score")),
  read_record_sheets(files[["weathering"]], "weathering", c("total data", "total IQR", "total z-score"))
)

# Dataset S6 contains one pooled pCO2 operational record in its sole sheet.
s6 <- suppressMessages(read_excel(files[["pCO2"]], sheet = "all data", col_names = TRUE))
pco2_iqr_cols <- which(grepl("^pCO2 variability", names(s6)))
pco2_iqr_col <- max(pco2_iqr_cols)
pco2_records <- tibble(
  category = "pCO2",
  record_id = "pCO2::all data",
  Age_Ma = suppressWarnings(as.numeric(s6[[pco2_iqr_col + 1L]])) / 1000,
  IQR = suppressWarnings(as.numeric(s6[[pco2_iqr_col]]))
) |>
  filter(is.finite(Age_Ma), is.finite(IQR))
records <- bind_rows(records, pco2_records)

if (anyDuplicated(records[c("record_id", "Age_Ma")])) {
  stop("Duplicate operational-record/bin IQR rows detected in formal workbooks")
}

# The global bootstrap uses all operational records from S2-S6. The S5 v4
# deletion is therefore respected directly, while the observed Global EV line
# remains the already-calculated series supplied by the latest S1 workbook.
record_sources <- list(
  global = records,
  temperature = filter(records, category == "temperature"),
  d18O = filter(records, category == "d18O"),
  d13C = filter(records, category == "d13C"),
  weathering = filter(records, category == "weathering")
)

bootstrap_panel <- function(panel_id, observed_bins, record_data) {
  observed_fit <- npreg::ss(
    x = observed_bins$Age_Ma,
    y = observed_bins$EV,
    spar = spar_value
  )
  observed_pred <- predict(observed_fit, x = prediction_grid, se.fit = FALSE)

  cluster_ids <- sort(unique(record_data$record_id))
  n_clusters <- length(cluster_ids)
  record_data <- record_data |>
    mutate(cluster_index = match(record_id, cluster_ids))

  boot_pred <- matrix(
    NA_real_, nrow = length(prediction_grid), ncol = bootstrap_iterations
  )
  set.seed(random_seed)

  for (iteration in seq_len(bootstrap_iterations)) {
    sampled <- sample.int(n_clusters, size = n_clusters, replace = TRUE)
    weights <- tabulate(sampled, nbins = n_clusters)
    boot_bins <- record_data |>
      mutate(cluster_weight = weights[cluster_index]) |>
      filter(cluster_weight > 0) |>
      group_by(Age_Ma) |>
      summarise(EV = weighted.mean(IQR, w = cluster_weight), .groups = "drop") |>
      arrange(Age_Ma)

    fit <- npreg::ss(x = boot_bins$Age_Ma, y = boot_bins$EV, spar = spar_value)
    pred <- predict(fit, x = prediction_grid, se.fit = FALSE)
    boot_pred[, iteration] <- as.numeric(pred$y)

    if (iteration %% 500L == 0L) {
      message(panel_id, ": completed ", iteration, " / ", bootstrap_iterations)
    }
  }

  list(
    envelope = tibble(
      panel = panel_id,
      Age_Ma = prediction_grid,
      observed_npreg_fit = as.numeric(observed_pred$y),
      bootstrap_p2_5 = apply(boot_pred, 1, quantile, probs = 0.025, na.rm = TRUE, names = FALSE),
      bootstrap_p97_5 = apply(boot_pred, 1, quantile, probs = 0.975, na.rm = TRUE, names = FALSE),
      spar = spar_value
    ),
    metadata = tibble(
      panel = panel_id,
      observed_bin_count = nrow(observed_bins),
      operational_record_clusters = n_clusters,
      record_bin_IQR_rows = nrow(record_data),
      bootstrap_iterations = bootstrap_iterations,
      random_seed = random_seed,
      spar = spar_value
    )
  )
}

panel_order <- c("global", "temperature", "d18O", "d13C", "weathering")
envelope_cache <- file.path(output_dir, "Fig_S6_npreg_cluster_bootstrap_envelopes.csv")
metadata_cache <- file.path(output_dir, "Fig_S6_bootstrap_metadata.csv")
if (file.exists(envelope_cache) && file.exists(metadata_cache)) {
  cached_envelopes <- read_csv(envelope_cache, show_col_types = FALSE)
  cached_metadata <- read_csv(metadata_cache, show_col_types = FALSE)
  cache_is_complete <- all(panel_order %in% cached_envelopes$panel) &&
    all(panel_order %in% cached_metadata$panel) &&
    all(cached_metadata$bootstrap_iterations == bootstrap_iterations) &&
    all(cached_metadata$random_seed == random_seed) &&
    all(abs(cached_metadata$spar - spar_value) < 1e-12)
} else {
  cache_is_complete <- FALSE
}

if (cache_is_complete) {
  message("Reusing the completed 5,000-iteration bootstrap envelope cache")
  results <- lapply(panel_order, function(id) {
    list(
      envelope = filter(cached_envelopes, panel == id),
      metadata = filter(cached_metadata, panel == id)
    )
  })
  names(results) <- panel_order
} else {
  results <- lapply(panel_order, function(id) {
    bootstrap_panel(id, observed[[id]], record_sources[[id]])
  })
  names(results) <- panel_order
}

envelopes <- bind_rows(lapply(results, `[[`, "envelope"))
metadata <- bind_rows(lapply(results, `[[`, "metadata"))
observed_out <- bind_rows(lapply(panel_order, function(id) {
  mutate(observed[[id]], panel = id, .before = 1)
}))
record_out <- records |>
  arrange(category, record_id, Age_Ma)

write_csv(observed_out, file.path(output_dir, "Fig_S6_observed_50kyr_EV_from_formal_datasets.csv"))
write_csv(envelopes, file.path(output_dir, "Fig_S6_npreg_cluster_bootstrap_envelopes.csv"))
write_csv(metadata, file.path(output_dir, "Fig_S6_bootstrap_metadata.csv"))
write_csv(record_out, file.path(output_dir, "Fig_S6_record_level_IQR_bootstrap_source.csv"))

panel_specs <- list(
  global = list(letter = "A", label1 = "Global EV", label2 = NULL, side = 2, col = "black", fill = "#D9D9D9", lty = 1, ylim = c(0.10, 0.72), ticks = seq(0.1, 0.7, by = 0.1), axis_cex = 0.64),
  temperature = list(letter = "B", label1 = "Temperature", label2 = "variability", side = 4, col = "#D7191C", fill = "#FBB4AE", lty = 2, ylim = c(-1.10, 1.00), ticks = seq(-1.1, 1.0, by = 0.1), axis_cex = 0.36),
  d18O = list(letter = "C", label1 = expression(delta^{18}*O), label2 = "variability", side = 2, col = "#F46D43", fill = "#FDD9C9", lty = 2, ylim = c(0.00, 1.00), ticks = seq(0.0, 1.0, by = 0.1), axis_cex = 0.52),
  d13C = list(letter = "D", label1 = expression(delta^{13}*C), label2 = "variability", side = 4, col = "#F89432", fill = "#FDE2C5", lty = 2, ylim = c(0.00, 0.60), ticks = seq(0.0, 0.6, by = 0.1), axis_cex = 0.64),
  weathering = list(letter = "E", label1 = "Weathering", label2 = "variability", side = 2, col = "#FDCB3A", fill = "#FFF3C9", lty = 2, ylim = c(-0.20, 1.70), ticks = seq(-0.2, 1.7, by = 0.1), axis_cex = 0.38)
)

draw_events <- function(ylim) {
  rect(66.90, ylim[1], 66.75, ylim[2], col = adjustcolor("#B9DCF2", alpha.f = 0.38), border = NA)
  rect(66.40, ylim[1], 65.55, ylim[2], col = adjustcolor("#F6C6CB", alpha.f = 0.30), border = NA)
  rect(66.02, ylim[1], 65.98, ylim[2], col = adjustcolor("#C98E60", alpha.f = 0.25), border = NA)
  abline(v = seq(65, 67, by = 0.5), col = "grey84", lty = 3, lwd = 0.55)
  abline(v = 66.00, col = "#8B4513", lty = 3, lwd = 1.15)
}

draw_fig_s6 <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  layout(matrix(seq_along(panel_order), ncol = 1), heights = c(1.45, 1, 1, 1, 1.08))
  par(
    oma = c(3.3, 5.0, 4.1, 5.0),
    family = "sans", mgp = c(2.0, 0.48, 0), tcl = -0.22,
    lend = "round", las = 1, cex = 0.78
  )

  for (i in seq_along(panel_order)) {
    id <- panel_order[[i]]
    spec <- panel_specs[[id]]
    env <- results[[id]]$envelope
    par(mar = c(if (i == length(panel_order)) 0.7 else 0.05, 0.25, if (i == 1) 1.4 else 0.05, 0.25))
    plot(
      NA_real_, NA_real_, xlim = c(67.1, 64.7), ylim = spec$ylim,
      xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "", bty = "n"
    )
    draw_events(spec$ylim)
    polygon(
      c(env$Age_Ma, rev(env$Age_Ma)),
      c(env$bootstrap_p2_5, rev(env$bootstrap_p97_5)),
      col = adjustcolor(spec$fill, alpha.f = 0.88), border = NA
    )
    lines(env$Age_Ma, env$observed_npreg_fit, col = spec$col, lwd = 2.0, lty = spec$lty)

    axis(spec$side, at = spec$ticks, labels = format(spec$ticks, nsmall = 1), cex.axis = spec$axis_cex, lwd = 0, lwd.ticks = 0.75)
    if (is.null(spec$label2)) {
      mtext(spec$label1, side = spec$side, line = 2.25, las = 0, font = 2, cex = 0.84)
    } else {
      mtext(spec$label1, side = spec$side, line = 1.95, las = 0, font = 2, cex = 0.80)
      mtext(spec$label2, side = spec$side, line = 3.00, las = 0, font = 2, cex = 0.80)
    }
    mtext(spec$letter, side = if (spec$side == 2) 2 else 4, line = 4.15, las = 1, font = 2, cex = 0.96)

    segments(67.1, spec$ylim[1], 67.1, spec$ylim[2], lwd = 0.8)
    segments(64.7, spec$ylim[1], 64.7, spec$ylim[2], lwd = 0.8)
    if (i == 1) {
      segments(67.1, spec$ylim[2], 64.7, spec$ylim[2], lwd = 0.8)
      axis(3, at = seq(65, 67, by = 0.5), labels = format(seq(65, 67, by = 0.5), trim = TRUE), cex.axis = 0.82, line = 1.45, padj = -0.25)
      axis(3, at = seq(64.8, 67.1, by = 0.1), labels = FALSE, tcl = -0.14, lwd = 0, lwd.ticks = 0.65)
      usr <- par("usr")
      yspan <- diff(usr[3:4])
      ygeo <- usr[4] + 0.035 * yspan
      segments(67.1, ygeo, 64.7, ygeo, xpd = NA, lwd = 0.8)
      segments(66.0, ygeo, 66.0, usr[4], xpd = NA, lwd = 0.8)
      text(66.55, ygeo + 0.020 * yspan, "Cretaceous", xpd = NA, cex = 0.78)
      text(65.50, ygeo + 0.020 * yspan, "Paleogene", xpd = NA, cex = 0.78)

      text(67.03, 0.675, "Antarctic\nglaciation", col = "#4898D3", font = 2, cex = 0.72, adj = c(0, 0.5))
      segments(66.96, 0.655, 66.84, 0.61, col = "#4898D3", lwd = 0.8)
      text(65.53, 0.685, "Deccan eruptions", col = "#C92C2C", font = 2, cex = 0.72, adj = c(0, 0.5))
      segments(65.55, 0.665, 65.83, 0.62, col = "#C92C2C", lwd = 0.8)
      text(65.86, 0.645, "Chicxulub impact", col = "#7A3E14", font = 2, cex = 0.69, adj = c(0, 0.5))
      segments(65.88, 0.625, 66.00, 0.60, col = "#7A3E14", lwd = 0.8)

      text(66.98, 0.54, "Phase 1", font = 4, cex = 0.76, adj = c(0, 0.5))
      arrows(66.87, 0.51, 66.66, 0.57, length = 0.08, lwd = 2.0)
      text(66.46, 0.64, "Phase 2", font = 4, cex = 0.76, adj = c(0, 0.5))
      arrows(66.36, 0.625, 66.16, 0.705, length = 0.08, lwd = 2.0)
    }
    if (i == length(panel_order)) {
      segments(67.1, spec$ylim[1], 64.7, spec$ylim[1], lwd = 0.8)
      axis(1, at = seq(65, 67, by = 0.5), labels = format(seq(65, 67, by = 0.5), trim = TRUE), cex.axis = 0.82)
      axis(1, at = seq(64.8, 67.1, by = 0.1), labels = FALSE, tcl = -0.14, lwd = 0, lwd.ticks = 0.65)
      mtext("Age (Ma)", side = 1, line = 2.15, font = 2, cex = 0.94)
    }
  }
}

base_name <- file.path(output_dir, "Fig_S6_revised_npreg_spar02_cluster_bootstrap")

png(paste0(base_name, ".png"), width = 1800, height = 3000, res = 400, bg = "white")
draw_fig_s6()
dev.off()

tiff(paste0(base_name, ".tiff"), width = 1800, height = 3000, res = 400, compression = "lzw", bg = "white")
draw_fig_s6()
dev.off()

pdf(paste0(base_name, ".pdf"), width = 4.5, height = 7.5, useDingbats = FALSE)
draw_fig_s6()
dev.off()

svglite::svglite(paste0(base_name, ".svg"), width = 4.5, height = 7.5, bg = "white")
draw_fig_s6()
dev.off()

message("Revised Fig. S6 created in: ", output_dir)
message("All five panels used npreg::ss(spar = 0.2); 5000 record-cluster bootstrap iterations; seed = 20260820")
