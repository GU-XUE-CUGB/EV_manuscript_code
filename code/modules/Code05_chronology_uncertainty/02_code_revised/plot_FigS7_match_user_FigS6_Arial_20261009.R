#!/usr/bin/env Rscript
# Formatting only; reuse completed chronology summaries and saved npreg curves.
suppressPackageStartupMessages({library(grid); library(svglite); library(ragg); library(systemfonts)})
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
font_family <- "Arial"
font_paths <- systemfonts::match_fonts(c("Arial", "Arial"), weight = c("normal", "bold"))$path
stopifnot(all(grepl("Arial", basename(font_paths))), capabilities("aqua"))
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
input <- file.path(root, "reference_figures", "Chronological_sensitivity_Fig_S7")
output <- file.path(root, "04_figures", "Fig_S7_match_user_FigS6_Arial_20261009")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(output, "Fig_S7_revised_match_FigS6_Arial")

# Coordinates measured from the supplied 3072 x 3120 px, 300-dpi S6.
W <- 3072; H <- 3120; left <- 200; right <- 3052
boundary <- 66.016
cols <- c(observed = "#737373", median = "#2B6F9F",
          interval = "#AFC6D9", boundary = "#7B3F28")
read_panel <- function(prefix, coverage_name) {
  paths <- file.path(input, c(paste0(prefix, "_plot_data.csv"),
                             paste0(prefix, "_npreg_predictions.csv"), coverage_name))
  result <- lapply(paths, read.csv, check.names = FALSE)
  names(result) <- c("ev", "smooth", "coverage")
  stopifnot(nrow(result$ev) == 49, nrow(result$coverage) == 49,
            all(diff(result$ev$Age_Ma) > 0),
            max(abs(result$ev$Age_Ma - result$coverage$Age_Ma)) < 1e-10,
            all(result$ev$simulation_p2_5 <= result$ev$simulation_median),
            all(result$ev$simulation_median <= result$ev$simulation_p97_5),
            all(result$coverage$coverage >= 0 & result$coverage$coverage <= 1),
            all(is.finite(as.matrix(result$ev))),
            max(result$ev$simulation_p97_5) <= 0.90)
  result$paths <- paths
  result
}
A <- read_panel("Fig_S7_chronological_sensitivity",
                "Fig_S7A_chronology_perturbation_coverage.csv")
B <- read_panel("Fig_S7B_conservative_method_based_chronological_sensitivity",
                "Fig_S7B_chronology_perturbation_coverage.csv")
stopifnot(max(abs(A$ev$observed_EV - B$ev$observed_EV)) < 1e-12)
inputs <- c(A$paths, B$paths)
input_hashes <- tools::md5sum(inputs)

px <- function(x) unit(x / W, "npc")
py <- function(y) unit(1 - y / H, "npc")
age_x <- function(age) left + (67.2 - age) / 2.4 * (right - left)
line_px <- function(x, y, colour = "black", width = 0.8, lty = 1) {
  grid.lines(px(x), py(y), gp = gpar(col = colour, lwd = width, lty = lty))
}
text_px <- function(label, x, y, size = 10, bold = FALSE,
                    colour = "black", just = "left", rot = 0) {
  Encoding(label) <- "UTF-8"
  grid.text(label, x = px(x), y = py(y), just = just, rot = rot,
            gp = gpar(fontfamily = font_family, fontsize = size,
                      fontface = if (bold) "bold" else "plain", col = colour))
}
panel <- function(d, offset, letter, subtitle, boundary_label = FALSE) {
  top <- 127 + offset; bottom <- 1127 + offset
  covtop <- 1172 + offset; covbottom <- 1353 + offset
  ev_y <- function(value) bottom - (value - 0.15) / 0.75 * (bottom - top)
  cov_y <- function(value) covbottom - value * (covbottom - covtop)
  text_px(letter, 0, 26 + offset, 14.3, TRUE)
  text_px(subtitle, 86, 26 + offset, 14.3, TRUE)

  # Only the completed interval, observed trajectory, median and saved smooth.
  pushViewport(viewport(x = px(left), y = py(bottom),
                        width = px(right - left), height = unit((bottom - top) / H, "npc"),
                        just = c("left", "bottom"), clip = "on"))
  nx <- function(age) (67.2 - age) / 2.4
  ny <- function(value) (value - 0.15) / 0.75
  e <- d$ev
  grid.polygon(x = c(nx(e$Age_Ma), rev(nx(e$Age_Ma))),
               y = c(ny(e$simulation_p2_5), rev(ny(e$simulation_p97_5))),
               gp = gpar(fill = cols["interval"], col = NA, alpha = 0.48))
  grid.lines(nx(e$Age_Ma), ny(e$observed_EV),
             gp = gpar(col = cols["observed"], lwd = 0.85))
  grid.lines(nx(e$Age_Ma), ny(e$simulation_median),
             gp = gpar(col = cols["median"], lwd = 1.20))
  grid.points(nx(e$Age_Ma), ny(e$observed_EV), pch = 16,
              size = unit(1.9, "mm"), gp = gpar(col = cols["observed"]))
  grid.lines(nx(d$smooth$Age_Ma), ny(d$smooth$npreg_fitted),
             gp = gpar(col = "black", lwd = 2.0))
  popViewport()

  # Match the reference's plain axes, 0.1 EV ticks and bold titles.
  line_px(c(left, left), c(top, bottom), width = 0.85)
  line_px(c(left, right), rep(1146 + offset, 2), width = 0.85)
  for (v in seq(0.2, 0.9, 0.1)) {
    line_px(c(left - 16, left), rep(ev_y(v), 2), width = 0.85)
    text_px(sprintf("%.1f", v), left - 27, ev_y(v), 10, just = "right")
  }
  text_px("Global EV", 35, (top + bottom) / 2, 12, TRUE, just = "centre", rot = 90)

  cov <- d$coverage
  grid.polygon(px(c(age_x(cov$Age_Ma[1]), age_x(cov$Age_Ma),
                    age_x(tail(cov$Age_Ma, 1)))),
               py(c(covbottom, cov_y(cov$coverage), covbottom)),
               gp = gpar(fill = cols["interval"], col = NA, alpha = 0.72))
  line_px(age_x(cov$Age_Ma), cov_y(cov$coverage), cols["median"], 1.2)
  line_px(c(left, left), c(covtop, covbottom), width = 0.85)
  line_px(c(left, right), rep(covbottom, 2), width = 0.85)
  for (v in c(0, 0.25, 0.5, 0.75, 1)) {
    line_px(c(left - 16, left), rep(cov_y(v), 2), width = 0.85)
    if (v %in% c(0, 0.5, 1))
      text_px(sprintf("%d%%", round(v * 100)), left - 27, cov_y(v), 10, just = "right")
  }
  text_px("Chronology", left + 20, covtop + 35, 12, TRUE)
  text_px("perturbation coverage (%)", left + 20, covtop + 88, 12, TRUE)

  for (age in round(seq(64.8, 67.2, 0.1), 1)) {
    major <- abs(age * 2 - round(age * 2)) < 1e-8
    line_px(rep(age_x(age), 2), c(covbottom, covbottom + if (major) 16 else 8), width = 0.8)
    if (major) text_px(sprintf("%.1f", age), age_x(age), covbottom + 41, 10, just = "centre")
  }
  text_px("Age (Ma)", (left + right) / 2, covbottom + 105, 12, TRUE, just = "centre")
  line_px(rep(age_x(boundary), 2), c(top, covbottom), cols["boundary"], 1.7, "55")
  if (boundary_label)
    text_px("K–Pg boundary", age_x(boundary), top - 35, 12, TRUE,
            cols["boundary"], "centre")
}

draw <- function() {
  grid.newpage()
  panel(A, 0, "A", "Documented chronological uncertainty", TRUE)
  panel(B, 1540, "B", "Conservative method-based chronological uncertainty")
  y <- 3095
  line_px(c(594, 654), c(y, y), cols["observed"], 0.85)
  grid.points(px(624), py(y), pch = 16, size = unit(1.9, "mm"), gp = gpar(col = cols["observed"]))
  text_px("Observed Global EV", 677, y, 10)
  line_px(c(1134, 1184), c(y, y), cols["median"], 1.2)
  text_px("Simulation median", 1207, y, 10)
  line_px(c(1630, 1680), c(y, y), "black", 2)
  text_px("npreg, spar = 0.2", 1703, y, 10)
  grid.rect(x = px(2103 + 34), y = py(y), width = px(68), height = unit(57 / H, "npc"),
            gp = gpar(fill = cols["interval"], col = NA, alpha = 0.48))
  text_px("95% simulation interval", 2195, y, 10)
}

# Editable vectors and raster exports; source summaries are read-only.
quartz(type = "pdf", file = paste0(stem, ".pdf"), width = W / 300, height = H / 300,
       family = font_family)
draw(); dev.off()
svglite(paste0(stem, ".svg"), width = W / 300, height = H / 300,
        system_fonts = list(sans = font_family), fix_text_size = FALSE)
draw(); dev.off()
agg_png(paste0(stem, ".png"), width = W, height = H, units = "px", res = 300)
draw(); dev.off()
agg_tiff(paste0(stem, ".tiff"), width = W, height = H, units = "px", res = 300,
         compression = "lzw")
draw(); dev.off()
stopifnot(identical(input_hashes, tools::md5sum(inputs)))
message("S7 formatting matched to supplied S6; six input files unchanged.")
message("Saved npreg predictions reused; no fitting or simulations performed.")
message("Arial font files: ", paste(font_paths, collapse = "; "))
message(stem)
