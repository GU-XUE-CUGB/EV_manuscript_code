#!/usr/bin/env Rscript
# Adopted window-equal background reference; no EV refitting/resampling.
suppressPackageStartupMessages({library(grid); library(svglite); library(ragg); library(systemfonts)})
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
input <- file.path(root, "plot_snapshot_inputs", "04_figures", "Global_EV_GAM")
out <- file.path(root, "outputs", "plots", "Fig_2A_background_comparison_Arial_20261009")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input, c("Global_EV_50kyr_raw_values.csv",
                           "Global_EV_npreg_spar02_predictions.csv",
                           "Global_EV_npreg_spar02_cluster_bootstrap_envelope.csv"))
hashes <- tools::md5sum(paths)
raw <- read.csv(paths[1]); fit <- read.csv(paths[2]); ci <- read.csv(paths[3])
stopifnot(nrow(raw) == 49, nrow(fit) == nrow(ci),
          max(abs(fit$Age_Ma - ci$Age_Ma)) < 1e-12)
pre <- raw$Global_EV[raw$Age_Ma >= 66.95 - 1e-8 & raw$Age_Ma <= 67.2 + 1e-8]
post <- raw$Global_EV[raw$Age_Ma >= 64.8 - 1e-8 & raw$Age_Ma <= 65.45 + 1e-8]
stopifnot(length(pre) == 6, length(post) == 14)
values <- c(Bin_equal_weight = mean(c(pre, post)),
            Window_equal_weight = mean(c(mean(pre), mean(post))))
stopifnot(round(values[1], 5) == 0.35944, round(values[2], 5) == 0.34503)
font <- "Arial"
stopifnot(all(grepl("Arial", basename(systemfonts::match_fonts(c(font, font),
          weight = c("normal", "bold"))$path))), capabilities("aqua"))
W <- 560; H <- 360
L <- 65; R <- 546; T <- 36; B <- 311
x <- function(age) L + (67.2 - age) / 2.4 * (R - L)
y <- function(ev) B - ev / 0.85 * (B - T)
xx <- function(v) unit(v / W, "npc")
yy <- function(v) unit(1 - v / H, "npc")
text <- function(label, a, b, size = 10, bold = FALSE, just = "centre", col = "black", rot = 0) {
  Encoding(label) <- "UTF-8"
  grid.text(label, xx(a), yy(b), just = just, rot = rot,
            gp = gpar(fontfamily = font, fontsize = size,
                      fontface = if (bold) "bold" else "plain", col = col))
}
line <- function(a, b, col = "black", lwd = 0.85, lty = 1) {
  grid.lines(xx(a), yy(b), gp = gpar(col = col, lwd = lwd, lty = lty))
}
draw <- function(reference) {
  grid.newpage()
  grid.polygon(xx(c(x(ci$Age_Ma), rev(x(ci$Age_Ma)))),
               yy(c(y(ci$bootstrap_p2_5), rev(y(ci$bootstrap_p97_5)))),
               gp = gpar(fill = "#D9D9D9", col = NA, alpha = 0.50))
  line(c(L, R), rep(y(reference), 2), "#999999", 0.95, "dotted")
  grid.points(xx(x(raw$Age_Ma)), yy(y(raw$Global_EV)), pch = 16,
              size = unit(1.8, "mm"), gp = gpar(col = "#8C8C8C"))
  line(x(fit$Age_Ma), y(fit$npreg_fitted_spar_0_2), lwd = 2.0)
  line(rep(x(66.016), 2), c(T, B), "#82491E", 1.25, "dashed")
  text("K–Pg boundary", x(66.016), 17, 10, TRUE, col = "#82491E")
  text(sprintf("Background = %.5f", reference), R - 7, T + 13,
       8.8, TRUE, just = "right", col = "#8C8C8C")
  line(c(L, R, R, L, L), c(T, T, B, B, T))
  for (v in seq(0, 0.8, 0.1)) {
    line(c(L - 5, L), rep(y(v), 2))
    text(sprintf("%.1f", v), L - 9, y(v), 9.5, just = "right")
  }
  for (age in round(seq(64.8, 67.2, 0.1), 1)) {
    major <- abs(age * 2 - round(age * 2)) < 1e-8
    line(rep(x(age), 2), c(B, B + if (major) 5 else 2.5))
    if (major) text(sprintf("%.1f", age), x(age), B + 17, 9.5)
  }
  text("A", 12, T, 15, TRUE)
  text("EV index", 19, (T + B) / 2, 12, TRUE, rot = 90)
  text("Age (Ma)", (L + R) / 2, H - 10, 12, TRUE)
}
selected <- if ("--all" %in% commandArgs(TRUE)) names(values) else "Window_equal_weight"
for (name in selected) {
  stem <- file.path(out, paste0("Fig_2A_background_", name, "_", sprintf("%.5f", values[name])))
  quartz(type = "pdf", file = paste0(stem, ".pdf"), width = W / 72, height = H / 72, family = font)
  draw(values[name]); dev.off()
  svglite(paste0(stem, ".svg"), width = W / 72, height = H / 72,
          system_fonts = list(sans = font), fix_text_size = FALSE)
  draw(values[name]); dev.off()
  agg_png(paste0(stem, ".png"), width = W / 72, height = H / 72, units = "in", res = 400)
  draw(values[name]); dev.off()
}
stopifnot(identical(hashes, tools::md5sum(paths)))
message("Adopted background: equal weights for the pre and post windows (0.34503).")
message("Standalone Fig. 2A reference panel; complete original Fig. 2 unchanged.")
message("Observed EV, npreg predictions and bootstrap interval reused without changes.")
