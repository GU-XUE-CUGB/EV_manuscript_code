#!/usr/bin/env Rscript
# Diagram only: no input reconstruction, statistical analysis or resampling.
suppressPackageStartupMessages({library(grid); library(svglite); library(ragg); library(systemfonts)})
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
out <- file.path(root, "04_figures", "Fig_S1_concise_workflow_Arial_20261009")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
stem <- file.path(out, "Fig_S1_EV_workflow_concise_Arial")
font <- "Arial"
font_paths <- systemfonts::match_fonts(c(font, font), weight = c("normal", "bold"))$path
stopifnot(all(grepl("Arial", basename(font_paths))), capabilities("aqua"))
W <- 640; H <- 1008
palette <- list(green = c("#E2F4EE", "#21836C", "#105747"),
                neutral = c("#F1EFE7", "#76756D", "#30302C"),
                purple = c("#EEEDFD", "#7163D2", "#483D94"),
                orange = c("#FBEDE7", "#B75430", "#843B1C"),
                gold = c("#FBF0DB", "#A67127", "#70450E"))
txt <- function(label, x, y, size = 10, bold = FALSE, color = "#333333", just = "centre", rot = 0) {
  Encoding(label) <- "UTF-8"
  # Compose super/subscripts from Arial digits, without fallback fonts.
  if (grepl("[¹³⁸₂₃₁]", label)) {
    rows <- strsplit(label, "\n", fixed = TRUE)[[1]]
    digit <- c("¹" = "1", "³" = "3", "⁸" = "8", "₂" = "2", "₃" = "3", "₁" = "1")
    for (j in seq_along(rows)) {
      chars <- strsplit(rows[j], "", fixed = TRUE)[[1]]
      tokens <- character(); kinds <- logical(); buffer <- ""
      for (ch in chars) {
        if (ch %in% names(digit)) {
          if (nzchar(buffer)) {tokens <- c(tokens, buffer); kinds <- c(kinds, FALSE); buffer <- ""}
          tokens <- c(tokens, ch); kinds <- c(kinds, TRUE)
        } else buffer <- paste0(buffer, ch)
      }
      if (nzchar(buffer)) {tokens <- c(tokens, buffer); kinds <- c(kinds, FALSE)}
      sizes <- ifelse(kinds, size * 0.72, size)
      labels <- tokens
      labels[kinds] <- unname(digit[tokens[kinds]])
      widths <- vapply(seq_along(labels), function(i) {
        systemfonts::string_width(labels[i], family = font, bold = bold,
                                  size = sizes[i], res = 720, include_bearing = FALSE) / 10
      }, numeric(1))
      cursor <- x - if (just == "centre") sum(widths) / 2 else 0
      for (i in seq_along(labels)) {
        offset <- if (!kinds[i]) 0 else if (tokens[i] %in% c("¹", "³", "⁸")) size * 0.34 else -size * 0.18
        grid.text(labels[i], unit(cursor, "native"),
                  unit(y + ((length(rows) + 1) / 2 - j) * size * 1.2 + offset, "native"),
                  just = "left", gp = gpar(fontfamily = font, fontsize = sizes[i],
                  fontface = if (bold) "bold" else "plain", col = color))
        cursor <- cursor + widths[i]
      }
    }
    return(invisible(NULL))
  }
  grid.text(label, unit(x, "native"), unit(y, "native"), just = just, rot = rot,
            gp = gpar(fontfamily = font, fontsize = size,
                      fontface = if (bold) "bold" else "plain", col = color))
}
box <- function(x, y, w, h, title, body = NULL, pal = "neutral", title_size = 11, body_size = 9.2) {
  p <- palette[[pal]]
  grid.roundrect(x = unit(x, "native"), y = unit(y, "native"),
                 width = unit(w, "native"), height = unit(h, "native"),
                 r = unit(2, "pt"), gp = gpar(fill = p[1], col = p[2], lwd = 0.9))
  if (is.null(body)) txt(title, x, y, title_size, TRUE, p[3]) else {
    lines <- strsplit(body, "\n", fixed = TRUE)[[1]]
    gap <- body_size * 1.30
    total <- title_size + 7 + length(lines) * gap
    start <- y + total / 2 - title_size / 2
    txt(title, x, start, title_size, TRUE, p[3])
    for (i in seq_along(lines))
      txt(lines[i], x, start - title_size / 2 - 7 - (i - 0.5) * gap, body_size, color = p[3])
  }
}
connector <- function(x, y, arrow_end = FALSE) {
  grid.lines(unit(x, "native"), unit(y, "native"),
             arrow = if (arrow_end) arrow(type = "closed", length = unit(4.5, "pt")) else NULL,
             gp = gpar(col = "#999999", lwd = 1.0))
}
stage <- function(n, top, bottom) {
  connector(c(19, 19), c(top, bottom))
  txt(paste("STAGE", n), 7, (top + bottom) / 2, 9.5, TRUE, "#999999", rot = 90)
}

draw <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, W), yscale = c(0, H)))
  stage(1, 999, 791)
  stage(2, 752, 587)
  stage(3, 552, 365)
  stage(4, 327, 31)

  txt("Data compilation", 338, 993, 14.5, TRUE, "#242424")
  proxy_x <- seq(87, 579, length.out = 5)
  proxy_labels <- c("Temperature\nrecords", "δ¹³C\nrecords", "δ¹⁸O\nrecords",
                    "Weathering\nisotope data", "Atmospheric\npCO₂ data")
  for (i in seq_along(proxy_x)) {
    box(proxy_x[i], 945, 112, 55, proxy_labels[i], pal = "green", title_size = 11)
    connector(rep(proxy_x[i], 2), c(917.5, 900))
  }
  connector(range(proxy_x), c(900, 900))
  connector(c(338, 338), c(900, 882), TRUE)
  box(338, 853, 550, 57, "Quality filtering and operational records",
      "Source, units, material/method, technical replicates and age-model checks\nOperational record: proxy × site/section × material/method",
      body_size = 9.0)
  connector(c(338, 338), c(824.5, 814), TRUE)
  box(338, 789, 550, 47, "Proxy redundancy removal",
      "Remove five δ¹⁸O-derived temperature records and derived Wang et al. pCO₂\nRetain independent records and newly compiled measurements",
      body_size = 9.0)

  connector(c(338, 338), c(765.5, 752), TRUE)
  txt("Standardization & binning", 338, 735, 14.5, TRUE, "#242424")
  box(181, 689, 270, 61, "Record-wide z-standardization",
      "z = (measurement − record mean) / sample SD\nUse the retained full record; SD has n − 1 denominator",
      pal = "purple", title_size = 10.5, body_size = 8.5)
  box(487, 689, 270, 61, "Time binning",
      "50 / 100 / 150 / 200 kyr, with fixed bin definitions\nMain reconstruction: 49 bins, 67.2–64.8 Ma",
      pal = "purple", title_size = 10.5, body_size = 8.7)
  connector(c(181, 181), c(658.5, 642), TRUE)
  connector(c(487, 487), c(658.5, 642), TRUE)
  box(334, 614, 576, 56, "Within-record, within-bin variability: IQR",
      "IQR = Q₃ − Q₁ of standardized measurements (linear quartiles)\nContributing record-bin: n ≥ 2 measurements and IQR > 0",
      pal = "purple", body_size = 9.3)

  connector(c(334, 334), c(586, 560), TRUE)
  txt("Environmental variability (EV) reconstruction", 338, 544, 14.1, TRUE, "#242424")
  xs <- c(134, 338, 542)
  connector(c(134, 542), c(522, 522))
  for (x in xs) connector(c(x, x), c(522, 513), TRUE)
  box(134, 477, 182, 71, "Proxy-specific EV",
      "Mean contributing record IQR\nwithin each proxy category",
      pal = "orange", title_size = 10.7, body_size = 9.0)
  box(338, 477, 182, 71, "Latitude-specific EV",
      "Mean record IQR within low\n(≤30°) / mid-high (>30°) latitudes",
      pal = "orange", title_size = 10.7, body_size = 8.6)
  box(542, 477, 182, 71, "Global EV",
      "Arithmetic mean of all\ncontributing record-level IQRs\nNot a mean of proxy means",
      pal = "orange", title_size = 10.7, body_size = 8.6)
  for (x in xs) connector(c(x, x), c(441.5, 430))
  connector(c(134, 542), c(430, 430))
  connector(c(338, 338), c(430, 419), TRUE)
  box(338, 386, 550, 65, "Bin-level npreg smoothing and bootstrap envelope",
      "Fit bin-level EV only: spar = 0.2 (50/100 kyr); 0.3 (150/200 kyr)\nLatitude curves: identical spar = 0.25 for both groups\nGlobal / proxy: 95% record-cluster bootstrap envelope (5,000 draws)",
      body_size = 8.9)

  connector(c(338, 338), c(353.5, 332), TRUE)
  txt("Robustness testing & uncertainty sensitivity", 338, 317, 14.0, TRUE, "#242424")
  # Parallel overview; execution details remain in Methods and code.
  box(181, 276, 270, 58, "Temporal resolution (Fig. S3)",
      "50 / 100 / 150 / 200 kyr bins\n5,000 record-cluster bootstrap iterations",
      pal = "gold", title_size = 10.4, body_size = 8.9)
  box(487, 276, 270, 58, "Sampling associations (Fig. S4)",
      "EV vs record, measurement and site counts\nSpearman: raw series and first differences",
      pal = "gold", title_size = 10.4, body_size = 8.9)
  box(181, 209, 270, 58, "Fixed-record rarefaction (Fig. S5A)",
      "3 / 5 / 7 records per eligible bin\n10,000 rarefaction iterations",
      pal = "gold", title_size = 10.4, body_size = 8.9)
  box(487, 209, 270, 58, "Joint rarefaction (Fig. S5B)",
      "5 records/bin; 3 / 4 / 5 measurements/record-bin\n10,000 rarefaction iterations",
      pal = "gold", title_size = 10.4, body_size = 8.9)
  box(181, 142, 270, 58, "Measurement uncertainty (Fig. S6)",
      "Documented and conservative method-based scenarios\n10,000 simulations; chronology fixed",
      pal = "gold", title_size = 10.4, body_size = 8.9)
  box(487, 142, 270, 58, "Chronological uncertainty (Fig. S7)",
      "Documented and conservative method-based scenarios\n10,000 simulations; measurement values fixed",
      pal = "gold", title_size = 10.4, body_size = 8.9)
  box(181, 75, 270, 58, "Joint uncertainty (Fig. S8)",
      "Conservative measurement + chronology\n10,000 simulations; joint propagation",
      pal = "gold", title_size = 10.4, body_size = 8.9)
  box(487, 75, 270, 58, "Latitude sensitivity (Figs. S14–S15)",
      "Site-balanced rarefaction and leave-one-site-out\n10,000 rarefaction iterations",
      pal = "gold", title_size = 10.0, body_size = 8.9)
  txt("95% rarefaction / simulation intervals: 2.5th–97.5th percentiles",
      338, 28, 8.5, color = "#555555")
  txt("Main sampling and uncertainty sensitivity analyses use 50-kyr EV.",
      338, 13, 8.5, color = "#555555")
  popViewport()
}

quartz(type = "pdf", file = paste0(stem, ".pdf"), width = W / 72, height = H / 72, family = font)
draw(); dev.off()
svglite(paste0(stem, ".svg"), width = W / 72, height = H / 72,
        system_fonts = list(sans = font), fix_text_size = FALSE)
draw(); dev.off()
agg_png(paste0(stem, ".png"), width = W / 72, height = H / 72, units = "in", res = 400)
draw(); dev.off()
agg_tiff(paste0(stem, ".tiff"), width = W / 72, height = H / 72, units = "in", res = 400, compression = "lzw")
draw(); dev.off()
message("Workflow diagram only; Arial verified; original files unchanged.")
message(stem)

