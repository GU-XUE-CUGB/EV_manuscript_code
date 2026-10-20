#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(npreg)
  library(patchwork)
  library(jsonlite)
  library(svglite)
  library(png)
  library(grid)
})

# 1. Confirmed inputs
script_arg <- sub("^--file=", "", gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE)[grepl(
  "^--file=", gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE))])
work <- dirname(normalizePath(script_arg))
latitude_file <- file.path(work, "output", "descriptive_EV.csv")
rarefaction_file <- file.path(work, "output", "site_balanced_rarefaction_results.json")
fig5_source <- file.path(work, "assets", "Fig5_ABC_confirmed.png")
out <- file.path(work, "output", "figures")
stopifnot(all(file.exists(c(latitude_file, rarefaction_file, fig5_source))))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

lat <- read.csv(latitude_file)[, c("Age_Ma", "low_EV", "low_records", "mid_high_EV", "mid_high_records")]
names(lat) <- c("Age_Ma", "Low_EV", "Low_records", "Mid_high_EV", "Mid_high_records")
lat[] <- lapply(lat, as.numeric)
stopifnot(nrow(lat) == 49L, all(diff(lat$Age_Ma) < 0),
          sum(is.finite(lat$Low_EV)) == 48L, sum(is.finite(lat$Mid_high_EV)) == 49L)
rare <- jsonlite::fromJSON(rarefaction_file)
stopifnot(rare$metadata$iterations == 10000L,
          sum(rare$bin_summary$eligible) == 20L,
          rare$phase_summary$eligible_bins[rare$phase_summary$phase == "Pre"] == 0L,
          rare$phase_summary$eligible_bins[rare$phase_summary$phase == "Phase 1"] == 1L,
          rare$phase_summary$eligible_bins[rare$phase_summary$phase == "Phase 2"] == 6L)
phase2 <- rare$phase_summary[rare$phase_summary$phase == "Phase 2", ]
stopifnot(abs(phase2$median - 0.025217029776706267) < 1e-12,
          abs(phase2$p2_5 + 0.13841051298398155) < 1e-12,
          abs(phase2$p97_5 - 0.20283406864998446) < 1e-12)
counts <- rare$bin_summary[, c("midpoint_Ma", "low_sites", "mid_high_sites")]
stopifnot(max(abs(counts$midpoint_Ma - lat$Age_Ma)) < 1e-9)
boundary <- 66.016
low_color <- "#FFC247"
high_color <- "#A67ECE"
boundary_color <- "#8C5632"

# 2. Figure S14: current EV estimates and observed geographic coverage
ev_long <- rbind(
  data.frame(Age_Ma = lat$Age_Ma, EV = lat$Low_EV, Group = "Low latitude"),
  data.frame(Age_Ma = lat$Age_Ma, EV = lat$Mid_high_EV, Group = "Mid-high latitude")
)
ev_long$Group <- factor(ev_long$Group, levels = c("Low latitude", "Mid-high latitude"))
ev_long <- ev_long[order(ev_long$Group, -ev_long$Age_Ma), ]
ev_long$Run <- ave(is.finite(ev_long$EV), ev_long$Group,
                   FUN = function(z) cumsum(!z))
fit_low <- npreg::ss(lat$Age_Ma[is.finite(lat$Low_EV)],
                     lat$Low_EV[is.finite(lat$Low_EV)], spar = 0.25)
fit_high <- npreg::ss(lat$Age_Ma, lat$Mid_high_EV, spar = 0.25)
age_grid <- seq(64.8, 67.2, by = 0.005)
fits <- rbind(
  data.frame(Age_Ma = age_grid,
             EV = as.numeric(predict(fit_low, x = age_grid)$y), Group = "Low latitude"),
  data.frame(Age_Ma = age_grid,
             EV = as.numeric(predict(fit_high, x = age_grid)$y), Group = "Mid-high latitude")
)
fits$Group <- factor(fits$Group, levels = levels(ev_long$Group))
coverage <- rbind(
  data.frame(Age_Ma = counts$midpoint_Ma, Sites = counts$low_sites, Group = "Low latitude"),
  data.frame(Age_Ma = counts$midpoint_Ma, Sites = counts$mid_high_sites,
             Group = "Mid-high latitude")
)
coverage$Group <- factor(coverage$Group, levels = levels(ev_long$Group))
palette <- c("Low latitude" = low_color, "Mid-high latitude" = high_color)
geologic_x <- scale_x_reverse(limits = c(67.2, 64.8),
                              breaks = c(67, 66.5, 66, 65.5, 65),
                              minor_breaks = seq(64.8, 67.2, by = 0.1),
                              expand = c(0, 0))
original_theme <- theme_bw(base_family = "sans", base_size = 9) +
  theme(panel.grid.major = element_line(colour = "#EDEDED", linewidth = 0.22),
        panel.grid.minor = element_line(colour = "#F5F5F5", linewidth = 0.13),
        panel.border = element_rect(colour = "#1F1F1F", linewidth = 0.30),
        axis.title = element_text(face = "bold", size = 9),
        axis.text = element_text(colour = "#222222", size = 8),
        legend.title = element_blank(),
        legend.position = "top", legend.text = element_text(size = 8),
        plot.margin = margin(3, 5, 3, 3, unit = "pt"))

s14_top <- ggplot() +
  geom_vline(xintercept = boundary, colour = boundary_color,
             linetype = "longdash", linewidth = 0.55) +
  geom_line(data = ev_long[is.finite(ev_long$EV), ],
            aes(Age_Ma, EV, colour = Group, group = interaction(Group, Run)),
            linewidth = 0.33, alpha = 0.68, show.legend = FALSE) +
  geom_point(data = ev_long[is.finite(ev_long$EV), ],
             aes(Age_Ma, EV, colour = Group), size = 0.75, alpha = 0.80,
             show.legend = FALSE) +
  geom_line(data = fits, aes(Age_Ma, EV, colour = Group), linewidth = 1.12) +
  scale_colour_manual(values = palette) + geologic_x +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, by = 0.2),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "EV index") + original_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
s14_bottom <- ggplot(coverage, aes(Age_Ma, Sites, colour = Group)) +
  geom_vline(xintercept = boundary, colour = boundary_color,
             linetype = "longdash", linewidth = 0.55) +
  geom_line(linewidth = 0.56) + geom_point(size = 0.88) +
  scale_colour_manual(values = palette) + geologic_x +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, by = 2),
                     expand = c(0, 0)) +
  labs(x = "Age (Ma)", y = "Independent sites") + original_theme +
  theme(legend.position = "none")
s14 <- s14_top / s14_bottom + plot_layout(heights = c(3.4, 1))

# 3. Figure S15: phase summaries and leave-one-site-out sensitivity
bins <- rare$bin_summary[rare$bin_summary$eligible, ]
stopifnot(nrow(bins) == 20L,
          all(is.finite(bins$median)),
          all(bins$p2_5 <= bins$median & bins$median <= bins$p97_5))
s15_a <- ggplot(bins, aes(midpoint_Ma, median)) +
  geom_hline(yintercept = 0, colour = "#777777", linetype = "dashed",
             linewidth = 0.40) +
  geom_vline(xintercept = boundary, colour = boundary_color,
             linetype = "longdash", linewidth = 0.55) +
  geom_errorbar(aes(ymin = p2_5, ymax = p97_5), colour = "#7E7E7E",
                width = 0.019, linewidth = 0.38) +
  geom_point(colour = high_color, size = 1.65, shape = 16) + geologic_x +
  scale_y_continuous(limits = c(-0.5, 1.25), breaks = seq(-0.5, 1.0, by = 0.25),
                     expand = c(0, 0)) +
  labs(x = "Age (Ma)", y = expression(Delta*"EV (mid-high - low)")) +
  annotate("text", x = 67.17, y = 1.13, label = "A", fontface = "bold",
           hjust = 0, size = 3.6) + original_theme +
  theme(legend.position = "none")
phases <- rare$phase_summary
phases$phase <- factor(phases$phase, levels = c("Pre", "Phase 1", "Phase 2"))
s15_b <- ggplot(phases, aes(phase, median)) +
  geom_hline(yintercept = 0, colour = "#777777", linetype = "dashed",
             linewidth = 0.40) +
  geom_errorbar(data = phases[is.finite(phases$median), ],
                aes(ymin = p2_5, ymax = p97_5), colour = "#7E7E7E",
                width = 0.12, linewidth = 0.55) +
  geom_point(data = phases[is.finite(phases$median), ],
             colour = high_color, size = 2.15) +
  annotate("text", x = 1, y = 0.40, label = "0 eligible bins\nNot estimable", size = 2.7) +
  annotate("text", x = 2, y = -0.32, label = "1 bin only\nNo phase-level inference", size = 2.3) +
  annotate("text", x = 3, y = 0.32, label = "6 eligible bins", size = 2.8) +
  scale_y_continuous(limits = c(-0.5, 1.25), breaks = seq(-0.5, 1.0, by = 0.25),
                     expand = c(0, 0)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = expression(Delta*"EV (mid-high - low)")) +
  original_theme + theme(legend.position = "none")
s15_phase <- s15_b +
  labs(title = "A  Phase-level contrasts") +
  theme(plot.title = element_text(face = "bold", size = 10))
phase2_ages <- bins$midpoint_Ma[bins$midpoint_Ma >= 66.05 & bins$midpoint_Ma <= 66.40]
phase2_sites <- unique(rare$site_contributions$canonical_site[
  rare$site_contributions$age_Ma %in% phase2_ages])
loo <- rare$loo_summary[rare$loo_summary$phase == "Phase 2" &
                         rare$loo_summary$removed_site %in% phase2_sites, ]
loo <- loo[order(loo$removed_site), ]
stopifnot(nrow(loo) == 25L)
loo$Label <- paste0(sub(",.*", "", loo$removed_site), " (", loo$eligible_bins_after_removal, ")")
loo$Label <- factor(loo$Label, levels = rev(loo$Label))
loo$Mark <- ifelse(loo$removed_site == "ODP 1209", "ODP 1209", "Other site")
s15_loo <- ggplot(loo, aes(median, Label)) +
  geom_vline(xintercept = 0, colour = "#777777", linetype = "dashed", linewidth = 0.4) +
  geom_segment(aes(x = p2_5, xend = p97_5, yend = Label), colour = "#7E7E7E", linewidth = 0.5) +
  geom_point(aes(shape = Mark), colour = high_color, size = 1.7) +
  scale_shape_manual(values = c("ODP 1209" = 18, "Other site" = 16)) +
  labs(x = expression(Delta*"EV (mid-high - low)"), y = NULL,
       title = "B  Leave-one-site-out: Phase 2",
       subtitle = "Sites removed (eligible bins remaining)") + original_theme +
  theme(legend.position = "none", axis.text.y = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7), panel.grid.minor = element_blank())
s15 <- s15_phase / s15_loo + plot_layout(heights = c(1, 3.0))

save_gg <- function(plot, name, width, height) {
  stem <- file.path(out, name)
  ggsave(paste0(stem, ".pdf"), plot = plot, device = grDevices::pdf,
         width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(paste0(stem, ".svg"), plot = plot, device = svglite::svglite,
         width = width, height = height, units = "in")
  ggsave(paste0(stem, ".png"), plot = plot, device = ragg::agg_png,
         width = width, height = height, units = "in", dpi = 400)
}
save_gg(s14, "Fig_S14_latitude_EV_final", 4.54, 4.15)
save_gg(s15, "Fig_S15_latitude_robustness_final", 4.54, 7.2)

# 4. Figure 5: preserve the confirmed A/B/C image; append a vector D panel
original_art <- png::readPNG(fig5_source)
stopifnot(identical(dim(original_art)[1:2], c(2832L, 2116L)))
W <- 1058; H0 <- 1416; HD <- 410; H <- H0 + HD
x_age <- function(z) 611.75 + (66 - z) * 335.75
y_delta <- function(z) 1740 - (z + 0.5) / 1.75 * 240
line_native <- function(x, y, color, lwd = 1.5, lty = "solid")
  grid.lines(x = unit(x, "native"), y = unit(y, "native"),
             gp = gpar(col = color, lwd = lwd, lty = lty, lineend = "round"))
text_native <- function(label, x, y, fontsize, fontface = "plain", just = "centre")
  grid.text(label, x = unit(x, "native"), y = unit(y, "native"),
            just = just, gp = gpar(fontfamily = "sans", fontsize = fontsize,
                                  fontface = fontface, col = "#1A1A1A"))
draw_fig5 <- function() {
  grid.newpage()
  grid.raster(original_art, x = 0.5,
              y = unit((HD + H0 / 2) / H, "npc"),
              width = unit(1, "npc"), height = unit(H0 / H, "npc"),
              interpolate = FALSE)
  pushViewport(viewport(xscale = c(0, W), yscale = c(H, 0)))
  text_native("D", 12, 1444, 21, "bold", "left")
  text_native("Site-balanced latitudinal EV contrast", 164, 1451, 14, "bold", "left")
  text_native("Median and 95% rarefaction interval", 1028, 1480, 9, "plain", "right")
  for (window in list(c(66.925, 66.425), c(66.425, 66.025))) {
    grid.rect(x = unit(mean(x_age(window)), "native"), y = unit(1620, "native"),
              width = unit(abs(diff(x_age(window))), "native"), height = unit(240, "native"),
              gp = gpar(fill = "#FDECEF", col = NA))
  }
  text_native("Phase 1", x_age(66.675), 1511, 9)
  text_native("Phase 2", x_age(66.225), 1511, 9)
  line_native(c(x_age(67.2), x_age(64.8)), rep(y_delta(0), 2), "#888888", 1, "dashed")
  line_native(rep(x_age(boundary), 2), c(1500,1740), boundary_color, 2, "longdash")
  for (i in seq_len(nrow(bins))) {
    xx <- x_age(bins$midpoint_Ma[i])
    line_native(rep(xx, 2), y_delta(c(bins$p2_5[i], bins$p97_5[i])), "#7E7E7E", 1.1)
    for (z in c(bins$p2_5[i], bins$p97_5[i]))
      line_native(c(xx-3,xx+3), rep(y_delta(z),2), "#7E7E7E", 1.1)
    grid.points(x=unit(xx,"native"), y=unit(y_delta(bins$median[i]),"native"),
                pch=16, size=unit(1.8,"mm"), gp=gpar(col=high_color))
  }
  line_native(c(x_age(67.2),x_age(64.8)), c(1740,1740), "#1A1A1A", 1.4)
  line_native(rep(x_age(67.2),2), c(1500,1740), "#1A1A1A", 1.4)
  for (z in seq(64.8,67.2,by=.1))
    line_native(rep(x_age(z),2), c(1740,1745), "#1A1A1A", .8)
  for (z in c(67,66.5,66,65.5,65)) {
    line_native(rep(x_age(z),2), c(1740,1749), "#1A1A1A", 1.1)
    text_native(format(z,trim=TRUE), x_age(z), 1766, 13)
  }
  for (z in c(-.5,0,.5,1)) {
    line_native(c(x_age(67.2)-6,x_age(67.2)), rep(y_delta(z),2), "#1A1A1A", 1)
    text_native(sprintf("%.1f",z), x_age(67.2)-12,y_delta(z),11,just="right")
  }
  grid.text(expression(Delta*"EV (mid-high - low)"), x=unit(130,"native"),
            y=unit(1620,"native"),rot=90,gp=gpar(fontsize=12,fontface="bold",fontfamily="sans"))
  text_native("Age (Ma)", 611.75, 1806, 16, "bold")
  popViewport()
}
fig5_stem <- file.path(out, "Fig_5_latitude_contrast_final")
grDevices::pdf(paste0(fig5_stem, ".pdf"), width = W / 144,
               height = H / 144, useDingbats = FALSE)
draw_fig5(); dev.off()
svglite::svglite(paste0(fig5_stem, ".svg"), width = W / 144,
                 height = H / 144, bg = "white")
draw_fig5(); dev.off()
ragg::agg_png(paste0(fig5_stem, ".png"), width = W * 2,
              height = H * 2, res = 288)
draw_fig5(); dev.off()

cat("Figures generated without running rarefaction: ", out, "\n", sep = "")
