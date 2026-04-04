# Example: Run Pearsen analysis using Global EV and biotic turnover rates

library(ggplot2)
library(dplyr)
library(zoo)
library(broom)
library(purrr)

# =============================================================================
# 1. CLIMATE DATA: Global EV
# Source: Table_S1__K-Pg_global_EV.xlsx, 'total IQR' sheet
#
# Two resolutions are used, matched to each taxon's rate bin width:
#   Elasmobranchii (~50kyr bins): 'Global EV.3' / 'Midpoint (ka).3'  → n = 49
#   Bivalvia & Nannofossils (~100kyr bins): 'Global EV.2' / 'Midpoint (ka).2' → n = 25
#
# Both series have full coverage (no missing bins), no interpolation needed
# within the EV time axis itself.
# =============================================================================

# 50kyr bin EV — for Elasmobranchii
time_ev_50 <- c(64.80, 64.85, 64.90, 64.95, 65.00, 65.05, 65.10, 65.15,
                65.20, 65.25, 65.30, 65.35, 65.40, 65.45, 65.50, 65.55,
                65.60, 65.65, 65.70, 65.75, 65.80, 65.85, 65.90, 65.95,
                66.00, 66.05, 66.10, 66.15, 66.20, 66.25, 66.30, 66.35,
                66.40, 66.45, 66.50, 66.55, 66.60, 66.65, 66.70, 66.75,
                66.80, 66.85, 66.90, 66.95, 67.00, 67.05, 67.10, 67.15,
                67.20)

ev_50 <- c(0.236800, 0.253634, 0.246739, 0.269819, 0.211499,
           0.164941, 0.449895, 0.339298, 0.230257, 0.312379,
           0.299465, 0.452992, 0.285010, 0.297225, 0.267897,
           0.317576, 0.259260, 0.360594, 0.296934, 0.258625,
           0.460565, 0.426274, 0.298131, 0.465425, 0.632332,
           0.696408, 0.513247, 0.583953, 0.592345, 0.564910,
           0.415765, 0.515287, 0.417777, 0.281256, 0.503640,
           0.536776, 0.595411, 0.455077, 0.301400, 0.448441,
           0.372686, 0.328310, 0.231403, 0.308152, 0.357173,
           0.244886, 0.191892, 0.341171, 0.388731)

# 100kyr bin EV — for Bivalvia and Calcareous Nannofossils
time_ev_100 <- c(64.80, 64.90, 65.00, 65.10, 65.20, 65.30, 65.40, 65.50,
                 65.60, 65.70, 65.80, 65.90, 66.00, 66.10, 66.20, 66.30,
                 66.40, 66.50, 66.60, 66.70, 66.80, 66.90, 67.00, 67.10,
                 67.20)

ev_100 <- c(0.278299, 0.303629, 0.260387, 0.383262, 0.307016,
            0.419718, 0.420686, 0.414123, 0.316894, 0.418716,
            0.513620, 0.435601, 0.754867, 0.617985, 0.704774,
            0.688383, 0.441404, 0.458406, 0.561418, 0.353562,
            0.404050, 0.261582, 0.356962, 0.291659, 0.448724)

# =============================================================================
# 2. BIOLOGICAL RATE DATA
# =============================================================================

# --- Elasmobranchii (n = 49, ~50kyr bins) ---
# Source: E_speciation_rate.csv, E_extinction_rate.csv, E_net_rate.csv,
#         E__rate_time.csv
df_elas <- data.frame(
  time_rate = c(64.810683, 64.860692, 64.910701, 64.960710, 65.010719,
                65.060728, 65.110737, 65.160746, 65.210755, 65.260764,
                65.310773, 65.360782, 65.410791, 65.460800, 65.510809,
                65.560818, 65.610827, 65.660836, 65.710845, 65.760854,
                65.810863, 65.860872, 65.910881, 65.960890, 66.010899,
                66.060908, 66.110917, 66.160927, 66.210936, 66.260945,
                66.310954, 66.360963, 66.410972, 66.460981, 66.510990,
                66.560999, 66.611008, 66.661017, 66.711026, 66.761035,
                66.811044, 66.861053, 66.911062, 66.961071, 67.011080,
                67.061089, 67.111098, 67.161107, 67.211116),
  speciation = c(0.213093, 0.216426, 0.217815, 0.218840, 0.219003,
                 0.220514, 0.227480, 0.246492, 0.437353, 0.512906,
                 0.611321, 0.646538, 0.899599, 1.000040, 1.000040,
                 1.001713, 1.008846, 1.011243, 1.015021, 1.015021,
                 1.014961, 0.916589, 0.841276, 0.841276, 0.841276,
                 0.797777, 0.728498, 0.589965, 0.154168, 0.138281,
                 0.137711, 0.120275, 0.111150, 0.110718, 0.110718,
                 0.110718, 0.107570, 0.107570, 0.107428, 0.107428,
                 0.107428, 0.107431, 0.107445, 0.107445, 0.107465,
                 0.107465, 0.107526, 0.107547, 0.107609),
  extinction = c(0.148022, 0.148022, 0.147975, 0.147939, 0.147939,
                 0.147939, 0.147954, 0.147954, 0.147954, 0.147931,
                 0.147931, 0.147979, 0.147979, 0.152997, 0.159141,
                 0.159141, 0.159141, 0.165859, 0.178716, 1.118476,
                 1.802119, 1.817298, 1.946004, 1.946260, 1.946260,
                 1.948998, 1.949111, 1.949111, 1.949111, 1.949111,
                 1.675583, 1.596995, 1.596907, 1.401776, 1.401776,
                 1.339231, 0.177863, 0.177863, 0.170245, 0.170245,
                 0.170245, 0.170245, 0.170163, 0.169041, 0.169028,
                 0.169028, 0.168995, 0.168992, 0.168980),
  net        = c( 0.119090,  0.154497,  0.166038,  0.171018,  0.172200,
                  0.195583,  0.253426,  0.380954,  0.430306,  0.448767,
                  0.504815,  0.523346,  0.799101,  0.800750,  0.583981,
                  0.586050,  0.593146,  0.367756,  0.217063,  0.032841,
                  -0.622899, -0.773813, -1.142271, -1.142545, -1.142545,
                  -1.226364, -1.301260, -1.419641, -1.519404, -1.617488,
                  -1.190022, -1.182899, -1.270239, -1.047327, -1.047327,
                  -1.001813, -0.430932, -0.430937, -0.130647, -0.130647,
                  -0.130647, -0.130167, -0.121367, -0.068709, -0.068518,
                  -0.068501, -0.068028, -0.067947, -0.067021)
)

# --- Bivalvia (n = 25, ~100kyr bins) ---
# Source: B_speciation_rate.csv, B_extinction_rate.csv, B_net_rate.csv,
#         B__rate_time.csv
df_biva <- data.frame(
  time_rate = c(64.792569, 64.892703, 64.992838, 65.092972, 65.193106,
                65.293241, 65.393375, 65.493509, 65.593643, 65.693778,
                65.793912, 65.894046, 65.994181, 66.094315, 66.194449,
                66.294584, 66.394718, 66.494852, 66.594986, 66.695121,
                66.795255, 66.895389, 66.995524, 67.095658, 67.195792),
  speciation = c(0.233803, 0.238854, 0.242317, 0.247810, 0.255345,
                 0.317666, 0.439343, 0.602097, 0.684834, 0.691013,
                 0.691013, 0.690178, 0.588693, 0.588693, 0.100234,
                 0.098629, 0.095247, 0.094770, 0.094718, 0.094710,
                 0.094710, 0.094704, 0.094706, 0.094706, 0.094714),
  extinction = c(0.219499, 0.219362, 0.219333, 0.219333, 0.219325,
                 0.220008, 0.223387, 0.223387, 0.238097, 2.038715,
                 2.039182, 2.290239, 2.290239, 2.290239, 1.382503,
                 0.546581, 0.397762, 0.292792, 0.292792, 0.159844,
                 0.157188, 0.155514, 0.151730, 0.150994, 0.150582),
  net        = c( 0.038189,  0.076577,  0.092347,  0.118138,  0.144167,
                  0.204099,  0.126307,  0.311004,  0.011188, -0.945948,
                  -1.067259, -1.516614, -1.676550, -1.676104, -1.148518,
                  -0.776735, -0.614289, -0.419128, -0.420140, -0.162140,
                  -0.146922, -0.135405, -0.088226, -0.081187, -0.077408)
)

# --- Calcareous Nannofossils (n = 25, 100kyr bins) ---
# Source: N_s/e/n___100kyr_64_8-67_2_predict_.csv, N__rate_time__predict_.csv
df_nano <- data.frame(
  time_rate = c(64.8, 64.9, 65.0, 65.1, 65.2, 65.3, 65.4, 65.5, 65.6,
                65.7, 65.8, 65.9, 66.0, 66.1, 66.2, 66.3, 66.4, 66.5,
                66.6, 66.7, 66.8, 66.9, 67.0, 67.1, 67.2),
  speciation = c(33.877959, 34.240131, 34.412947, 34.405595, 34.227263,
                 33.887136, 33.394403, 32.758250, 31.987864, 31.092434,
                 30.081145, 28.963186, 27.747742, 26.444002, 25.061153,
                 23.608381, 22.094874, 20.529819, 18.922802, 17.285009,
                 15.628023, 13.963427, 12.302807, 10.657746,  9.039828),
  extinction = c( 7.528582, 10.596456, 13.987250, 17.633330, 21.467061,
                  25.420810, 29.426939, 33.417816, 37.325805, 41.083272,
                  44.622581, 47.876097, 50.776187, 53.255215, 55.245547,
                  56.679547, 57.489581, 57.608014, 56.989894, 55.680995,
                  53.749778, 51.264699, 48.294216, 44.906788, 41.170873),
  net        = c( 26.349376,  23.643675,  20.425697,  16.772265,  12.760202,
                  8.466326,   3.967464,  -0.659566,  -5.337941,  -9.990838,
                  -14.541436, -18.912911, -23.028445, -26.811213, -30.184394,
                  -33.071166, -35.394707, -37.078195, -38.067092, -38.395986,
                  -38.121755, -37.301272, -35.991409, -34.249042, -32.131045)
)

# =============================================================================
# 3. INTERPOLATION: align Global EV to each taxon's rate time axis
#    Elasmobranchii → 50kyr EV (time_ev_50 / ev_50)
#    Bivalvia & Nannofossils → 100kyr EV (time_ev_100 / ev_100)
#    EV is interpolated onto each taxon's true (non-uniform) rate time axis.
# =============================================================================

interpolate_ev <- function(time_target, time_ev, ev) {
  zoo_ev    <- zoo(ev, order.by = time_ev)
  approx_ev <- na.approx(zoo_ev, xout = time_target, rule = 2)
  return(coredata(approx_ev))
}

df_elas$ev_aligned <- interpolate_ev(df_elas$time_rate, time_ev_50,  ev_50)
df_biva$ev_aligned <- interpolate_ev(df_biva$time_rate, time_ev_100, ev_100)
df_nano$ev_aligned <- interpolate_ev(df_nano$time_rate, time_ev_100, ev_100)

df_elas$group <- "Elasmobranchii"
df_biva$group <- "Bivalvia"
df_nano$group <- "Calcareous Nannofossils"

# =============================================================================
# 4. RESHAPE TO LONG FORMAT
# =============================================================================

rate_types  <- c("speciation", "extinction", "net")
rate_labels <- c("Speciation Rate", "Extinction Rate", "Net Diversification Rate")
taxon_order <- c("Elasmobranchii", "Bivalvia", "Calcareous Nannofossils")

df_all <- bind_rows(df_elas, df_biva, df_nano) %>%
  tidyr::pivot_longer(cols      = all_of(rate_types),
                      names_to  = "rate_type",
                      values_to = "rate") %>%
  mutate(rate_type = factor(rate_type, levels = rate_types,
                            labels = rate_labels),
         group     = factor(group, levels = taxon_order))

# Standardize rate within each taxon x rate_type group
df_all <- df_all %>%
  group_by(group, rate_type) %>%
  mutate(rate_standardized = as.numeric(scale(rate))) %>%
  ungroup()

# =============================================================================
# 5. LINEAR FIT + PEARSON CORRELATION FOR EACH TAXON x RATE_TYPE CELL
# =============================================================================

fit_cell <- function(df_sub) {
  fit  <- lm(rate_standardized ~ ev_aligned, data = df_sub)
  pred <- predict(fit, interval = "confidence")
  df_sub <- cbind(df_sub, pred)
  df_sub$fit_values <- fit$fitted.values
  ct <- cor.test(df_sub$rate_standardized, df_sub$ev_aligned, method = "pearson")
  df_sub$r_label <- sprintf("r = %.3f\np = %s",
                            ct$estimate,
                            format.pval(ct$p.value, digits = 2))
  df_sub
}

df_plot <- df_all %>%
  group_by(group, rate_type) %>%
  group_modify(~ fit_cell(.x)) %>%
  ungroup()

# =============================================================================
# 6. PRINT CORRELATION TABLE
# =============================================================================

cat("\n=== Pearson Correlation: Global EV (50kyr) vs Biological Rates ===\n\n")
cat(sprintf("%-28s %-26s %s\n", "Taxon", "Rate Type", "r (p-value)"))
cat(strrep("-", 72), "\n")

df_all %>%
  group_by(group, rate_type) %>%
  summarise(
    r = cor(rate_standardized, ev_aligned, method = "pearson"),
    p = cor.test(rate_standardized, ev_aligned, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  arrange(factor(group, levels = taxon_order), rate_type) %>%
  rowwise() %>%
  mutate(line = sprintf("%-28s %-26s r = %7.4f  (p = %s)",
                        as.character(group),
                        as.character(rate_type),
                        r,
                        format.pval(p, digits = 3))) %>%
  pull(line) %>%
  cat(sep = "\n")

# =============================================================================
# 7. PLOT: 3 (rate types) x 3 (taxa) facet grid
# =============================================================================

label_df <- df_plot %>%
  group_by(group, rate_type) %>%
  slice(1) %>%
  ungroup() %>%
  select(group, rate_type, r_label)

taxon_colors <- c("Elasmobranchii"         = "#E63946",
                  "Bivalvia"                = "#457B9D",
                  "Calcareous Nannofossils" = "#2A9D8F")

p <- ggplot(df_plot,
            aes(x = ev_aligned, y = rate_standardized,
                color = group, fill = group)) +
  geom_point(alpha = 0.60, size = 1.8) +
  geom_line(aes(y = fit_values), linewidth = 0.85) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_text(data = label_df,
            aes(label = r_label),
            x = Inf, y = Inf,
            hjust = 1.08, vjust = 1.3,
            size = 3.0, color = "black",
            inherit.aes = FALSE) +
  facet_grid(rate_type ~ group, scales = "free_y") +
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values  = taxon_colors) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "grey60"),
    panel.spacing    = unit(0.8, "lines")
  ) +
  labs(
    x        = "Global EV (50kyr bin, interpolated to rate time axis)",
    y        = "Standardized Rate",
    title    = "Pearson Correlation: Global EV vs Biological Rates",
    subtitle = "64.8-67.2 Ma  |  rows: rate type  |  columns: taxon"
  )

print(p)

ggsave("pearson_global_ev_all_rates.pdf",
       plot   = p,
       width  = 11,
       height = 9,
       device = "pdf")

cat("\nFigure saved: pearson_global_ev_all_rates.pdf\n")