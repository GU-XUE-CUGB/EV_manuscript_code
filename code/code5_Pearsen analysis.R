# Example: Run Pearsen analysis using Global EV and biotic turnover rates

library(ggplot2)
library(dplyr)
library(zoo)
library(broom)
library(purrr)

# =============================================================================
# 1. CLIMATE DATA AND BIOLOGICAL RATE DATA
# =============================================================================
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1])))
file_path <- file.path(here, "modules/Code08_EV_proxy_biota_Pearson/inputs/current_variability_biotic_pairs.csv")
environment_path <- file.path(here, "modules/Code08_EV_proxy_biota_Pearson/inputs/current_environmental_variability_inputs.csv")
data <- subset(read.csv(file_path), EV_version=="Final revised" & predictor=="Global EV")
environment <- subset(read.csv(environment_path), predictor=="Global EV")
output_path <- file.path(here, "results/generated/code5")
dir.create(output_path, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# 2. BIOLOGICAL RATE DATA
# =============================================================================
prepare_data <- function(taxon, width) {
  d <- data[data$taxon==taxon, ]
  stopifnot(all(d$EV_bin_width_kyr==width))
  e <- environment[environment$EV_bin_width_kyr==width & is.finite(environment$variability), ]
  e <- e[order(e$age_Ma), ]
  eligible <- function(age) {
    if (any(abs(e$age_Ma-age)<1e-8)) return(TRUE)
    lo <- e$age_Ma[e$age_Ma<age]; hi <- e$age_Ma[e$age_Ma>age]
    length(lo)>0 && length(hi)>0 && min(hi)-max(lo)<=width/1000+1e-8
  }
  d <- d[vapply(d$rate_age_Ma, eligible, logical(1)) &
           is.finite(d$raw_rate) & is.finite(d$aligned_variability), ]
  out <- d[d$rate_type=="Speciation", c("rate_age_Ma", "aligned_variability", "raw_rate")]
  names(out) <- c("time_rate", "ev_aligned", "speciation")
  for (variable in c("Extinction", "Net diversification")) {
    a <- d[d$rate_type==variable, ]
    idx <- match(out$time_rate, a$rate_age_Ma)
    stopifnot(!anyNA(idx), !anyDuplicated(a$rate_age_Ma))
    out[[if(variable=="Extinction")"extinction" else "net"]] <- a$raw_rate[idx]
  }
  out[order(out$time_rate, decreasing=TRUE), ]
}
df_elas <- prepare_data("Elasmobranchs", 50)
df_biva <- prepare_data("Bivalves", 100)
df_nano <- prepare_data("Calcareous nannofossils", 100)
stopifnot(nrow(df_elas)==48, nrow(df_biva)==24, nrow(df_nano)==25)

# =============================================================================
# 3. AGE ALIGNMENT
# =============================================================================
# Reuse the final aligned observations; no new interpolation or endpoint filling.

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

cat("\n=== Pearson Correlation: Matched-resolution Global EV vs Biological Rates ===\n\n")
cat(sprintf("%-28s %-26s %s\n", "Taxon", "Rate Type", "r (p-value)"))
cat(strrep("-", 72), "\n")

df_all %>%
  group_by(group, rate_type) %>%
  summarise(
    r = cor(rate_standardized, ev_aligned, method = "pearson"),
    p = cor.test(rate_standardized, ev_aligned, method = "pearson")$p.value,
    n = sum(is.finite(rate_standardized) & is.finite(ev_aligned)),
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
    x        = "Global EV (matched 50- or 100-kyr resolution)",
    y        = "Standardized Rate",
    title    = "Pearson Correlation: Global EV vs Biological Rates",
    subtitle = "64.8-67.2 Ma  |  rows: rate type  |  columns: taxon"
  )

if (interactive() || "--plot" %in% commandArgs(TRUE)) {
print(p)

ggsave(file.path(output_path, "pearson_global_ev_all_rates.pdf"),
       plot   = p,
       width  = 11,
       height = 9,
       device = "pdf")

cat("\nFigure saved: pearson_global_ev_all_rates.pdf\n")
}
# Full proxy-specific and mean-state Pearson outputs
for (script in c("Code08_EV_proxy_biota_Pearson/Supplementary_Code08.R",
                 "Code09_mean_state_biota_Pearson/Supplementary_Code09.R",
                 "Code08_EV_proxy_biota_Pearson/Supplementary_pCO2_variability.R")) {
  extra <- if(grepl("Code08",script,fixed=TRUE))"--current-inputs" else character()
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(here, "modules", script)),extra))
  stopifnot(status == 0)
}
