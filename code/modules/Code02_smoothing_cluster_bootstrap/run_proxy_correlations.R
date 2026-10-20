#!/usr/bin/env Rscript
# Ancillary final Figure S11: 12 Pearson tests from the exact final bin inputs.
# Workbook extraction is the input handoff; O/P cell provenance is packaged.
argv <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
root <- dirname(normalizePath(script))
current_mode <- "--current-inputs" %in% commandArgs(TRUE)
input <- if(current_mode)file.path(root,"../Code08_EV_proxy_biota_Pearson/inputs/current_bin_level_EV_inputs.csv") else file.path(root,"inputs/final_bin_level_EV_inputs.csv")
wide <- read.csv(input)
reference <- read.csv(file.path(root, "reference_results", "EV_proxy_Pearson_12_tests_final.csv"))
phases <- list(c("Full interval",64800,67200), c("Phase 1",66450,66900), c("Phase 2",66050,66400))
proxies <- c("temperature", "d18O", "d13C", "weathering")
rows <- pairs <- exclusions <- list()
for (phase in phases) for (proxy in proxies) {
  lo <- as.numeric(phase[2]); hi <- as.numeric(phase[3])
  window <- wide[wide$midpoint_ka >= lo & wide$midpoint_ka <= hi,]
  valid <- window[complete.cases(window[c("global_EV", proxy)]),]
  valid$proxy_EV <- valid[[proxy]]
  test <- cor.test(valid$proxy_EV, valid$global_EV, method="pearson", alternative="two.sided")
  fit <- lm(global_EV ~ proxy_EV, data=valid)
  rows[[length(rows)+1L]] <- data.frame(phase=phase[1], proxy=proxy,
    young_bound_Ma=lo/1000, old_bound_Ma=hi/1000, n_pairs=nrow(valid),
    pearson_r=unname(test$estimate), p_value=test$p.value,
    slope=unname(coef(fit)[2]), intercept=unname(coef(fit)[1]),
    method="Pearson two-sided, unadjusted")
  pairs[[length(pairs)+1L]] <- data.frame(phase=phase[1], proxy=proxy,
    midpoint_ka=valid$midpoint_ka, age_Ma=valid$age_Ma,
    global_EV=valid$global_EV, proxy_EV=valid$proxy_EV)
  missing <- window$midpoint_ka[is.na(window[[proxy]])]
  if (length(missing)) exclusions[[length(exclusions)+1L]] <- data.frame(
    phase=phase[1], proxy=proxy, midpoint_ka=missing,
    reason="No stored proxy bin-level estimate; no interpolation")
}
results <- do.call(rbind, rows)
joined <- merge(results, reference, by=c("phase", "proxy"), suffixes=c("_new", "_ref"))
stopifnot(nrow(joined) == 12L)
if (!current_mode) stopifnot(all(joined$n_pairs_new == joined$n_pairs_ref))
errors <- vapply(c("pearson_r", "p_value", "slope", "intercept"), function(v)
  max(abs(joined[[paste0(v,"_new")]] - joined[[paste0(v,"_ref")]])), numeric(1))
if (!current_mode) stopifnot(max(errors) < 1e-12)
out <- file.path(root, "outputs", if(current_mode)"proxy_correlations_FigS11_current_dataset" else "proxy_correlations_FigS11")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write.csv(results, file.path(out, "EV_proxy_Pearson_12_tests_recomputed.csv"), row.names=FALSE)
write.csv(do.call(rbind, pairs), file.path(out, "paired_data_all_12_tests.csv"), row.names=FALSE)
write.csv(do.call(rbind, exclusions), file.path(out, "excluded_bins_by_test.csv"), row.names=FALSE)
write.csv(data.frame(quantity=names(errors), maximum_absolute_error=errors), file.path(out,"reference_checks.csv"), row.names=FALSE)
writeLines(c("output_status=newly_recomputed_deterministic", "figure_mapping=final Figure S11",
  "source_handoff=final_bin_level_EV_inputs.csv from read-only formal workbook O/P cells",
  "method=two-sided Pearson; unadjusted P; exact midpoint joins; no interpolation or smoothing",
  "plot_mapping=x proxy-specific EV; y Global EV"), file.path(out, "run_status.txt"))
print(results)
