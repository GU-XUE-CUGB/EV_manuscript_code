#!/usr/bin/env Rscript
# Calculate npreg fits and record-cluster envelopes from record trajectories.
suppressPackageStartupMessages(library(npreg))
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- dirname(normalizePath(script))
cli <- commandArgs(TRUE)
validate_only <- "--validate-only" %in% cli
count_arg <- cli[grepl("^--iterations=", cli)]
iterations <- if (length(count_arg)) as.integer(sub("^--iterations=", "", count_arg[1])) else 5000L
stopifnot(is.finite(iterations), iterations >= 2)
label <- if (validate_only) "validation" else if (iterations == 5000L) "formal_5000" else paste0("smoke_", iterations)
out <- file.path(root, "outputs", label)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
iqrs <- read.csv(file.path(root, "inputs", "record_level_IQR_all_bin_widths_revised_v4_n2.csv"), stringsAsFactors=FALSE)
reference <- file.path(root, "reference_results")
ref_s3 <- read.csv(file.path(reference, "Fig_S3_final_npreg_bootstrap_envelopes.csv"))
ref_bins <- read.csv(file.path(reference, "Fig_S3_revised_observed_Global_EV.csv"))
ref_global <- read.csv(file.path(reference, "Global_EV_npreg_spar02_predictions.csv"))
ref_proxy <- read.csv(file.path(reference, "Fig_S6_npreg_cluster_bootstrap_envelopes.csv"))
RNGkind(kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")
random_seed <- 20260820L

fit_panel <- function(source, observed, grid, spar, id) {
  source <- source[order(source$record_id, source$Age_Ma),]
  observed <- observed[order(observed$Age_Ma),]
  stopifnot(!anyDuplicated(source[c("record_id", "Age_Ma")]))
  fit <- npreg::ss(x=observed$Age_Ma, y=observed$EV, spar=spar)
  prediction <- as.numeric(predict(fit, x=grid, se.fit=FALSE)$y)
  envelope <- data.frame(panel=id, Age_Ma=grid, observed_npreg_fit=prediction, spar=spar)
  clusters <- sort(unique(source$record_id))
  n_clusters <- length(clusters)
  source$record_index <- match(source$record_id, clusters)
  if (!validate_only) {
    boot <- matrix(NA_real_, nrow=length(grid), ncol=iterations)
    available <- integer(iterations)
    set.seed(random_seed) # reset independently for every panel
    for (i in seq_len(iterations)) {
      sampled <- sample.int(n_clusters, size=n_clusters, replace=TRUE)
      weights <- tabulate(sampled, nbins=n_clusters)
      source$cluster_weight <- weights[source$record_index]
      included <- source[source$cluster_weight > 0,]
      by_age <- split(included, included$Age_Ma)
      bins <- do.call(rbind, lapply(by_age, function(g) data.frame(
        Age_Ma=g$Age_Ma[1], EV=weighted.mean(g$IQR, w=g$cluster_weight))))
      bins <- bins[order(bins$Age_Ma),]
      available[i] <- nrow(bins)
      bootstrap_fit <- npreg::ss(x=bins$Age_Ma, y=bins$EV, spar=spar)
      boot[,i] <- as.numeric(predict(bootstrap_fit, x=grid, se.fit=FALSE)$y)
    }
    envelope$bootstrap_median <- apply(boot, 1, median, na.rm=TRUE)
    envelope$bootstrap_p2_5 <- apply(boot, 1, quantile, probs=.025, na.rm=TRUE, names=FALSE)
    envelope$bootstrap_p97_5 <- apply(boot, 1, quantile, probs=.975, na.rm=TRUE, names=FALSE)
  } else available <- integer()
  metadata <- data.frame(panel=id, observed_bin_count=nrow(observed),
    operational_record_clusters=n_clusters, record_bin_IQR_rows=nrow(source),
    iterations=if (validate_only) 0 else iterations, seed=random_seed, spar=spar,
    prediction_grid_points=length(grid),
    minimum_available_bins=if (length(available)) min(available) else NA_integer_)
  list(envelope=envelope, metadata=metadata, observed=observed)
}

checks <- list(); results <- list()
for (scheme in c("200", "150", "100", "50")) {
  spar <- if (scheme %in% c("200", "150")) .3 else .2
  source <- iqrs[as.character(iqrs$scheme) == scheme & iqrs$alignment == "standard" & iqrs$IQR != 0,
    c("record_id", "bin_center_Ma", "IQR")]
  names(source)[2] <- "Age_Ma"
  observed <- aggregate(IQR ~ Age_Ma, data=source, FUN=mean)
  names(observed)[2] <- "EV"
  locked <- ref_bins[ref_bins$bin_width_kyr == as.integer(scheme),]
  locked <- locked[order(locked$Age_Ma),]
  stopifnot(nrow(observed) == nrow(locked), max(abs(observed$Age_Ma - locked$Age_Ma)) < 1e-12)
  error <- max(abs(observed$EV - locked$Global_EV))
  stopifnot(error < 1e-12)
  grid <- seq(min(observed$Age_Ma), max(observed$Age_Ma), length.out=481L)
  result <- fit_panel(source, observed, grid, spar, paste0("S3_", scheme))
  existing <- ref_s3[ref_s3$bin_width_kyr == as.integer(scheme),]
  existing <- existing[order(existing$Age_Ma),]
  fit_error <- max(abs(result$envelope$observed_npreg_fit - existing$observed_npreg_fit))
  stopifnot(fit_error < 1e-9)
  checks[[scheme]] <- data.frame(panel=paste0("S3_", scheme), observed_mean_max_error=error, observed_fit_max_error=fit_error)
  if (!validate_only && iterations == 5000L) {
    for (column in c("bootstrap_median", "bootstrap_p2_5", "bootstrap_p97_5")) {
      stopifnot(max(abs(result$envelope[[column]] - existing[[column]])) < 1e-9)
    }
  }
  results[[paste0("S3_", scheme)]] <- result
  if (scheme == "50") {
    global_grid <- seq(min(observed$Age_Ma), max(observed$Age_Ma), length.out=481L)
    global_result <- fit_panel(source, observed, global_grid, .2, "Global_50kyr")
    stopifnot(max(abs(global_result$envelope$observed_npreg_fit - ref_global$npreg_fitted_spar_0_2)) < 1e-9)
    results[["Global_50kyr"]] <- global_result
    pre <- observed$EV[observed$Age_Ma >= 66.95-1e-8 & observed$Age_Ma <= 67.2+1e-8]
    post <- observed$EV[observed$Age_Ma >= 64.8-1e-8 & observed$Age_Ma <= 65.45+1e-8]
    background <- mean(c(mean(pre), mean(post)))
    stopifnot(length(pre)==6L, length(post)==14L, abs(background-0.34503257945945265)<1e-11)
    write.csv(data.frame(background_equal_window_mean=background, pre_mean=mean(pre), post_mean=mean(post)),
      file.path(out, "background_reference.csv"), row.names=FALSE)
  }
}

# Proxy panels use the supplied workbook-derived IQR and bin-mean tables,
# including the zero-IQR rows in those inputs.
proxy_source <- read.csv(file.path(root, "inputs", "Fig_S6_record_level_IQR_bootstrap_source.csv"), stringsAsFactors=FALSE)
proxy_bins <- read.csv(file.path(root, "inputs", "Fig_S6_observed_50kyr_EV_from_formal_datasets.csv"), stringsAsFactors=FALSE)
for (panel in c("global", "temperature", "d18O", "d13C", "weathering")) {
  source <- if (panel == "global") proxy_source else proxy_source[proxy_source$category == panel,]
  observed <- proxy_bins[proxy_bins$panel == panel, c("Age_Ma", "EV")]
  existing <- ref_proxy[ref_proxy$panel == panel,]
  existing <- existing[order(existing$Age_Ma),]
  grid <- seq(64.8, 67.2, length.out=481L)
  result <- fit_panel(source, observed, grid, .2, paste0("Ancillary_proxy_", panel))
  error <- max(abs(result$envelope$observed_npreg_fit - existing$observed_npreg_fit))
  stopifnot(error < 1e-9)
  if (!validate_only && iterations == 5000L) {
    for (column in c("bootstrap_p2_5", "bootstrap_p97_5")) {
      stopifnot(max(abs(result$envelope[[column]] - existing[[column]])) < 1e-9)
    }
  }
  checks[[paste0("proxy_", panel)]] <- data.frame(panel=paste0("Ancillary_proxy_", panel),
    observed_mean_max_error=NA_real_, observed_fit_max_error=error)
  results[[paste0("Ancillary_proxy_", panel)]] <- result
}

envelopes <- do.call(rbind, lapply(results, `[[`, "envelope"))
metadata <- do.call(rbind, lapply(results, `[[`, "metadata"))
write.csv(envelopes, file.path(out, "recomputed_npreg_predictions_and_envelopes.csv"), row.names=FALSE)
write.csv(metadata, file.path(out, "run_metadata.csv"), row.names=FALSE)
write.csv(do.call(rbind, checks), file.path(out, "deterministic_reference_checks.csv"), row.names=FALSE)
writeLines(c(paste0("output_status=", label), "historical_outputs_are_reference_only=true",
  "final_S3_spar: 50/100=.2, 150/200=.3", "bootstrap: complete operational record trajectories; 5000 formal iterations",
  "RNG: Mersenne-Twister/Inversion/Rejection; seed 20260820 reset each panel",
  "global and S3: exclude IQR=0; ancillary proxy execution: preserve exact workbook source including zeros",
  "formal historical 5000-iteration confidence envelopes were not rerun during package audit"),
  file.path(out, "run_status.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
print(metadata)
