#!/usr/bin/env Rscript
# Plot-only entry point. Exact final S3 drawing code is evaluated from the
# preserved source. Completed numerical inputs are explicit reference caches.
suppressPackageStartupMessages({library(readr);library(dplyr);library(svglite)})
argv <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
root <- dirname(dirname(normalizePath(script)))
output_dir <- file.path(root,"outputs","plots","Fig_S3_final_mixed_spar_cluster_bootstrap")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
source <- readLines(file.path(root,"provenance","Code2_FigS3_final_mixed_spar_cluster_bootstrap.R"),warn=FALSE)
scheme_start <- grep("^schemes <- list",source)[1]
scheme_end <- grep("^expected_bin_counts <-",source)[1]-1L
eval(parse(text=source[scheme_start:scheme_end]),envir=.GlobalEnv)
final_envelopes <- read_csv(file.path(root,"reference_results","Fig_S3_final_npreg_bootstrap_envelopes.csv"),show_col_types=FALSE)
observed_all <- read_csv(file.path(root,"reference_results","Fig_S3_revised_observed_Global_EV.csv"),show_col_types=FALSE)
final_envelopes_by_scheme <- observed_by_scheme <- list()
for (scheme in c("200","150","100","50")) {
  final_envelopes_by_scheme[[scheme]] <- filter(final_envelopes,bin_width_kyr==as.integer(scheme))
  observed_by_scheme[[scheme]] <- filter(observed_all,bin_width_kyr==as.integer(scheme))
}
plot_start <- grep("^envelope_min <-",source)[1]
stopifnot(length(plot_start)==1L,!is.na(plot_start))
eval(parse(text=source[plot_start:length(source)]),envir=.GlobalEnv)
writeLines("Plot-only export from historical locked numerical inputs; no 5000-iteration resampling was run.",
  file.path(output_dir,"PLOT_INPUT_STATUS.txt"))
