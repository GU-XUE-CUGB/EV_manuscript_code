#!/usr/bin/env Rscript
# Draw Fig. S10 proxy panels.
# Preserve exact original drawing styles and completed envelope inputs.
suppressPackageStartupMessages({library(readr);library(dplyr);library(svglite)})
argv <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
root <- dirname(dirname(normalizePath(script)))
output_dir <- file.path(root,"outputs","plots","Ancillary_FigS10_proxy_bootstrap")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
envelopes <- read_csv(file.path(root,"reference_results","Fig_S6_npreg_cluster_bootstrap_envelopes.csv"),show_col_types=FALSE)
panel_order <- c("global","temperature","d18O","d13C","weathering")
results <- lapply(panel_order,function(id)list(envelope=filter(envelopes,panel==id)))
names(results) <- panel_order
source <- readLines(file.path(root,"provenance","Code_FigS6_revised_npreg_cluster_bootstrap.R"),warn=FALSE)
start <- grep("^panel_specs <- list",source)[1]
stopifnot(length(start)==1L,!is.na(start))
eval(parse(text=source[start:length(source)]),envir=.GlobalEnv)
writeLines("Plot-only former Fig S6 (final ancillary S10): locked historical confidence envelopes; no bootstrap rerun.",
 file.path(output_dir,"PLOT_INPUT_STATUS.txt"))
