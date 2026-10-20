# Load necessary library
library(npreg)
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1])))
cli <- commandArgs(TRUE)
input <- file.path(here, "modules/Code02_smoothing_cluster_bootstrap/plot_snapshot_inputs/04_figures/Global_EV_GAM/Global_EV_50kyr_raw_values.csv")
dataset <- read.csv(input)
dataset <- data.frame(Midpoint=dataset$Age_Ma, IQR=dataset$Global_EV)
stopifnot(nrow(dataset)==49, all(is.finite(as.matrix(dataset))))

# Apply the ss function to the dataset, specifying the 'Midpoint' and 'IQR' columns, and spar parameter
mod.ss <- with(dataset, ss(Midpoint, IQR, spar = 0.2))

# View the result
mod.ss

# Plot the result
if (interactive() || "--plot" %in% cli) plot(mod.ss)

# Full record-cluster bootstrap and temporal-resolution fits
if ("--cluster-bootstrap" %in% cli) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(file.path(here, "modules/Code02_smoothing_cluster_bootstrap/run.R")))
  stopifnot(status == 0)
}
if ("--proxy-correlations" %in% cli) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(here, "modules/Code02_smoothing_cluster_bootstrap/run_proxy_correlations.R")),"--current-inputs"))
  stopifnot(status == 0)
}
