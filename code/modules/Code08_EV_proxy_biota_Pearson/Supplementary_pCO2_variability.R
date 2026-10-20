# Additional existing pCO2 column in the current Fig. 3 and Dataset S9.
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1])))
current_mode <- '--current-inputs' %in% commandArgs(TRUE)
paired_file <- if(current_mode)"inputs/current_variability_biotic_pairs.csv" else "inputs/variability_biotic_paired_data_final.csv"
environment_file <- if(current_mode)"inputs/current_environmental_variability_inputs.csv" else "inputs/final_environmental_variability_inputs.csv"
output <- file.path(here,if(current_mode)"results/generated/current_dataset" else "results/generated")
data <- subset(read.csv(file.path(here, paired_file)),
               EV_version=="Final revised" & predictor=="pCO2 variability")
environment <- subset(read.csv(file.path(here, environment_file)),
                      predictor=="pCO2 variability" & is.finite(variability))
results <- pairs <- list()
for (taxon in unique(data$taxon)) for (variable in unique(data$rate_type)) {
  a <- data[data$taxon==taxon & data$rate_type==variable, ]
  width <- if(taxon=="Elasmobranchs")50 else 100
  stopifnot(all(a$EV_bin_width_kyr==width))
  e <- environment[environment$EV_bin_width_kyr==width, ]
  keep <- vapply(a$rate_age_Ma, function(age) {
    if (any(abs(e$age_Ma-age)<1e-8)) return(TRUE)
    lo <- e$age_Ma[e$age_Ma<age]; hi <- e$age_Ma[e$age_Ma>age]
    length(lo)>0 && length(hi)>0 && min(hi)-max(lo)<=width/1000+1e-8
  }, logical(1)) & is.finite(a$aligned_variability) & is.finite(a$raw_rate)
  a <- a[keep, ]; pairs[[length(pairs)+1]] <- a
  test <- cor.test(a$aligned_variability, a$raw_rate, method="pearson", alternative="two.sided")
  results[[length(results)+1]] <- data.frame(Taxon=taxon, Biological_variable=variable,
    Environmental_metric="pCO2 variability", Resolution_kyr=width, n=nrow(a),
    Pearson_r=unname(test$estimate), Nominal_P=test$p.value)
}
results <- do.call(rbind, results)
reference <- read.csv(file.path(here, "reference/pCO2_variability_current_S9.csv"))
check <- merge(results, reference, by=c("Taxon", "Biological_variable"), suffixes=c("", "_S9"))
stopifnot(nrow(check)==9)
if (!current_mode) stopifnot(all(check$n==check$n_S9),
          max(abs(check$Pearson_r-check$Pearson_r_S9))<5.1e-7,
          max(abs(check$Nominal_P-check$Nominal_P_S9))<5.1e-7)
dir.create(output, recursive=TRUE, showWarnings=FALSE)
write.csv(results, file.path(output, "pCO2_variability_Pearson_9.csv"), row.names=FALSE)
write.csv(do.call(rbind,pairs), file.path(output, "pCO2_variability_eligible_pairs.csv"), row.names=FALSE)
write.csv(check, file.path(output, "pCO2_variability_S9_comparison.csv"), row.names=FALSE)
print(results)
