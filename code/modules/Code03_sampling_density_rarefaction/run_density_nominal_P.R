#!/usr/bin/env Rscript
# Final S4 nominal P labels: same two-sided Spearman cor.test(exact=FALSE)
# as plot_FigS4_rectangular_rho_P_20261008.R lines 268--303.
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
script <- sub("^--file=", "", args[grepl("^--file=",args)])
root <- dirname(normalizePath(script))
inventory <- read.csv(file.path(root,"reference_results","audited_50kyr_sampling_composition_table.csv"))
inventory <- inventory[order(inventory$bin_center_Ma,decreasing=TRUE),]
reference <- read.csv(file.path(root,"plot_snapshot_inputs","04_figures","Fig_S4_S11_style_rectangular_rho_P_20261008","Spearman_rho_nominal_P_existing_CI.csv"))
predictors <- c(contributing_operational_record_count="contributing_record_count",
 raw_measurement_count_contributing_record_bins_primary="raw_measurement_count_contributing_record_bins",
 raw_measurement_count_all_retained_secondary="raw_measurement_count_all_retained",
 contributing_canonical_site_count="contributing_canonical_site_count")
rows <- list()
for (scale in c("raw_bins","first_difference_older_to_younger")) for (predictor in names(predictors)) {
  x <- inventory[[predictors[[predictor]]]]; y <- inventory$observed_global_EV
  if (scale=="first_difference_older_to_younger") { x<-diff(x);y<-diff(y) }
  ct <- cor.test(x,y,method="spearman",alternative="two.sided",exact=FALSE)
  rows[[length(rows)+1L]] <- data.frame(predictor=predictor,scale=scale,n_pairs=length(x),
    rho=unname(ct$estimate),nominal_P=ct$p.value)
}
results <- do.call(rbind,rows)
joined <- merge(results,reference,by=c("predictor","scale"),suffixes=c("_new","_ref"))
stopifnot(nrow(joined)==8L,all(joined$n_pairs_new==joined$n_pairs_ref),
  max(abs(joined$rho_new-joined$rho_ref))<1e-12,
  max(abs(joined$nominal_P_new-joined$nominal_P_ref))<1e-12)
out <- file.path(root,"outputs","density_nominal_P")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(results,file.path(out,"Spearman_rho_nominal_P_recomputed.csv"),row.names=FALSE)
writeLines("Newly recomputed nominal Spearman P labels; block-bootstrap CIs are separately calculated by run.py.",
  file.path(out,"RUN_STATUS.txt"))
print(results)
