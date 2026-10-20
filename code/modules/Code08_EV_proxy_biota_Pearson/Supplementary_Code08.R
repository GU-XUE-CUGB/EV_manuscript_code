# Final full-transition Pearson analysis, adapted from stage_proxy_Pearson_20260922.R.
# Run: Rscript /path/to/Code08_EV_proxy_biota_Pearson/Supplementary_Code08.R
# Previously aligned EV values are locked input; this script creates no interpolation.
args <- gsub("~+~", " ", commandArgs(trailingOnly=FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
cfg <- list(paired=file.path(here,'inputs/variability_biotic_paired_data_final.csv'),
            environmental=file.path(here,'inputs/final_environmental_variability_inputs.csv'),
            reference=file.path(here,'reference/Stage_proxy_Pearson_results.csv'),
            output=file.path(here,'results/generated'))
current_mode <- '--current-inputs' %in% commandArgs(TRUE)
if (current_mode) {
  cfg$paired <- file.path(here,'inputs/current_variability_biotic_pairs.csv')
  cfg$environmental <- file.path(here,'inputs/current_environmental_variability_inputs.csv')
  cfg$output <- file.path(here,'results/generated/current_dataset')
}
dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
hash <- tools::md5sum(unlist(cfg[1:3]))
metrics <- c('Global EV','Temperature variability','d18O variability','d13C variability','Weathering variability')
d <- subset(read.csv(cfg$paired),EV_version=='Final revised' & predictor %in% metrics)
env <- read.csv(cfg$environmental)
results <- audit <- paired <- list()
for(taxon in unique(d$taxon)) for(variable in unique(d$rate_type)) for(metric in metrics) {
  a <- d[d$taxon==taxon & d$rate_type==variable & d$predictor==metric,]
  width <- unique(a$EV_bin_width_kyr)
  stopifnot(length(width)==1,width==if(taxon=='Elasmobranchs')50 else 100)
  e <- env[env$predictor==metric & env$EV_bin_width_kyr==width & is.finite(env$variability),]
  e <- e[order(e$age_Ma),]; stopifnot(!anyDuplicated(e$age_Ma))
  reason <- vapply(a$rate_age_Ma,function(age) {
    if(any(abs(e$age_Ma-age)<1e-8)) return('Retained: existing exact-age bin')
    lo <- e$age_Ma[e$age_Ma<age]; hi <- e$age_Ma[e$age_Ma>age]
    if(!length(lo)||!length(hi)) return('Excluded: existing endpoint carry')
    if(min(hi)-max(lo)>width/1000+1e-8) return('Excluded: existing interpolation across missing bins')
    'Retained: existing alignment between adjacent observed bins'
  },character(1))
  keep <- startsWith(reason,'Retained') & is.finite(a$raw_rate) & is.finite(a$aligned_variability)
  audit[[length(audit)+1]] <- data.frame(Taxon=taxon,Biological_variable=variable,
    Environmental_metric=metric,Resolution_kyr=width,Actual_age_Ma=a$rate_age_Ma,Retained=keep,Reason=reason)
  a <- a[keep,]; paired[[length(paired)+1]] <- a
  fit <- cor.test(a$aligned_variability,a$raw_rate,method='pearson',alternative='two.sided')
  results[[length(results)+1]] <- data.frame(Taxon=taxon,Biological_variable=variable,
    Stage='Full_transition',Environmental_metric=metric,Resolution_kyr=width,n=nrow(a),
    Pearson_r=unname(fit$estimate),Nominal_P=fit$p.value,
    Actual_age_range_Ma=sprintf('%.8f-%.8f',max(a$rate_age_Ma),min(a$rate_age_Ma)),
    Status='Nominal P; no serial-dependence or multiplicity correction')
}
r <- do.call(rbind,results); stopifnot(nrow(r)==45)
ref <- subset(read.csv(cfg$reference),Stage=='Full_transition')
key <- c('Taxon','Biological_variable','Stage','Environmental_metric')
z <- merge(r,ref,by=key,suffixes=c('_recomputed','_reference'))
stopifnot(nrow(z)==45,identical(hash,tools::md5sum(unlist(cfg[1:3]))))
if (!current_mode) stopifnot(all(z$n_recomputed==z$n_reference),
          max(abs(z$Pearson_r_recomputed-z$Pearson_r_reference))<1e-12,
          max(abs(z$Nominal_P_recomputed-z$Nominal_P_reference))<1e-12)
write.csv(r,file.path(cfg$output,'Full_transition_Pearson_45.csv'),row.names=FALSE)
write.csv(do.call(rbind,paired),file.path(cfg$output,'Final_eligible_EV_biota_pairs.csv'),row.names=FALSE)
write.csv(do.call(rbind,audit),file.path(cfg$output,'Existing_pair_eligibility_audit.csv'),row.names=FALSE)
write.csv(z,file.path(cfg$output,'Reference_comparison.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(cfg$output,'R_sessionInfo.txt'))
print(r[,c(1,2,4,5,6,7,8)],row.names=FALSE)
