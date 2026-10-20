# Existing age-matched values only; no new interpolation or reconstruction.
args <- gsub("~+~", " ", commandArgs(trailingOnly=FALSE), fixed=TRUE)
root <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
src<-file.path(root,'inputs')
out<-file.path(root,'results/ancillary_existing_stages')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files<-file.path(src,c('variability_biotic_paired_data_final.csv','final_environmental_variability_inputs.csv'))
hash<-tools::md5sum(files)
metrics<-c('Global EV','Temperature variability','d18O variability','d13C variability','Weathering variability')
d<-subset(read.csv(files[1]),EV_version=='Final revised' & predictor %in% metrics)
env<-read.csv(files[2]);results<-list();audit<-list()
for(taxon in unique(d$taxon))for(variable in unique(d$rate_type))for(metric in metrics){
 a<-d[d$taxon==taxon & d$rate_type==variable & d$predictor==metric,]
 width<-unique(a$EV_bin_width_kyr);stopifnot(length(width)==1,width==if(taxon=='Elasmobranchs')50 else 100)
 e<-env[env$predictor==metric & env$EV_bin_width_kyr==width & is.finite(env$variability),];e<-e[order(e$age_Ma),]
 stopifnot(!anyDuplicated(e$age_Ma))
 reason<-vapply(a$rate_age_Ma,function(age){
  if(any(abs(e$age_Ma-age)<1e-8))return('Retained: existing exact-age bin')
  lo<-e$age_Ma[e$age_Ma<age];hi<-e$age_Ma[e$age_Ma>age]
  if(!length(lo)||!length(hi))return('Excluded: existing endpoint carry')
  if(min(hi)-max(lo)>width/1000+1e-8)return('Excluded: existing interpolation across missing bins')
  'Retained: existing alignment between adjacent observed bins'
 },character(1))
 keep<-startsWith(reason,'Retained') & is.finite(a$raw_rate) & is.finite(a$aligned_variability)
 audit[[length(audit)+1]]<-data.frame(Taxon=taxon,Biological_variable=variable,Environmental_metric=metric,Resolution_kyr=width,Actual_age_Ma=a$rate_age_Ma,Retained=keep,Reason=reason)
 a<-a[keep,]
 stages<-list(Pre_boundary=a$rate_age_Ma>66,Boundary_early_recovery=a$rate_age_Ma<=66 & a$rate_age_Ma>=65.5,Post_recovery_background=a$rate_age_Ma<65.5,Full_transition=rep(TRUE,nrow(a)))
 stopifnot(all(Reduce('+',stages[1:3])==1))
 for(stage in names(stages)){
  z<-a[stages[[stage]],];n<-nrow(z);r<-p<-NA_real_
  valid<-n>=3 && sd(z$raw_rate)>0 && sd(z$aligned_variability)>0
  if(valid){fit<-cor.test(z$aligned_variability,z$raw_rate,method='pearson');r<-unname(fit$estimate);p<-fit$p.value}
  results[[length(results)+1]]<-data.frame(Taxon=taxon,Biological_variable=variable,Stage=stage,Environmental_metric=metric,Resolution_kyr=width,n=n,Pearson_r=r,Nominal_P=p,Actual_age_range_Ma=if(n)sprintf('%.8f-%.8f',max(z$rate_age_Ma),min(z$rate_age_Ma))else NA_character_,Status=if(valid)'Nominal P; no serial-dependence or multiplicity correction'else 'Not estimable: n<3 or constant series')
 }
}
r<-do.call(rbind,results);stopifnot(nrow(r)==180)
write.csv(r,file.path(out,'Stage_proxy_Pearson_results.csv'),row.names=FALSE,na='NA')
x<-subset(r,Biological_variable=='Extinction')
wide<-reshape(x[,c('Taxon','Stage','Resolution_kyr','Environmental_metric','Pearson_r','n')],idvar=c('Taxon','Stage','Resolution_kyr'),timevar='Environmental_metric',direction='wide')
write.csv(wide,file.path(out,'Extinction_proxy_stage_summary.csv'),row.names=FALSE,na='NA')
write.csv(do.call(rbind,audit),file.path(out,'Existing_pair_eligibility_audit.csv'),row.names=FALSE)
writeLines(c(paste(names(hash),hash),
 'Reuse prior matched environmental values unchanged. No biological or environmental interpolation is performed in this script.',
 'Retain exact observed environmental bin matches or existing matches between two adjacent observed bins. Exclude pre-existing endpoint carry and interpolation spanning missing environmental bins, identically for all five metrics.',
 'Stage membership uses actual rate_age_Ma: >66; 65.5<=age<=66; <65.5. Full uses all eligible existing pairs.',
 'Consequently n and Global EV correlations can differ from earlier diagnostics that retained endpoint carry. Metrics may also have different available n. No formal test of differences between correlations is performed.',
 'The 66.0 stage cutoff is not the physical K-Pg age 66.016 Ma.',
 'Calcareous nannofossil biological series is a 100-kyr spline grid derived from approximately 1-Myr source data; nominal P is not independent high-resolution evidence.',
 'Only two-sided Pearson tests; no Spearman, cointegration, VECM, bootstrap, simulation or figures.'),file.path(out,'Methods_and_provenance.txt'))
stopifnot(identical(hash,tools::md5sum(files)))
reference <- read.csv(file.path(root,'reference/Stage_proxy_Pearson_results.csv'))
comparison <- merge(r,reference,by=c('Taxon','Biological_variable','Stage','Environmental_metric'),
                    suffixes=c('_recomputed','_reference'))
stopifnot(nrow(comparison)==180,all(comparison$n_recomputed==comparison$n_reference),
          identical(is.na(comparison$Pearson_r_recomputed),is.na(comparison$Pearson_r_reference)),
          max(abs(comparison$Pearson_r_recomputed-comparison$Pearson_r_reference),na.rm=TRUE)<1e-12,
          max(abs(comparison$Nominal_P_recomputed-comparison$Nominal_P_reference),na.rm=TRUE)<1e-12)
write.csv(comparison,file.path(out,'Reference_comparison_180.csv'),row.names=FALSE)
print(wide,row.names=FALSE)
