# Mean states and full-transition Pearson tests, adapted from update_mean_state_FigS13_20260922.R.
# Run: Rscript /path/to/Code09_mean_state_biota_Pearson/Supplementary_Code09.R
# Core Fig S13: global plus four proxy classes (45 tests). Separate pCO2: nine ancillary tests.
args <- gsub("~+~", " ", commandArgs(trailingOnly=FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
cfg <- list(raw=file.path(here,'inputs/retained_measurements_with_zscores_v4.csv'),
            biology=file.path(here,'inputs/biology'),reference=file.path(here,'reference'),
            output=file.path(here,'results/generated'))
dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
metrics <- c('Global mean environmental state','Temperature mean state','d13C mean state','d18O mean state','Weathering mean state','pCO2 mean state')
categories <- c('global','temperature','d13C','d18O','weathering','pCO2')
rates <- c('Speciation','Extinction','Net diversification')
sources <- list(Elasmobranchs=c('Elasmobranchs_time.csv','Elasmobranchs_Speciation.csv','Elasmobranchs_Extinction.csv','Elasmobranchs_Net_diversification.csv'),
 Bivalves=c('Bivalves_time.csv','Bivalves_Speciation.csv','Bivalves_Extinction.csv','Bivalves_Net_diversification.csv'),
 'Calcareous nannofossils'=c('Nannofossils_time.csv','Nannofossils_Speciation.csv','Nannofossils_Extinction.csv','Nannofossils_Net_diversification.csv'))
paths <- c(cfg$raw,unlist(lapply(sources,function(x)file.path(cfg$biology,x))))
hash <- tools::md5sum(paths)
m <- read.csv(cfg$raw,check.names=FALSE)
stopifnot(nrow(m)==10881,length(unique(m$record_id))==70,all(is.finite(m$value)),all(is.finite(m$age_Ma)))
m$z_new <- ave(m$value,m$record_id,FUN=function(v)(v-mean(v))/sd(v))
stopifnot(max(abs(m$z_new-m$z_score_recomputed))<1e-9)
states <- list()
for(i in seq_along(metrics)) {
  a <- if(i==1)m else m[m$category==categories[i],]
  for(w in c(50,100)) for(mid in seq(64800,67200,by=w)) {
    # Literal established bin arithmetic in Ma; retain floating boundary behavior.
    start <- if(w==50)64.775 else 64.75
    bi <- round((mid-64800)/w); lo <- start+bi*(w/1000); hi <- start+(bi+1)*(w/1000)
    keep <- a$age_Ma>=lo & a$age_Ma<hi; v <- a$z_new[keep]
    states[[length(states)+1]] <- data.frame(Mean_environmental_metric=metrics[i],
      Resolution_kyr=w,Midpoint_Ma=mid/1000,Mean_state=if(length(v))mean(v)else NA_real_,
      Measurement_n=length(v),Record_n=length(unique(a$record_id[keep])))
  }
}
states <- do.call(rbind,states)
for(w in c(50,100)) {
  actual <- subset(states,Resolution_kyr==w & Mean_environmental_metric%in%metrics[1:5])
  ref <- read.csv(file.path(cfg$reference,paste0('Revised_mean_states_',w,'kyr.csv')))
  rownames(actual) <- rownames(ref) <- NULL
  stopifnot(nrow(actual)==nrow(ref),
            all(actual$Mean_environmental_metric==ref$Mean_environmental_metric),
            all(actual$Resolution_kyr==ref$Resolution_kyr),
            max(abs(actual$Midpoint_Ma-ref$Midpoint_Ma))<1e-12,
            identical(is.na(actual$Mean_state),is.na(ref$Mean_state)),
            max(abs(actual$Mean_state-ref$Mean_state),na.rm=TRUE)<1e-12,
            all(actual$Measurement_n==ref$Measurement_n),all(actual$Record_n==ref$Record_n))
}
match_state <- function(age,e,w) {
  j <- which(abs(e$Midpoint_Ma-age)<1e-8)
  if(length(j)) return(c(value=e$Mean_state[j[1]],lo=e$Midpoint_Ma[j[1]],hi=e$Midpoint_Ma[j[1]]))
  lo <- which(e$Midpoint_Ma<age); hi <- which(e$Midpoint_Ma>age)
  if(!length(lo)||!length(hi)) return(c(value=NA,lo=NA,hi=NA))
  l <- tail(lo,1); h <- hi[1]
  if(e$Midpoint_Ma[h]-e$Midpoint_Ma[l]>w/1000+1e-8 || !all(is.finite(e$Mean_state[c(l,h)])))
    return(c(value=NA,lo=NA,hi=NA))
  c(value=approx(e$Midpoint_Ma[c(l,h)],e$Mean_state[c(l,h)],xout=age,rule=1)$y,
    lo=e$Midpoint_Ma[l],hi=e$Midpoint_Ma[h])
}
res <- pairs <- list()
for(taxon in names(sources)) {
  v <- lapply(file.path(cfg$biology,sources[[taxon]]),function(f)scan(f,sep=',',quiet=TRUE))
  age <- v[[1]]; stopifnot(!anyDuplicated(age),all(lengths(v)==length(age)))
  w <- if(taxon=='Elasmobranchs')50 else 100
  for(metric in metrics) {
    e <- states[states$Mean_environmental_metric==metric & states$Resolution_kyr==w,]
    e <- e[order(e$Midpoint_Ma),]; matched <- t(vapply(age,match_state,numeric(3),e=e,w=w))
    for(k in 1:3) {
      y <- v[[k+1]]; z <- as.numeric(scale(y)); keep <- is.finite(matched[,1]) & is.finite(z)
      ct <- cor.test(matched[keep,1],z[keep],method='pearson',alternative='two.sided')
      res[[length(res)+1]] <- data.frame(Taxon=taxon,Biological_variable=rates[k],
        Mean_environmental_metric=metric,Resolution_kyr=w,n=sum(keep),Pearson_r=unname(ct$estimate),
        Nominal_P=ct$p.value,Actual_age_range_Ma=sprintf('%.8f-%.8f',max(age[keep]),min(age[keep])))
      pairs[[length(pairs)+1]] <- data.frame(Taxon=taxon,Biological_variable=rates[k],
        Mean_environmental_metric=metric,Resolution_kyr=w,Age_Ma=age,Biological_rate=y,
        Standardized_rate=z,Matched_mean_state=matched[,1],Environmental_lower_midpoint=matched[,2],
        Environmental_upper_midpoint=matched[,3],Included=keep)
    }
  }
}
res <- do.call(rbind,res); pairs <- do.call(rbind,pairs)
core <- subset(res,Mean_environmental_metric%in%metrics[1:5]); stopifnot(nrow(core)==45,nrow(res)==54)
ref <- read.csv(file.path(cfg$reference,'Mean_state_biota_Pearson_results.csv'))
keys <- c('Taxon','Biological_variable','Mean_environmental_metric')
comparison <- merge(core,ref,by=keys,suffixes=c('_recomputed','_reference'))
stopifnot(nrow(comparison)==45,all(comparison$n_recomputed==comparison$n_reference),
          max(abs(comparison$Pearson_r_recomputed-comparison$Pearson_r_reference))<1e-12,
          max(abs(comparison$Nominal_P_recomputed-comparison$Nominal_P_reference))<1e-12,
          identical(hash,tools::md5sum(paths)))
# The nine pCO2 values in the final S9 export are verified, never replayed as calculations.
co2ref <- read.csv(file.path(cfg$reference,'pCO2_mean_state_final_9.csv'))
co2cmp <- merge(subset(res,Mean_environmental_metric==metrics[6]),co2ref,by=keys,suffixes=c('_recomputed','_reference'))
stopifnot(nrow(co2cmp)==9,all(co2cmp$n_recomputed==co2cmp$n_reference),
          max(abs(co2cmp$Pearson_r_recomputed-co2cmp$Pearson_r_reference))<1e-12,
          max(abs(co2cmp$Nominal_P_recomputed-co2cmp$Nominal_P_reference))<1e-12)
write.csv(states,file.path(cfg$output,'Recomputed_mean_states_all_6.csv'),row.names=FALSE,na='NA')
write.csv(core,file.path(cfg$output,'Mean_state_biota_Pearson_results.csv'),row.names=FALSE)
write.csv(res,file.path(cfg$output,'Ancillary_complete_mean_state_Pearson_54.csv'),row.names=FALSE)
write.csv(pairs,file.path(cfg$output,'Mean_state_biota_paired_observations_all_6.csv'),row.names=FALSE,na='NA')
write.csv(comparison,file.path(cfg$output,'Core_reference_comparison.csv'),row.names=FALSE)
write.csv(co2cmp,file.path(cfg$output,'pCO2_reference_comparison.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(cfg$output,'R_sessionInfo.txt'))
print(res[,c(1:7)],row.names=FALSE)
