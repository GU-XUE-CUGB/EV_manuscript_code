# 1. Packages and primary model functions
args <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1])))
source(file.path(here,"code6_ADF_Johansen_VECM.R"))

# 2. Input data
file_path <- file.path(here,"modules/Code08_EV_proxy_biota_Pearson/inputs/current_variability_biotic_pairs.csv")
data <- subset(read.csv(file_path),EV_version=="Final revised" & predictor=="Global EV")

# 3. ADF and Johansen screening
results <- list()
for(taxon in unique(data$taxon)) for(variable in unique(data$rate_type)) {
  d <- data[data$taxon==taxon & data$rate_type==variable & is.finite(data$aligned_variability),]
  d <- d[order(d$rate_age_Ma,decreasing=TRUE),]
  stopifnot(all(diff(d$rate_age_Ma)<0))
  ev_order <- integration_order(d$aligned_variability)
  rate_order <- integration_order(d$raw_rate)
  row <- data.frame(Taxon=taxon,Biological_variable=variable,n=nrow(d),
    EV_resolution=if(taxon=="Elasmobranchs")50 else 100,
    EV_integration=ev_order,Rate_integration=rate_order,
    Trace_r0=NA_real_,Trace_r_le1=NA_real_,rank_5pct=NA_integer_,rank_10pct=NA_integer_)
  if(ev_order=="I(1)-compatible" && rate_order=="I(1)-compatible") {
    y <- cbind(EV=d$aligned_variability,Rate=d$raw_rate)
    j <- ca.jo(y,type="trace",ecdet="none",K=2,spec="longrun")
    row$Trace_r0 <- j@teststat[2];row$Trace_r_le1 <- j@teststat[1]
    row$rank_5pct <- rank_at(j,"5pct");row$rank_10pct <- rank_at(j,"10pct")
  }
  results[[length(results)+1]] <- row
}

# 4. Output
results <- do.call(rbind,results)
write.csv(results,file.path(output,"Johansen_all_rates_current_dataset.csv"),row.names=FALSE)
print(results)
