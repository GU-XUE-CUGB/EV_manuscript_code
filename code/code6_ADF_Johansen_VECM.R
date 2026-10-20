# Standalone extinction models using current Dataset S1 EV values.
# Use the supplied biological ages, rates and model settings.
# Run: Rscript code6_ADF_Johansen_VECM.R
args <- gsub("~+~", " ", commandArgs(trailingOnly=FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
output <- file.path(here,'results/generated/code6_current_dataset')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
# 1. Packages
library(urca)

# 2. Input data
elasmo <- data.frame(
  Age_Ma = c(
    67.161106965768795, 67.111097943322093, 67.061088920875406, 67.011079898428704, 66.961070875982003, 66.911061853535301,
    66.8610528310886, 66.811043808641898, 66.761034786195196, 66.711025763748594, 66.661016741301907, 66.611007718855205,
    66.560998696408504, 66.510989673961802, 66.460980651515101, 66.410971629068399, 66.360962606621698, 66.310953584174996,
    66.260944561728394, 66.210935539281706, 66.160926516835005, 66.110917494388303, 66.060908471941602, 66.0108994494949,
    65.960890427048199, 65.910881404601497, 65.860872382154795, 65.810863359708193, 65.760854337261506, 65.710845314814804,
    65.660836292368103, 65.610827269921401, 65.5608182474747, 65.510809225027998, 65.460800202581297, 65.410791180134595,
    65.360782157687893, 65.310773135241305, 65.260764112794604, 65.210755090347902, 65.160746067901201, 65.110737045454499,
    65.060728023007798, 65.010719000561096, 64.960709978114394, 64.910700955667707, 64.860691933221105, 64.810682910774403
  ),
  Global_EV = c(
    0.34776866674282803, 0.25387889431777649, 0.249164599760579, 0.3176766551489425, 0.36728769662804711, 0.25590071657307389,
    0.3668258691388625, 0.33849482218343191, 0.41851287187092823, 0.33531051750415342, 0.43203885703251299, 0.55668716799247198,
    0.51383798291795402, 0.47115855177572052, 0.308595171484496, 0.37116415570289341, 0.45184069249079301, 0.45348930226172141,
    0.52296088508696736, 0.54676823230845317, 0.57480057892402003, 0.5403799049346435, 0.65357703295627112, 0.62469515359352479,
    0.50832176439927057, 0.34298637580210972, 0.38890626284315399, 0.4554065791901013, 0.32984404745080109, 0.3202636831760573,
    0.34140709779352851, 0.33929076417484277, 0.35401315667882599, 0.32102628689682489, 0.33587847497941048, 0.27918043076083471,
    0.4507639666888883, 0.40533209988576208, 0.43116302636828868, 0.30438292364882918, 0.34205551201903089, 0.51623562015493996,
    0.3664183852374257, 0.38797240868934457, 0.40629271835756359, 0.37095493966002191, 0.36155600869747068, 0.36462488948932359
  ),
  Extinction_rate = c(
    0.16899212294930599, 0.16899548072069101, 0.169027777442706, 0.169027777442706, 0.16904082143364799, 0.170162603327777,
    0.17024465515064499, 0.17024465515064499, 0.17024465515064499, 0.17024465515064499, 0.177863466229981, 0.177863466229981,
    1.3392314743548901, 1.4017759881335701, 1.4017759881335701, 1.5969067229186, 1.59699466927849, 1.67558343563992,
    1.9491105912795801, 1.9491105912795801, 1.9491105912795801, 1.9491105912795801, 1.9489977437986401, 1.9462604915021899,
    1.9462604915021899, 1.9460040449135301, 1.8172981651860101, 1.8021194715606099, 1.11847583504244, 0.17871583247963899,
    0.16585872424296499, 0.15914053498304001, 0.15914053498304001, 0.15914053498304001, 0.15299730242783599, 0.147979165780219,
    0.147979165780219, 0.14793092698417201, 0.14793092698417201, 0.14795354810026501, 0.14795354810026501, 0.14795354810026501,
    0.14793905805635099, 0.14793905805635099, 0.14793905805635099, 0.147974997033848, 0.14802221455319201, 0.14802221455319201
  )
)
bivalve <- data.frame(
  Age_Ma = c(
    67.195792212694798, 67.095657916481002, 66.995523620267207, 66.895389324053397, 66.795255027839602, 66.695120731625806,
    66.594986435411997, 66.494852139198201, 66.394717842984406, 66.294583546770596, 66.194449250556701, 66.094314954342906,
    65.994180658129096, 65.894046361915301, 65.793912065701505, 65.693777769487696, 65.5936434732739, 65.493509177060105,
    65.393374880846295, 65.2932405846325, 65.193106288418704, 65.092971992204895, 64.992837695991, 64.892703399777204
  ),
  Global_EV = c(
    0.38224678476614199, 0.307201416280385, 0.3559935489808041, 0.31516440774368298, 0.44162368163752458, 0.36537801367285472,
    0.52610191209248092, 0.45054345141240848, 0.45494433594452471, 0.68640386210972326, 0.68106613846801023, 0.65252131196023455,
    0.76576484236027542, 0.44804134973229981, 0.50861183577136182, 0.41707399807513029, 0.3687001622386446, 0.47160569335645142,
    0.38663454419130661, 0.43264566834511659, 0.40249971542578128, 0.481351977077007, 0.39818257767945869, 0.4504431834761991
  ),
  Extinction_rate = c(
    0.15058198604892301, 0.150993584568614, 0.15172991693173499, 0.15551358713274099, 0.157188483605116, 0.15984380704295301,
    0.29279218521230899, 0.29279218521230899, 0.39776236598585601, 0.54658067341239802, 1.38250326345032, 2.2902388266938698,
    2.2902388266938698, 2.2902388266938698, 2.03918180717847, 2.0387145695070599, 0.23809690358917501, 0.223387267722684,
    0.223387267722684, 0.22000786910734599, 0.21932522709333499, 0.21933295057529301, 0.21933295057529301, 0.21936171488834499
  )
)

series <- list(
  Elasmobranchs = list(data = elasmo, EV_resolution = 50L),
  Bivalves = list(data = bivalve, EV_resolution = 100L)
)

# 3. ADF tests
integration_order <- function(x) {
  level <- ur.df(x, type = "none", selectlags = "AIC", lags = 1)
  first_difference <- ur.df(diff(x), type = "none", selectlags = "AIC", lags = 1)
  if (level@teststat[1] < level@cval[1, "5pct"]) return("I(0)")
  if (first_difference@teststat[1] < first_difference@cval[1, "5pct"]) return("I(1)-compatible")
  "Unresolved"
}

# 4. Johansen cointegration
rank_at <- function(j, threshold) {
  if (j@teststat[2] <= j@cval[2, threshold]) return(0L)
  if (j@teststat[1] <= j@cval[1, threshold]) return(1L)
  2L
}

# 5. VECM and ECT
details <- list()
analyse <- function(label, item) {
  d <- item$data
  stopifnot(all(diff(d$Age_Ma) < 0), !anyNA(d), nrow(d) >= 6L)
  y <- cbind(EV = d$Global_EV, Rate = d$Extinction_rate)
  if (integration_order(y[, "EV"]) != "I(1)-compatible" ||
      integration_order(y[, "Rate"]) != "I(1)-compatible") {
    stop(paste("I(1) prerequisite not met:", label))
  }
  j <- ca.jo(y, type = "trace", ecdet = "none", K = 2, spec = "longrun")
  rank_5 <- rank_at(j, "5pct")
  rank_10 <- rank_at(j, "10pct")
  if (rank_5 != 1L || rank_10 != 1L) stop(paste("Rank-1 prerequisite not met:", label))
  v <- cajorls(j, r = 1)
  stopifnot(identical(colnames(coef(v$rlm)), c("EV.d", "Rate.d")),
            abs(v$beta[1, 1] - 1) < 1e-12)
  sm <- summary(v$rlm)
  names(sm) <- colnames(coef(v$rlm))
  coefficient <- function(equation, column) sm[[equation]]$coefficients["ect1", column]
  beta <- as.vector(v$beta[, 1])
  forward <- directional_tests(y, beta, "Rate", "EV.dl1")
  reverse <- directional_tests(y, beta, "EV", "Rate.dl1")
  stopifnot(max(abs(forward$coefficients-coef(v$rlm)[,"Rate.d"]))<1e-9,
            max(abs(reverse$coefficients-coef(v$rlm)[,"EV.d"]))<1e-9,
            max(abs(forward$SE-sm[["Rate.d"]]$coefficients[,2]))<1e-9,
            max(abs(reverse$SE-sm[["EV.d"]]$coefficients[,2]))<1e-9,
            abs(forward$short_run["F"]-sm[["Rate.d"]]$coefficients["EV.dl1",3]^2)<1e-8,
            abs(reverse$short_run["F"]-sm[["EV.d"]]$coefficients["Rate.dl1",3]^2)<1e-8)
  adf <- do.call(rbind,lapply(colnames(y),function(variable) {
    do.call(rbind,lapply(c("Level","First difference"),function(transformation) {
      x <- if(transformation=="Level")y[,variable]else diff(y[,variable])
      a <- ur.df(x,type="none",selectlags="AIC",lags=1)
      data.frame(Taxon=label,Variable=variable,Transformation=transformation,
        Statistic=as.numeric(a@teststat[1]),CV_1pct=a@cval[1,"1pct"],
        CV_5pct=a@cval[1,"5pct"],CV_10pct=a@cval[1,"10pct"],
        Integration_status=integration_order(y[,variable]))
    }))
  }))
  trace <- data.frame(Taxon=label,EV_resolution_kyr=item$EV_resolution,n=nrow(d),
    Hypothesis=c("r = 0","r <= 1"),Trace_statistic=j@teststat[c(2,1)],
    CV_1pct=j@cval[c(2,1),"1pct"],CV_5pct=j@cval[c(2,1),"5pct"],
    CV_10pct=j@cval[c(2,1),"10pct"],rank_5pct=rank_5,rank_10pct=rank_10,
    type="trace",ecdet="none",K=2,spec="longrun")
  terms <- do.call(rbind,lapply(names(sm),function(equation) {
    a <- sm[[equation]]$coefficients
    data.frame(Taxon=label,Equation=equation,Term=rownames(a),
      Coefficient=a[,1],SE=a[,2],t_statistic=a[,3],P=a[,4],row.names=NULL)
  }))
  tests <- do.call(rbind,lapply(c("forward","reverse"),function(direction) {
    a <- if(direction=="forward")forward else reverse
    do.call(rbind,lapply(c("short_run","joint"),function(kind) {
      z <- a[[kind]]
      data.frame(Taxon=label,Direction=if(direction=="forward")"EV to Extinction"else "Extinction to EV",
        Test=kind,F=unname(z["F"]),df1=if(kind=="short_run")1 else 2,
        df2=a$df,P=unname(z["P"]))
    }))
  }))
  details[[label]] <<- list(ADF=adf,Johansen=trace,Coefficients=terms,Directional=tests,
    Beta=data.frame(Taxon=label,Term=rownames(v$beta),Coefficient=beta))
  data.frame(
    Taxon = label, n = nrow(d), EV_resolution = item$EV_resolution,
    rank_5pct = rank_5, rank_10pct = rank_10,
    ECT_Extinction = coefficient("Rate.d", 1),
    P_ECT_Extinction = coefficient("Rate.d", 4),
    ECT_EV = coefficient("EV.d", 1),
    P_ECT_EV = coefficient("EV.d", 4),
    ShortRun_EV_to_Extinction_F = forward$short_run["F"],
    P_EV_to_Extinction = forward$short_run["P"],
    ShortRun_Extinction_to_EV_F = reverse$short_run["F"],
    P_Extinction_to_EV = reverse$short_run["P"],
    Joint_EV_to_Extinction_F = forward$joint["F"],
    P_Joint_EV_to_Extinction = forward$joint["P"],
    Joint_Extinction_to_EV_F = reverse$joint["F"],
    P_Joint_Extinction_to_EV = reverse$joint["P"],
    row.names = NULL
  )
}

# 6. Directional tests
# Granger predictive information does not establish physical causation.
directional_tests <- function(y, beta, response_variable, cross_lag) {
  n <- nrow(y)
  level_lag2 <- y[seq_len(n - 2L), , drop = FALSE]
  change_lag1 <- y[2:(n - 1L), , drop = FALSE] - level_lag2
  change_now <- y[3:n, , drop = FALSE] - y[2:(n - 1L), , drop = FALSE]
  X <- cbind(
    ect1 = as.vector(level_lag2 %*% beta),
    constant = 1,
    EV.dl1 = change_lag1[, "EV"],
    Rate.dl1 = change_lag1[, "Rate"]
  )
  target <- change_now[, response_variable]
  full <- lm.fit(X, target)
  df <- nrow(X) - ncol(X)
  rss <- sum(full$residuals^2)
  test <- function(terms) {
    reduced <- lm.fit(X[, setdiff(colnames(X), terms), drop = FALSE], target)
    q <- length(terms)
    F_value <- ((sum(reduced$residuals^2) - rss) / q) / (rss / df)
    c(F = F_value, P = pf(F_value, q, df, lower.tail = FALSE))
  }
  list(short_run = test(cross_lag), joint = test(c("ect1", cross_lag)),
       coefficients=full$coefficients,
       SE=sqrt(diag(chol2inv(chol(crossprod(X))))*rss/df),df=df)
}

# 7. Output
summary_table <- do.call(rbind, lapply(names(series), function(name) analyse(name, series[[name]])))
rownames(summary_table) <- NULL
print(summary_table, row.names = FALSE, digits = 9)
write.csv(summary_table, file.path(output,"code10_summary_results.csv"), row.names = FALSE)
write.csv(elasmo,file.path(output,'Embedded_elasmobranchs_matched_pairs_48.csv'),row.names=FALSE)
write.csv(bivalve,file.path(output,'Embedded_bivalves_matched_pairs_24.csv'),row.names=FALSE)
for(kind in c("ADF","Johansen","Coefficients","Directional","Beta")) {
  table <- do.call(rbind,lapply(details,function(x)x[[kind]]))
  write.csv(table,file.path(output,paste0("Code10_",kind,".csv")),row.names=FALSE)
}
if(file.exists(file.path(here,'reference/code6_unified_results.csv'))) {
reference <- subset(read.csv(file.path(here,'reference/code6_unified_results.csv')),Biological_variable=='Extinction' & Taxon%in%c('Elasmobranchs','Bivalves'))
comparison <- merge(summary_table,reference,by='Taxon')
map <- c(ECT_Extinction='ECT_Rate',P_ECT_Extinction='P_ECT_Rate',ECT_EV='ECT_EV',P_ECT_EV='P_ECT_EV',
         P_EV_to_Extinction='P_EV_to_Rate',P_Extinction_to_EV='P_Rate_to_EV',P_Joint_EV_to_Extinction='P_Joint_EV_to_Rate')
stopifnot(nrow(comparison)==2L,all(comparison$n.x==comparison$n.y))
for(column in names(map)) {
  a <- if(column%in%names(reference))paste0(column,'.x')else column
  b <- if(map[column]%in%names(summary_table))paste0(map[column],'.y')else map[column]
  stopifnot(max(abs(comparison[[a]]-comparison[[b]]))<1e-9)
}
write.csv(comparison,file.path(output,'Extinction_reference_comparison.csv'),row.names=FALSE)
}
capture.output(sessionInfo(),file=file.path(output,'R_sessionInfo.txt'))
