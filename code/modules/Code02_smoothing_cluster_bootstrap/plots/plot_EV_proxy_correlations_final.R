#!/usr/bin/env Rscript
# Run via recalculate_EV_proxy_correlations_final.py after extracting final XLSX inputs.
# Same Pearson cor.test / lm framework as existing analyses; no smoothing.
suppressPackageStartupMessages({library(svglite)})
arg <- gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=','',arg[grepl('^--file=',arg)]))))
out <- file.path(root,'reference_results','current_dataset')
fig <- file.path(root,'outputs','plots','EV_proxy_correlations_current_dataset')
dir.create(fig, recursive=TRUE, showWarnings=FALSE)
d <- read.csv(file.path(out,'paired_data_all_12_tests.csv'))
r <- read.csv(file.path(out,'EV_proxy_Pearson_12_tests_recomputed.csv'))
# Plot the saved current estimates without recomputing or overwriting statistics.
colors <- c('#DF0028','#FA673D','#FFA044','#FFDA5A')
proxies <- c('temperature','d18O','d13C','weathering')
phases <- c('Full interval','Phase 1','Phase 2')
labs <- list(
  expression(bold(Temperature~variability)),
  expression(bold(delta^{18}*O~variability)),
  expression(bold(delta^{13}*C~variability)),
  expression(bold(Weathering~variability))
)
titles <- c('A   67.2-64.8 Ma','B   66.9-66.45 Ma (Phase 1)','C   66.4-66.05 Ma (Phase 2)')
draw <- function(family='Arial'){
 par(mfrow=c(4,3),mar=c(3.5,3.8,2.8,1),oma=c(.3,.3,1,0),family=family,mgp=c(2.2,.65,0),tcl=-.25,cex=.95)
 for(j in seq_along(proxies)) for(k in seq_along(phases)){
  z<-d[d$phase==phases[k]&d$proxy==proxies[j],];rr<-r[r$phase==phases[k]&r$proxy==proxies[j],]
  # y is explicitly Global EV; x is the corresponding proxy-specific mean IQR.
  xl<-c(min(.1,floor(min(z$proxy_EV)*10)/10),max(.8,ceiling(max(z$proxy_EV)*10)/10+.05))
  plot(z$proxy_EV,z$global_EV,pch=16,col=colors[j],cex=1.05,xlim=xl,ylim=c(0,1),xaxs='i',yaxs='i',
       xlab=labs[[j]],ylab=expression(bold(Global~EV)),font.lab=2,yaxt='n',bty='o')
  axis(2,at=seq(0,1,.2),las=1)
  xx<-range(z$proxy_EV);lines(xx,rr$intercept+rr$slope*xx,lty=3,lwd=1,col='#333333')
  text(xl[1]+.035*diff(xl),.965,labels=sprintf("Pearson's r = %.2f\nP = %.3g",rr$pearson_r,rr$p_value),adj=c(0,1),cex=.95)
  if(j==1)mtext(titles[k],side=3,line=1.6,adj=0,font=2,cex=1.03)
 }
}
stem<-file.path(fig,'Fig_EV_proxy_correlations_final_50kyr')
pdf(paste0(stem,'.pdf'),width=12,height=12.2,family='Helvetica',useDingbats=FALSE);draw('Helvetica');dev.off()
svglite(paste0(stem,'.svg'),width=12,height=12.2);draw();dev.off()
ragg::agg_png(paste0(stem,'.png'),width=12,height=12.2,units='in',res=600);draw();dev.off()
ragg::agg_png(file.path(fig,'preview.png'),width=12,height=12.2,units='in',res=120);draw();dev.off()
writeLines(capture.output(sessionInfo()),file.path(fig,'R_sessionInfo.txt'))
