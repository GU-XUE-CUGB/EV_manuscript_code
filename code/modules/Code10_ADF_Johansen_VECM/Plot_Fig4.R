# Fig. 4 redraw only: all three rate variables in A, extinction-only
# directional F tests in B. Reads previously finalized statistics.
suppressPackageStartupMessages({library(ggplot2);library(grid)})
if (.Platform$OS.type=='unix') invisible(try(Sys.setlocale('LC_CTYPE','en_US.UTF-8'),silent=TRUE))
args <- gsub("~+~", " ", commandArgs(trailingOnly=FALSE), fixed=TRUE)
here <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
cfg <- list(
  johansen=file.path(here,'../../results/code6_current_dataset/Johansen_all_rates_current_dataset.csv'),
  directional=file.path(here,'../../results/code6_current_dataset/Code10_Directional.csv'),
  out=file.path(here,'../../results/generated/figures_from_frozen_code6_snapshot')
)
stopifnot(file.exists(cfg$johansen),file.exists(cfg$directional))
before <- tools::md5sum(c(cfg$johansen,cfg$directional))
dir.create(cfg$out,recursive=TRUE,showWarnings=FALSE)
j <- read.csv(cfg$johansen)
direction <- read.csv(cfg$directional)
names(j)[names(j)=='EV_resolution'] <- 'EV_resolution_kyr'
names(j)[names(j)=='Trace_r0'] <- 'trace_r0'
names(j)[names(j)=='Trace_r_le1'] <- 'trace_r_le1'
direction <- data.frame(Taxon=direction$Taxon,Biological_variable='Extinction',
  EV_resolution_kyr=ifelse(direction$Taxon=='Elasmobranchs',50,100),
  test=ifelse(direction$Test=='short_run','Short-run cross-lag','Joint ECT + cross-lag'),
  equation=ifelse(direction$Direction=='EV to Extinction','Extinction change','Global EV change'),
  F_statistic=direction$F,P=direction$P,df_den=direction$df2)
taxa <- c('Elasmobranchs','Bivalves')
rates <- c('Speciation','Extinction','Net diversification')
j <- j[j$Taxon%in%taxa & j$Biological_variable%in%rates,]
j <- j[order(match(j$Biological_variable,rates),match(j$Taxon,taxa)),]
stopifnot(nrow(j)==6L,identical(as.character(j$Biological_variable),rep(rates,each=2L)),
          identical(as.character(j$Taxon),rep(taxa,3L)),
          identical(as.integer(j$EV_resolution_kyr),rep(c(50L,100L),3L)),
          all(is.finite(j$trace_r0)),all(is.finite(j$trace_r_le1)),
          identical(as.integer(j$rank_5pct),c(0L,2L,1L,1L,0L,0L)),
          identical(as.integer(j$rank_10pct),c(0L,2L,1L,1L,1L,0L)))
colors <- c(Elasmobranchs='#FAC3CF',Bivalves='#97CB79')
crit_colors <- c('1%'='#006400','5%'='#D62828','10%'='#173DDB')
theme_fig <- theme_classic(base_size=12.5,base_family='Helvetica')+
  theme(legend.position='none',axis.title=element_text(face='bold'),
        axis.text.x=element_text(face='bold',color='black',margin=margin(t=7)),
        axis.text.y=element_text(color='black'),
        plot.title=element_text(face='bold',hjust=.5,size=13.5,margin=margin(b=9)),
        axis.line=element_line(color='#666666',linewidth=.4),
        axis.ticks=element_line(color='#666666',linewidth=.4),
        plot.margin=margin(4,18,8,14))
trace_plot <- function(value,title,top,breaks,cv) {
  dd <- data.frame(rate=factor(j$Biological_variable,levels=rates),
                   taxon=factor(j$Taxon,levels=taxa),statistic=j[[value]])
  stopifnot(max(dd$statistic,cv)<top)
  ggplot(dd,aes(rate,statistic,fill=taxon))+
    geom_hline(yintercept=cv['1%'],linetype='dashed',color=crit_colors['1%'],linewidth=.55)+
    geom_hline(yintercept=cv['5%'],linetype='dashed',color=crit_colors['5%'],linewidth=.55)+
    geom_hline(yintercept=cv['10%'],linetype='dashed',color=crit_colors['10%'],linewidth=.55)+
    geom_col(position=position_dodge(width=.72),width=.32,color='#404040',linewidth=.35)+
    scale_fill_manual(values=colors,drop=FALSE)+
    scale_x_discrete(labels=c('Speciation','Extinction','Net diversif.'),drop=FALSE)+
    scale_y_continuous(limits=c(0,top),breaks=breaks,expand=expansion(mult=c(0,0)))+
    labs(x=NULL,y='Trace statistic',title=title)+theme_fig
}
p0 <- trace_plot('trace_r0','No cointegration (r = 0)',35,seq(0,35,5),
                 c('1%'=23.52,'5%'=17.95,'10%'=15.66))
p1 <- trace_plot('trace_r_le1','At most one relationship (r ≤ 1)',12,seq(0,12,3),
                 c('1%'=11.65,'5%'=8.18,'10%'=6.5))

# Only the two reliable extinction rank-1 models are eligible for B.
b <- direction[direction$Taxon%in%taxa & direction$Biological_variable=='Extinction' &
                 direction$test=='Short-run cross-lag',]
b <- b[order(match(b$Taxon,taxa),match(b$equation,c('Extinction change','Global EV change'))),]
stopifnot(nrow(b)==4L,identical(as.character(b$Taxon),rep(taxa,each=2L)),
          identical(as.character(b$equation),rep(c('Extinction change','Global EV change'),2L)),
          identical(as.integer(b$EV_resolution_kyr),c(50L,50L,100L,100L)),
          identical(as.integer(b$df_den),c(42L,42L,18L,18L)),
          b$P[1]<.05,b$P[2]>=.05,b$P[3]<.05,b$P[4]>=.05)
b$direction <- factor(rep(c('ΔEV(t−1) →\nΔExtinction(t)',
                            'ΔExtinction(t−1) →\nΔEV(t)'),2L),
                      levels=c('ΔEV(t−1) →\nΔExtinction(t)',
                               'ΔExtinction(t−1) →\nΔEV(t)'))
b$label <- vapply(b$P,function(p) sprintf('P = %s',formatC(p,format='f',digits=if(p<.01)4L else 3L)),character(1))
direction_plot <- function(taxon) {
  dd <- b[b$Taxon==taxon,]
  c0 <- unname(colors[taxon])
  shades <- setNames(c(c0,grDevices::adjustcolor(c0,alpha.f=.40)),levels(b$direction))
  # Neutral black outline and nominal P labels guide inference; F axes are
  # visual summaries of fixed-model tests, not directly comparable effect sizes.
  ggplot(dd,aes(direction,F_statistic,fill=direction))+
    geom_col(width=.48,color='#404040',linewidth=.35)+
    geom_text(aes(label=label),vjust=-.5,size=3.65,family='Helvetica')+
    scale_fill_manual(values=shades,drop=FALSE)+
    scale_y_continuous(limits=c(0,17.2),breaks=seq(0,15,5),
                       expand=expansion(mult=c(0,0)))+
    labs(x=NULL,y='F statistic',title=sprintf('%s extinction (%d kyr EV)',
       taxon,if(taxon=='Elasmobranchs')50L else 100L))+theme_fig+
    theme(axis.text.x=element_text(size=10.5,face='bold'),
          plot.title=element_text(size=13.5,face='bold',hjust=.5))
}
pb_e <- direction_plot(taxa[1]);pb_b <- direction_plot(taxa[2])
txt <- function(label,x,y,size=12,face='plain',just='center',col='#171717')
  grid.text(label,x=x,y=y,just=just,gp=gpar(fontfamily='Helvetica',fontface=face,
                                           fontsize=size,col=col))
draw <- function() {
  grid.newpage()
  txt('Global EV and marine turnover across the K–Pg transition',.5,.976,17,'bold')
  txt('A',.018,.923,21,'bold',just='left')
  txt('Johansen trace tests across biological rates',.5,.923,15,'bold')
  g0 <- ggplotGrob(p0);g1 <- ggplotGrob(p1)
  h <- unit.pmax(g0$heights,g1$heights);g0$heights <- h;g1$heights <- h
  pushViewport(viewport(x=.26,y=.727,width=.48,height=.335));grid.draw(g0);popViewport()
  pushViewport(viewport(x=.75,y=.727,width=.48,height=.335));grid.draw(g1);popViewport()
  txt('Group',.07,.543,10.8,'bold',just='left')
  grid.rect(x=.155,y=.543,width=.025,height=.018,gp=gpar(fill=colors[1],col='#404040'))
  txt('Elasmobranchs (50 kyr)',.175,.543,10.8,'bold',just='left')
  grid.rect(x=.42,y=.543,width=.025,height=.018,gp=gpar(fill=colors[2],col='#404040'))
  txt('Bivalves (100 kyr)',.44,.543,10.8,'bold',just='left')
  txt('Critical value',.61,.543,10.8,'bold',just='left')
  for(i in seq_along(crit_colors)) {
    x <- c(.725,.825,.925)[i]
    grid.lines(x=c(x,x+.028),y=c(.543,.543),gp=gpar(col=crit_colors[i],lty=2,lwd=1.5))
    txt(names(crit_colors)[i],x+.034,.543,10.8,just='left')
  }
  grid.lines(x=c(.035,.965),y=c(.515,.515),gp=gpar(col='#B6B6B6',lwd=.7))
  txt('B',.018,.487,21,'bold',just='left')
  txt('One-bin lagged prediction: extinction models only',.5,.487,15,'bold')
  g2 <- ggplotGrob(pb_e);g3 <- ggplotGrob(pb_b)
  h2 <- unit.pmax(g2$heights,g3$heights);g2$heights <- h2;g3$heights <- h2
  pushViewport(viewport(x=.26,y=.265,width=.48,height=.385));grid.draw(g2);popViewport()
  pushViewport(viewport(x=.75,y=.265,width=.48,height=.385));grid.draw(g3);popViewport()
  txt('Left bar: lagged ΔEV predicts subsequent ΔExtinction; right bar: reverse direction',
      .5,.048,10.3)
  txt('Nominal P values; predictive direction does not establish physical causation',
      .5,.023,9.6)
}
stem <- file.path(cfg$out,'Fig4_revised_all_rates_bar_direction')
if(capabilities('aqua')) {
  quartz(type='pdf',file=paste0(stem,'.pdf'),width=12,height=8.6,family='Helvetica')
} else {
  pdf(file=paste0(stem,'.pdf'),width=12,height=8.6,family='Helvetica',useDingbats=FALSE)
}
draw();dev.off()
ragg::agg_png(paste0(stem,'.png'),width=12,height=8.6,units='in',res=300,background='white')
draw();dev.off()
svglite::svglite(paste0(stem,'.svg'),width=12,height=8.6,bg='white')
draw();dev.off()
stopifnot(identical(before,tools::md5sum(c(cfg$johansen,cfg$directional))))
print(j[,c('Taxon','Biological_variable','trace_r0','trace_r_le1','rank_5pct','rank_10pct')],row.names=FALSE)
print(b[,c('Taxon','direction','F_statistic','df_den','P')],row.names=FALSE)
