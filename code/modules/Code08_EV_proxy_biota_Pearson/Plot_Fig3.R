# Expanded final Full-transition plot from plot_split_color07.R / plot_split_no_r.R.
# Reads the 45 recomputed Code08 full-transition results; no new statistical tests.
# Original submitted correlation figure is stored as Fig. 3 var._biota_pearson.
# Its caption defines three circle sizes: P<.005, P<.05, and nonsignificant.
library(ggplot2)
library(patchwork)
args <- gsub("~+~", " ", commandArgs(trailingOnly=FALSE), fixed=TRUE)
root <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
src <- c(file.path(root,'results/current_dataset/Full_transition_Pearson_45.csv'),
         file.path(root,'results/current_dataset/pCO2_variability_Pearson_9.csv'))
out <- file.path(root,'results/generated/figures_from_frozen_snapshot')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
h <- tools::md5sum(src)
d <- read.csv(src[1],fileEncoding='UTF-8-BOM',check.names=FALSE)
co2 <- read.csv(src[2]);co2$Stage <- 'Full_transition'
columns <- c('Taxon','Biological_variable','Stage','Environmental_metric','Pearson_r','Nominal_P','n')
d <- rbind(d[,columns],co2[,columns])
stopifnot(nrow(d)==54,!anyDuplicated(d[c('Taxon','Biological_variable','Stage','Environmental_metric')]),
          all(is.finite(d$Pearson_r)),all(is.finite(d$Nominal_P)))
d$Taxon <- factor(d$Taxon,levels=c('Elasmobranchs','Bivalves','Calcareous nannofossils'))
d$rate <- factor(d$Biological_variable,levels=c('Net diversification','Extinction','Speciation'))
d$metric <- factor(d$Environmental_metric,levels=c('Global EV','Temperature variability','d13C variability','d18O variability','Weathering variability','pCO2 variability'))
lev <- c('P >= 0.05','0.005 <= P < 0.05','P < 0.005')
d$significance <- factor(ifelse(d$Nominal_P<.005,lev[3],ifelse(d$Nominal_P<.05,lev[2],lev[1])),levels=lev)
d$label <- sprintf('%.2f',d$Pearson_r)
d$text_color <- ifelse(abs(d$Pearson_r)>=.65 & d$significance!=lev[1],'white','black')
stages <- c('Full_transition','Pre_boundary','Boundary_early_recovery','Post_recovery_background')
titles <- c('A. Full transition','B. Pre-boundary (>66.0 Ma)',
            'C. Boundary / early recovery (66.0-65.5 Ma)','D. Post-recovery background (<65.5 Ma)')
make_panel <- function(i){
 a <- d[d$Stage==stages[i],];stopifnot(nrow(a)==54)
 ggplot(a,aes(metric,rate))+
  geom_tile(fill='white',color='#B8B8B8',linewidth=.5)+
  geom_point(aes(size=significance,fill=Pearson_r),shape=21,color='white',stroke=0,show.legend=TRUE)+
  facet_grid(Taxon~.,switch='y')+
  coord_fixed(ratio=1)+
  scale_color_identity()+
  scale_size_manual(values=setNames(c(7.5,12,17),lev),limits=lev,breaks=lev,drop=FALSE,
       name='Nominal P',labels=expression(P>=0.05,0.005<=P~'<'~0.05,P<0.005))+
  scale_fill_gradient2(low='#2C8AC4',mid='white',high='#E63325',midpoint=0,
       limits=c(-.7,.7),oob=scales::squish,breaks=c(-.7,-.6,-.4,-.2,0,.2,.4,.6,.7),name='Strength of the correlation')+
  scale_x_discrete(expand=expansion(add=.5),labels=c('Global EV','Temp. var.',
       expression(bold(delta^{13}*C~var.)),expression(bold(delta^{18}*O~var.)),'Weath. var.','pCO2 var.'))+
  scale_y_discrete(expand=expansion(add=.5),labels=c('Net diversif.','Extinction','Speciation'))+
  labs(x=NULL,y=NULL,title=titles[i])+
  guides(fill=guide_colorbar(order=1,barheight=unit(8,'cm'),barwidth=unit(.45,'cm'),title.position='right',
            theme=theme(legend.title=element_text(angle=90,face='bold',size=10))),
         size=guide_legend(order=2,override.aes=list(fill='#999999',color='white',stroke=0),title.position='top'))+
  theme_bw(base_size=11,base_family='Helvetica')+
  theme(panel.grid=element_blank(),panel.border=element_rect(color='black',fill=NA,linewidth=.65),
    axis.ticks=element_blank(),axis.text=element_text(color='black'),
    axis.text.x=element_text(face='bold',angle=55,hjust=1,vjust=1,size=10),
    axis.text.y=element_text(size=10),strip.placement='outside',strip.background=element_blank(),
    strip.text.y.left=element_text(face='bold',angle=90,size=10,margin=margin(r=10)),
    panel.spacing.y=unit(0,'pt'),plot.title=element_text(face='bold',size=12,margin=margin(b=10)),
    plot.margin=margin(10,12,12,6),legend.position='right',legend.text=element_text(size=9),
    legend.title=element_text(face='bold',size=10),legend.key.height=unit(1.8,'cm'))
}
p <- make_panel(1)+labs(title='Full transition')
stem <- file.path(out,'Fig3_Full_transition_no_r_color07')
pdf(paste0(stem,'.pdf'),width=7.5,height=9,family='Helvetica',useDingbats=FALSE);print(p);dev.off()
ragg::agg_png(paste0(stem,'.png'),width=7.5,height=9,units='in',res=300,background='white');print(p);dev.off()
svglite::svglite(paste0(stem,'.svg'),width=7.5,height=9,bg='white');print(p);dev.off()
stopifnot(identical(h,tools::md5sum(src)))
writeLines(c('All 54 current full-transition r/P/n results reused; source MD5 unchanged.',paste(names(h),h),
 'Three original circle-size classes: P>=0.05; 0.005<=P<0.05; P<0.005. No fourth class or stars.',
 'Original order: Speciation, Extinction, Net diversification; Global EV, temperature, d13C, d18O, weathering.',
 'No r labels. Color limits -0.7 to +0.7. Original plot geometry, three P classes and order retained.'),file.path(out,'Plot_only_validation.txt'))
