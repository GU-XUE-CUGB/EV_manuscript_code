# Plot-only: reference geometry from the author's 1611 x 1671, 400-dpi Fig. 3 TIFF.
# All r/P/n are read unchanged from the completed mean-state correlation table.
library(grid)
args <- gsub("~+~", " ", commandArgs(trailingOnly=FALSE), fixed=TRUE)
root <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
src<-file.path(root,'results/Mean_state_biota_Pearson_results.csv')
out<-file.path(root,'results/figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(src,check.names=FALSE);h<-tools::md5sum(src)
taxa<-c('Elasmobranchs','Bivalves','Calcareous nannofossils')
rates<-c('Speciation','Extinction','Net diversification')
metrics<-c('Global mean environmental state','Temperature mean state','d13C mean state','d18O mean state','Weathering mean state')
stopifnot(nrow(d)==45,!anyDuplicated(d[c('Taxon','Biological_variable','Mean_environmental_metric')]))
# Keep all template geometry and type sizes; add only outer whitespace so the
# longer mean-state label and full middle P legend are not clipped.
W<-1651;H<-1716;dpi<-400
left<-333;right<-1141;top<-2;bottom<-1455
cw<-(right-left)/5;ch<-(bottom-top)/9
palette<-colorRamp(c('#2C8AC4','#FFFFFF','#E63325'),space='Lab')
color<-function(r){rgb(palette(pmax(0,pmin(1,(r+.7)/1.4))),maxColorValue=255)}
px<-function(x)unit(x/W,'npc');py<-function(y)unit(1-y/H,'npc')
line<-function(x,y,col='black',lwd=.45){grid.lines(px(x),py(y),gp=gpar(col=col,lwd=lwd))}
txt<-function(label,x,y,size=7.5,bold=FALSE,rot=0,just='centre')grid.text(label,px(x),py(y),rot=rot,just=just,gp=gpar(fontfamily='Helvetica',fontsize=size,fontface=if(bold)'bold' else 'plain'))
circle<-function(x,y,r,fill,border='#888888')grid.circle(px(x),py(y),r=unit(r/dpi,'in'),gp=gpar(fill=fill,col=border,lwd=.3))
rect<-function(x,y,w,hh,fill,col=NA,lwd=.4)grid.rect(px(x+w/2),py(y+hh/2),width=unit(w/W,'npc'),height=unit(hh/H,'npc'),gp=gpar(fill=fill,col=col,lwd=lwd))
curve<-function(p){t<-seq(0,1,length.out=25);z<-outer((1-t)^3,p[1,])+outer(3*(1-t)^2*t,p[2,])+outer(3*(1-t)*t^2,p[3,])+outer(t^3,p[4,]);line(z[,1],z[,2],lwd=.45)}
brace<-function(y1,y2){mid<-(y1+y2)/2
 curve(rbind(c(91,y1),c(82,y1),c(77,y1+5),c(77,y1+12)))
 line(c(77,77),c(y1+12,mid-8));line(c(77,61,77),c(mid-8,mid,mid+8))
 line(c(77,77),c(mid+8,y2-12))
 curve(rbind(c(77,y2-12),c(77,y2-5),c(82,y2),c(91,y2)))
}
draw<-function(){
 grid.newpage();rect(0,0,W,H,'white')
 for(j in 1:4)line(rep(left+j*cw,2),c(top,bottom),'#B8B8B8',.55)
 for(i in 1:8)line(c(left,right),rep(top+i*ch,2),'#B8B8B8',.55)
 rect(left,top,right-left,bottom-top,NA,'black',.65)
 for(t in 1:3){
  yy<-top+((t-1)*3+c(.5,2.5))*ch;brace(yy[1],yy[2])
  txt(c('Elasmobranchs','Bivalves','Nannofossils')[t],19,mean(yy),8,TRUE,90)
  for(k in 1:3){y<-top+((t-1)*3+k-.5)*ch;txt(c('Speciation','Extinction','Net diversif.')[k],108,y,7.5,just='left')
   for(j in 1:5){a<-d[d$Taxon==taxa[t]&d$Biological_variable==rates[k]&d$Mean_environmental_metric==metrics[j],];stopifnot(nrow(a)==1)
    radius<-if(a$Nominal_P<.005)65 else if(a$Nominal_P<.05)51.5 else 27.5
    circle(left+(j-.5)*cw,y,radius,color(a$Pearson_r))
   }
  }
 }
 labels<-list('Global mean','Temp.',expression(bold(delta^{13}*C)),expression(bold(delta^{18}*O)),'Weath.')
 for(j in 1:5)txt(labels[[j]],left+(j-.5)*cw+25,bottom+22,8,TRUE,55,just='right')
 # Original Fig. 3 colorbar spans -0.7 to +0.7; labels every 0.2, minor ticks every 0.1.
 x<-1250;y<-125;hh<-800;ww<-51;steps<-512
 for(i in 0:(steps-1))rect(x,y+i*hh/steps,ww,hh/steps+.05,color(.7-1.4*(i+.5)/steps))
 rect(x,y,ww,hh,NA,'black',.65)
 for(v in seq(-.7,.7,.1)){yy<-y+(.7-v)/1.4*hh;line(c(x+ww,x+ww+10),c(yy,yy),'black',.4)}
 for(v in seq(-.6,.6,.2))txt(sprintf('%.1f',if(abs(v)<1e-9)0 else v),x+ww+23,y+(.7-v)/1.4*hh,6.5,just='left')
 txt('Strength of the correlation',1433,y+hh/2,8,rot=90)
 for(i in 1:3){yy<-c(1054,1190,1346)[i];circle(1276,yy,c(27.5,51.5,65)[i],'#D6D6D6',NA)
  txt(list(expression(P>=0.05),expression(0.005<=P~'<'~0.05),expression(P<0.005))[[i]],1368,yy,6.5,just='left')}
}
stem<-file.path(out,'Fig_S13_mean_biota_Pearson_revised_Fig3_style')
pdf(file=paste0(stem,'.pdf'),width=W/dpi,height=H/dpi,family='Helvetica',useDingbats=FALSE);draw();dev.off()
ragg::agg_png(paste0(stem,'.png'),width=W,height=H,units='px',res=dpi,background='white');draw();dev.off()
svglite::svglite(paste0(stem,'.svg'),width=W/dpi,height=H/dpi,bg='white');draw();dev.off()
stopifnot(identical(h,tools::md5sum(src)))
writeLines(c('Plot-only; all 45 mean-state r/P/n values unchanged.',
 'Color mapping: -0.7 to +0.7, zero white, negative blue, positive red. Values beyond limits use endpoint colors only.',
 'Circle-size classes: P>=0.05; 0.005<=P<0.05; P<0.005. No r overlays.',
 'Internal geometry follows the provided 1611x1671 Fig. 3 TIFF at 400 dpi; 40 px right/45 px bottom outer padding avoids clipping. Square cells, gray circle outlines, taxon braces, side colorbar and three-size legend retained.',
 'Five existing mean-state columns retained; pCO2 is not added to the figure.',
 paste('Unchanged source MD5:',unname(h))),file.path(out,'Plot_verification.txt'))
