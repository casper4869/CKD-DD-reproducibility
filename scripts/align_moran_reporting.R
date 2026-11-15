suppressPackageStartupMessages({library(spdep);library(dplyr)})
source('scripts/config.R')
a<-core_dir;spdir<-spatial_dir
d<-read.csv(file.path(a,'S13_GWR_diagnostics.csv'));local<-read.csv(file.path(a,'S12_GWR_coefficients.csv'))
canon<-read.csv(file.path(spdir,'GWR_sensitivity_diagnostics.csv'))
for(i in which(d$spec=='primary')){row<-canon[canon$coordinate=='original_city_points'&canon$spec=='SDI_PM25'&canon$target==d$target[i],];d$Moran_p[i]<-row$Moran_p_plus1_two_sided}
set.seed(2026110411)
for(i in which(d$spec=='population')){
 z<-local[local$spec=='population'&local$target==d$target[i],];coords<-as.matrix(z[,c('lng','lat')]);nb<-knn2nb(knearneigh(coords,k=4,longlat=TRUE),sym=TRUE);lw<-nb2listw(nb,style='W')
 mc<-moran.mc(z$residual,lw,nsim=999,alternative='two.sided');observed<-tail(mc$res,1);simulated<-head(mc$res,-1)
 d$Moran_p[i]<-min(1,2*min((1+sum(simulated>=observed))/1000,(1+sum(simulated<=observed))/1000))
}
write.csv(d,file.path(a,'S13_GWR_diagnostics.csv'),row.names=FALSE)
cat('Moran plus-one probabilities aligned\n')
