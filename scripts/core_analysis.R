options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(vars);library(forecast);library(readr);library(jsonlite)})
select<-dplyr::select
source('scripts/config.R')
inp<-private; out<-core_dir
set.seed(20261104)
savecsv<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
ckd<-read.csv(file.path(inp,'CKD_ASIR_1992_2021.csv')); dd<-read.csv(file.path(inp,'DD_ASIR_1992_2021.csv'))
stopifnot(nrow(ckd)==18360,nrow(dd)==18360)
dat<-inner_join(ckd%>%select(location_id,location_name,sex_name,year,CKD=val), dd%>%select(location_id,sex_name,year,DD=val),by=c('location_id','sex_name','year'))
sdi<-read.csv(file.path(inp,'SDI.csv'))%>%pivot_longer(starts_with('val_'),names_prefix='val_',names_to='year',values_to='SDI')%>%mutate(year=as.integer(year))
dat<-left_join(dat,sdi,by=c('location_name','year'))
stopifnot(nrow(dat)==18360,!anyNA(dat),length(unique(dat$location_id))==204)
savecsv(dat,'S1_analysis_panel')
cat('Input matched: 204 locations, 30 years, 3 sex strata\n')
# Quartile display uses rank with location_id tie break; high-high is Q4/Q4 only.
qdat<-dat%>%arrange(location_id)%>%group_by(sex_name,year)%>%mutate(q_CKD=ntile(CKD,4),q_DD=ntile(DD,4),type=case_when(q_CKD==q_DD~'Same quartile',q_CKD>q_DD~'Higher CKD quartile',TRUE~'Higher DD quartile'),high_high=q_CKD==4&q_DD==4)%>%ungroup()
savecsv(qdat,'S2_quartile_country_year')
qsum<-qdat%>%group_by(sex_name,year,type)%>%summarise(n=n(),.groups='drop'); savecsv(qsum,'S3_quartile_counts')
thresholds<-bind_rows(lapply(c('Both','Female','Male'),function(s){d<-filter(dat,sex_name==s,year==2021); bind_rows(lapply(c(.70,.75,.80),function(th){n<-nrow(d); m<-ceiling((1-th)*n); a<-rank(-d$CKD,ties.method='first')<=m;b<-rank(-d$DD,ties.method='first')<=m;k<-sum(a&b); cor<-cor.test(d$CKD,d$DD,method='spearman',exact=FALSE);data.frame(sex=s,top_fraction=m/n,threshold_percentile=100*th,n=n,top_n=m,observed=k,expected=m*m/n,enrichment=k/(m*m/n),upper_tail_p=phyper(k-1,m,n-m,m,lower.tail=FALSE),lower_tail_p=phyper(k,m,n-m,m),rho=unname(cor$estimate),rho_p=cor$p.value)}))}))
thresholds$q_upper_BH<-p.adjust(thresholds$upper_tail_p,'BH');savecsv(thresholds,'S4_threshold_sensitivity')
cat('Quartile sensitivity complete\n')
# Descriptive within-country correlation across 19 age groups from 5-9 to 95+ years.
ac<-read.csv(file.path(inp,'CKD_age_2021.csv'));ad<-read.csv(file.path(inp,'DD_age_2021.csv'))
ages<-c('5-9 years','10-14 years',paste0(seq(15,90,5),'-',seq(19,94,5),' years'),'95+ years')
ag<-inner_join(ac%>%select(location_id,location_name,sex_name,age_name,CKD=val),ad%>%select(location_id,sex_name,age_name,DD=val),by=c('location_id','sex_name','age_name'))%>%filter(age_name%in%ages)
stopifnot(nrow(ag)==204*2*19)
savecsv(ag,'S5_age_source')
profile<-ag%>%group_by(location_id,location_name,sex_name)%>%summarise(n_age=n(),profile_r=cor(CKD,DD),.groups='drop');savecsv(profile,'S6_age_profile_descriptive')
agecors<-ag%>%group_by(sex_name,age_name)%>%summarise(n_country=n(),rho=cor(CKD,DD,method='spearman'),p=cor.test(CKD,DD,method='spearman',exact=FALSE)$p.value,.groups='drop')%>%group_by(sex_name)%>%mutate(q_BH=p.adjust(p,'BH'))%>%ungroup()%>%mutate(lower=NA_real_,upper=NA_real_)
# Paired country bootstrap preserves the covariance among age-band correlations.
B<-2000; agediff<-list()
for(s in c('Female','Male')){
 z<-filter(ag,sex_name==s); cx<-pivot_wider(z%>%select(location_id,age_name,CKD),names_from=age_name,values_from=CKD)%>%arrange(location_id);dy<-pivot_wider(z%>%select(location_id,age_name,DD),names_from=age_name,values_from=DD)%>%arrange(location_id)
 matx<-as.matrix(cx[,ages]);maty<-as.matrix(dy[,ages]);lo<-which(ages%in%paste0(seq(15,45,5),'-',seq(19,49,5),' years'));hi<-which(ages%in%c(paste0(seq(50,90,5),'-',seq(54,94,5),' years'),'95+ years'))
 stat<-function(ii){rr<-sapply(seq_along(ages),function(j)cor(matx[ii,j],maty[ii,j],method='spearman'));c(rr,mean(rr[hi])-mean(rr[lo]))}
 obs<-stat(seq_len(204)); bt<-replicate(B,stat(sample.int(204,204,TRUE)))
 ci<-t(apply(bt,1,quantile,c(.025,.975),na.rm=TRUE));agecors$lower[agecors$sex_name==s]<-ci[match(agecors$age_name[agecors$sex_name==s],ages),1];agecors$upper[agecors$sex_name==s]<-ci[match(agecors$age_name[agecors$sex_name==s],ages),2]
 agediff[[s]]<-data.frame(sex=s,difference_mean_rho=tail(obs,1),lower=tail(ci[,1],1),upper=tail(ci[,2],1),definition='Mean age-specific country Spearman rho: 50+ minus 15-49 years')
}
savecsv(agecors,'S7_age_correlations');savecsv(bind_rows(agediff),'S8_age_contrast')
cat('Age analysis complete\n')
# Differenced/log-level sensitivity only. Final primary raw-level, endogenous-SDI tests are in temporal_models.R.
adfp<-function(z)tryCatch(suppressWarnings(tseries::adf.test(z,k=1)$p.value),error=function(e)NA_real_)
kpssp<-function(z)tryCatch(suppressWarnings(tseries::kpss.test(z,null='Level',lshort=TRUE)$p.value),error=function(e)NA_real_)
eqtest<-function(fit,target,source){m<-fit$varresult[[target]]; nm<-paste0(source,'.l',seq_len(fit$p)); X<-as.matrix(fit$datamat[,-seq_len(fit$K),drop=FALSE]);id<-match(nm,colnames(X));stopifnot(!anyNA(id),identical(colnames(X),names(coef(m))));y<-fit$datamat[[target]];r<-lm.fit(X[,-id,drop=FALSE],y);F<-((sum(r$residuals^2)-sum(resid(m)^2))/length(id))/(sum(resid(m)^2)/df.residual(m));c(F=F,p=pf(max(0,F),length(id),df.residual(m),lower.tail=FALSE))}
grun<-function(d,spec){
 d<-arrange(d,year);if(spec=='pre2019')d<-filter(d,year<=2019);Y<-cbind(CKD=log(d$CKD),DD=log(d$DD));X<-matrix(d$SDI,ncol=1,dimnames=list(NULL,'SDI'));typ<-'both'
 if(spec!='levels'){Y<-diff(Y);X<-diff(X);typ<-'const'}
 p<-if(spec=='lag1')1L else if(spec=='lag2')2L else as.integer(VARselect(Y,lag.max=2,type=typ,exogen=X)$selection['SC(n)'])
 fit<-VAR(Y,p=p,type=typ,exogen=X);a<-eqtest(fit,'DD','CKD');b<-eqtest(fit,'CKD','DD')
 white<-tryCatch(serial.test(fit,lags.pt=4,type='PT.adjusted')$serial$p.value,error=function(e)NA_real_)
 pa<-sapply(1:2,function(j)adfp(Y[,j]));pk<-sapply(1:2,function(j)kpssp(Y[,j]));stable<-max(roots(fit))<1
 data.frame(location_id=d$location_id[1],location_name=d$location_name[1],spec=spec,n=nrow(d),lag=p,p_CKD_DD=a['p'],p_DD_CKD=b['p'],ADF_CKD=pa[1],ADF_DD=pa[2],KPSS_CKD=pk[1],KPSS_DD=pk[2],ADF_level_CKD=adfp(log(d$CKD)),ADF_level_DD=adfp(log(d$DD)),max_root=max(roots(fit)),serial_p=white,stable=stable,diagnostics_pass=all(pa<.05,pk>=.05,stable,white>=.05,na.rm=FALSE))
}
both<-filter(dat,sex_name=='Both'); gl<-list()
for(spec in c('primary','lag1','lag2','pre2019','levels')){
 rr<-bind_rows(lapply(split(both,both$location_id),function(d)tryCatch(grun(d,spec),error=function(e)data.frame(location_id=d$location_id[1],location_name=d$location_name[1],spec=spec,error=conditionMessage(e)))))
 pp<-c(rr$p_CKD_DD,rr$p_DD_CKD);pp[!is.finite(pp)]<-1;adj<-p.adjust(pp,'BH');by<-p.adjust(pp,'BY');rr$q_CKD_DD<-head(adj,204);rr$q_DD_CKD<-tail(adj,204);rr$qBY_CKD_DD<-head(by,204);rr$qBY_DD_CKD<-tail(by,204)
 rr$category<-ifelse(rr$q_CKD_DD<.05&rr$q_DD_CKD<.05,'Both predictive tests',ifelse(rr$q_CKD_DD<.05,'CKD to DD only',ifelse(rr$q_DD_CKD<.05,'DD to CKD only','Neither test')))
 gl[[spec]]<-rr;cat('Predictive tests',spec,'complete\n')
}
gres<-bind_rows(gl);savecsv(gres,'S9_predictive_tests_diagnostics');gsum<-gres%>%group_by(spec,category)%>%summarise(n=n(),diagnostic_pass=sum(diagnostics_pass,na.rm=TRUE),.groups='drop');savecsv(gsum,'S10_predictive_summary')
# GWR: recreate the documented SDI + PM2.5 model. Population is a separate sensitivity.
coord<-read.csv(file.path(reference_dir,'Coordinates.csv'));pm<-read.csv(file.path(inp,'PM25.csv'));pop<-read.csv(file.path(inp,'Population.csv'))
geo<-filter(both,year==2021)%>%left_join(coord,by='location_name')%>%left_join(pm%>%select(location_name,PM25=val_2021),by='location_name')%>%left_join(pop%>%select(location_name,population=val_2021),by='location_name')%>%arrange(location_id)
stopifnot(nrow(geo)==204,!anyNA(geo));geo$log_population<-log(geo$population)
for(v in c('CKD','DD','SDI','PM25','log_population'))geo[[paste0('z',v)]]<-as.numeric(scale(geo[[v]]))
savecsv(geo,'S11_spatial_inputs'); gg<-list(); gd<-list()
for(target in c('CKD','DD')) for(spec in c('primary','population')){
 source<-if(target=='CKD')'DD' else 'CKD';form<-as.formula(paste0('z',target,' ~ z',source,' + zSDI + zPM25',if(spec=='population')' + zlog_population' else ''))
 bw<-spgwr::gwr.sel(form,data=geo,coords=as.matrix(geo[,c('lng','lat')]),longlat=TRUE,verbose=FALSE)
 fit<-spgwr::gwr(form,data=geo,coords=as.matrix(geo[,c('lng','lat')]),bandwidth=bw,longlat=TRUE,hatmatrix=TRUE,se.fit=TRUE)
 rr<-as.data.frame(fit$SDF);coefname<-paste0('z',source);secol<-if(paste0(coefname,'_se_EDF')%in%names(rr))paste0(coefname,'_se_EDF') else paste0(coefname,'_se');se<-rr[[secol]]
 edf<-fit$results$edf; pvals<-2*pt(-abs(rr[[coefname]]/se),df=edf)
 kn<-spdep::knearneigh(as.matrix(geo[,c('lng','lat')]),k=4,longlat=TRUE);nb<-spdep::knn2nb(kn,sym=TRUE);lw<-spdep::nb2listw(nb,style='W'); mt<-spdep::moran.mc(rr$gwr.e,lw,nsim=999,alternative='two.sided')
 design<-model.matrix(form,geo); coords<-as.matrix(geo[,c('lng','lat')]);cn<-sapply(seq_len(nrow(geo)),function(i){dist<-sp::spDistsN1(coords,coords[i,],longlat=TRUE);w<-spgwr::gwr.Gauss(dist^2,bw);kappa(design*sqrt(w),exact=TRUE)})
 result<-data.frame(location_id=geo$location_id,location_name=geo$location_name,target=target,spec=spec,coefficient=rr[[coefname]],se=se,p=pvals,q_BH=p.adjust(pvals,'BH'),local_R2=rr$localR2,condition_number=cn,residual=rr$gwr.e,lng=geo$lng,lat=geo$lat)
 gg[[paste(target,spec)]]<-result;gd[[paste(target,spec)]]<-data.frame(target=target,spec=spec,n=204,bandwidth_km=bw,edf=edf,R2=1-fit$results$rss/sum((geo[[paste0('z',target)]]-mean(geo[[paste0('z',target)]]))^2),AICc=fit$results$AICc,local_R2_min=min(rr$localR2),local_R2_median=median(rr$localR2),local_R2_max=max(rr$localR2),negative=sum(rr[[coefname]]<0),positive=sum(rr[[coefname]]>0),BH_local_significant=sum(p.adjust(pvals,'BH')<.05),condition_number_over30=sum(cn>30),Moran_I=unname(mt$statistic),Moran_p=mt$p.value)
 cat('GWR',target,spec,'complete\n')
}
savecsv(bind_rows(gg),'S12_GWR_coefficients');savecsv(bind_rows(gd),'S13_GWR_diagnostics')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
saveRDS(list(dat=dat,qdat=qdat,ag=ag,profile=profile,agecors=agecors,agediff=bind_rows(agediff),gres=gres,geo=geo,gwr=bind_rows(gg),gwrdiag=bind_rows(gd)),file.path(out,'analysis_objects.rds'))
cat('ANALYSIS CORE COMPLETE\n')
