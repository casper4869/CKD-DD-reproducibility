options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(vars);library(dplyr)})
source('scripts/config.R');out<-temporal_dir
dat<-read.csv(file.path(core_dir,'S1_analysis_panel.csv'));dat<-dat[dat$sex_name=='Both',]
eq<-function(fit,target,source){X<-as.matrix(fit$datamat[,-seq_len(fit$K),drop=FALSE]);m<-fit$varresult[[target]];j<-match(paste0(source,'.l',seq_len(fit$p)),colnames(X));r<-lm.fit(X[,-j,drop=FALSE],fit$datamat[[target]]);F<-max(0,((sum(r$residuals^2)-sum(resid(m)^2))/length(j))/(sum(resid(m)^2)/df.residual(m)));pf(F,length(j),df.residual(m),lower.tail=FALSE)}
sp<-function(fit,pt=4,bg=1,typ='PT.adjusted')tryCatch(serial.test(fit,lags.pt=pt,lags.bg=bg,type=typ)$serial$p.value,error=function(e)NA_real_)
rr<-list();k<-0;pre<-list();leverage<-list()
for(d in split(dat,dat$location_id)){
 d<-d[order(d$year),];raw<-as.matrix(d[,c('CKD','DD','SDI')]);lg<-raw;lg[,1:2]<-log(lg[,1:2])
 fits<-list(original_raw_endo=VAR(raw,p=2,type='const'),TY_log_endo_p2_d1_const=VAR(lg,p=3,type='const'),TY_log_endo_p2_d1_both=VAR(lg,p=3,type='both'),TY_log_endo_p2_d2_const=VAR(lg,p=4,type='const'),TY_log_endo_p2_d2_both=VAR(lg,p=4,type='both'))
 for(nm in names(fits)){
  fit<-fits[[nm]];u<-resid(fit)
  k<-k+1;rr[[k]]<-data.frame(location_id=d$location_id[1],model=nm,PT4=if(fit$p<4)sp(fit,4) else NA,PT8=sp(fit,8),ES1=sp(fit,bg=1,typ='ES'),ES2=sp(fit,bg=2,typ='ES'),
   residual_lag1_CKD=cor(u[-1,1],u[-nrow(u),1]),residual_lag1_DD=cor(u[-1,2],u[-nrow(u),2]),max_root=max(roots(fit)))
 }
 fit<-fits[[1]];X<-as.matrix(fit$datamat[,-1:-3,drop=FALSE]);qq<-qr.Q(qr(X));h<-rowSums(qq*qq)
 leverage[[length(leverage)+1]]<-data.frame(location_id=d$location_id[1],location_name=d$location_name[1],year=tail(d$year,length(h)),leverage=h)
 for(ending in c(2019,2021)){
  dd<-d[d$year<=ending,];Y<-as.matrix(dd[,c('CKD','DD','SDI')]);p<-as.integer(VARselect(Y,lag.max=2,type='const')$selection['AIC(n)']);fit<-VAR(Y,p=p,type='const')
  pre[[length(pre)+1]]<-data.frame(location_id=d$location_id[1],location_name=d$location_name[1],ending=ending,p=p,p_CKD_DD=eq(fit,'DD','CKD'),p_DD_CKD=eq(fit,'CKD','DD'),serial_PT4=sp(fit,4),serial_ES1=sp(fit,bg=1,typ='ES'),max_root=max(roots(fit)))
 }
}
diag<-bind_rows(rr);write.csv(diag,file.path(out,'A17_residual_diagnostics_country.csv'),row.names=FALSE)
write.csv(diag%>%group_by(model)%>%summarise(n=n(),PT4_reject=sum(PT4<.05,na.rm=TRUE),PT8_reject=sum(PT8<.05,na.rm=TRUE),ES1_reject=sum(ES1<.05,na.rm=TRUE),ES2_reject=sum(ES2<.05,na.rm=TRUE),median_abs_residual_lag1_CKD=median(abs(residual_lag1_CKD)),median_abs_residual_lag1_DD=median(abs(residual_lag1_DD)),stable=sum(max_root<1),.groups='drop'),file.path(out,'A18_residual_diagnostics_summary.csv'),row.names=FALSE)
pre<-bind_rows(pre)
pre<-bind_rows(lapply(split(pre,pre$ending),function(z){q<-p.adjust(c(z$p_CKD_DD,z$p_DD_CKD),'BH');z$q_CKD_DD<-q[1:nrow(z)];z$q_DD_CKD<-tail(q,nrow(z));z}))
write.csv(pre,file.path(out,'A19_matched_model_pre_pandemic.csv'),row.names=FALSE)
write.csv(pre%>%group_by(ending)%>%summarise(n=n(),CKD_DD=sum(q_CKD_DD<.05),DD_CKD=sum(q_DD_CKD<.05),both=sum(q_CKD_DD<.05&q_DD_CKD<.05),neither=sum(q_CKD_DD>=.05&q_DD_CKD>=.05),p2=sum(p==2),serial_PT4_reject=sum(serial_PT4<.05,na.rm=TRUE),serial_ES1_reject=sum(serial_ES1<.05,na.rm=TRUE),stable=sum(max_root<1),.groups='drop'),file.path(out,'A20_matched_model_pre_pandemic_summary.csv'),row.names=FALSE)
lev<-bind_rows(leverage);write.csv(lev,file.path(out,'A21_original_leverage.csv'),row.names=FALSE)
write.csv(lev%>%group_by(year)%>%summarise(n=n(),median=median(leverage),min=min(leverage),max=max(leverage),over_point5=sum(leverage>.5),over_point8=sum(leverage>.8),over_point9=sum(leverage>.9),.groups='drop'),file.path(out,'A22_original_leverage_summary.csv'),row.names=FALSE)
