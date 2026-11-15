options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(vars);library(dplyr);library(forecast);library(urca);library(sandwich)})
select<-dplyr::select
source('scripts/config.R')
out<-temporal_dir
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dat<-read.csv(file.path(core_dir,'S1_analysis_panel.csv'))
dat<-dat[dat$sex_name=='Both',];ds<-split(dat,dat$location_id)
ds<-lapply(ds,function(d)d[order(d$year),])
stopifnot(length(ds)==204,all(sapply(ds,nrow)==30))
write_out<-function(z,n)write.csv(z,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
eqtest<-function(fit,target,source,p_test=fit$p){
 m<-fit$varresult[[target]]
 X<-as.matrix(fit$datamat[,-seq_len(fit$K),drop=FALSE])
 stopifnot(identical(colnames(X),names(coef(m))))
 id<-match(paste0(source,'.l',seq_len(p_test)),colnames(X))
 stopifnot(!anyNA(id));y<-fit$datamat[[target]]
 r<-lm.fit(X[,-id,drop=FALSE],y)
 rss<-sum(resid(m)^2);df<-df.residual(m)
 Fv<-max(0,((sum(r$residuals^2)-rss)/length(id))/(rss/df))
 b<-coef(m)[id];V<-vcov(m)[id,id,drop=FALSE]
 W<-as.numeric(crossprod(b,solve(V,b)))
 # vars stores an lm formula whose model.matrix() adds an extra intercept;
 # construct HC3 from the verified, actual VAR design rather than that formula.
 q<-qr(X);Q<-qr.Q(q);R<-qr.R(q);Rinv<-backsolve(R,diag(ncol(R)))
 leverage<-rowSums(Q*Q)
 Hpivot<-Rinv%*%crossprod(Q*as.numeric(resid(m)/(1-leverage)))%*%t(Rinv)
 HC3full<-matrix(0,ncol(X),ncol(X));HC3full[q$pivot,q$pivot]<-Hpivot
 H<-HC3full[id,id,drop=FALSE]
 WH<-tryCatch(as.numeric(crossprod(b,solve(H,b))),error=function(e)NA_real_)
 stopifnot(abs(W-length(id)*Fv)<max(1e-5,abs(W)*1e-5))
 c(F=Fv,p=pf(Fv,length(id),df,lower.tail=FALSE),W=W,pW=pchisq(W,length(id),lower.tail=FALSE),WHC3=WH,pHC3=pchisq(WH,length(id),lower.tail=FALSE))
}
serialp<-function(fit,lags=8)tryCatch(serial.test(fit,lags.pt=lags,type='PT.adjusted')$serial$p.value,error=function(e)NA_real_)
tests<-list();k<-0
for(d in ds){
 raw<-as.matrix(d[,c('CKD','DD','SDI')]);lg<-raw;lg[,1:2]<-log(lg[,1:2])
 pold<-as.integer(VARselect(raw,lag.max=2,type='const')$selection['AIC(n)'])
 fits<-list(original_raw_endo=VAR(raw,p=pold,type='const'),
            raw_endo_trend=VAR(raw,p=pold,type='both'),
            log_endo_const=VAR(lg,p=pold,type='const'),
            log_endo_trend=VAR(lg,p=pold,type='both'))
 for(nm in names(fits)){
  fit<-fits[[nm]];a<-eqtest(fit,'DD','CKD');b<-eqtest(fit,'CKD','DD')
  k<-k+1;tests[[k]]<-data.frame(location_id=d$location_id[1],location_name=d$location_name[1],model=nm,p_base=pold,p_fit=fit$p,dmax=0,
    p_CKD_DD=a['p'],p_DD_CKD=b['p'],joint_CKD=causality(fit,cause='CKD')$Granger$p.value,joint_DD=causality(fit,cause='DD')$Granger$p.value,
    pW_CKD_DD=a['pW'],pW_DD_CKD=b['pW'],pHC3_CKD_DD=a['pHC3'],pHC3_DD_CKD=b['pHC3'],
    serial_p=serialp(fit),max_root=max(roots(fit)),df_residual=df.residual(fit$varresult[[1]]),
    design_condition=kappa(scale(as.matrix(fit$datamat[,-seq_len(fit$K),drop=FALSE]),center=FALSE,scale=apply(as.matrix(fit$datamat[,-seq_len(fit$K),drop=FALSE]),2,function(x)sqrt(sum(x*x))))))
 }
 # Toda-Yamamoto sensitivity: fixed base p=2 (or p=1 for a lag sensitivity), all three variables endogenous.
 # Test only first p lag coefficients; the dmax extra lag coefficients remain unrestricted.
 for(dmax in 1:2)for(pt in 1:2)for(typ in c('const','both')){
  fit<-VAR(lg,p=pt+dmax,type=typ)
  a<-eqtest(fit,'DD','CKD',pt);b<-eqtest(fit,'CKD','DD',pt)
  k<-k+1;tests[[k]]<-data.frame(location_id=d$location_id[1],location_name=d$location_name[1],
   model=paste0('TY_log_endo_p',pt,'_d',dmax,'_',typ),p_base=pt,p_fit=pt+dmax,dmax=dmax,
   p_CKD_DD=a['pW'],p_DD_CKD=b['pW'],joint_CKD=NA,joint_DD=NA,pW_CKD_DD=a['pW'],pW_DD_CKD=b['pW'],
   pHC3_CKD_DD=a['pHC3'],pHC3_DD_CKD=b['pHC3'],serial_p=serialp(fit),max_root=max(roots(fit)),
   df_residual=df.residual(fit$varresult[[1]]),design_condition=kappa(scale(as.matrix(fit$datamat[,-seq_len(fit$K),drop=FALSE]),center=FALSE,scale=apply(as.matrix(fit$datamat[,-seq_len(fit$K),drop=FALSE]),2,function(x)sqrt(sum(x*x))))))
 }
 if(length(tests)%%240==0)cat('Granger rows',length(tests),'\n')
}
tt<-bind_rows(tests)
adj_models<-function(z){
 n<-nrow(z)
 for(suf in c('','HC3_')){
  c1<-if(suf=='')'p_CKD_DD' else 'pHC3_CKD_DD';c2<-if(suf=='')'p_DD_CKD' else 'pHC3_DD_CKD'
  pp<-c(z[[c1]],z[[c2]]);pp[!is.finite(pp)]<-1;qq<-p.adjust(pp,'BH');by<-p.adjust(pp,'BY')
  z[[paste0('q_',suf,'CKD_DD')]]<-qq[1:n];z[[paste0('q_',suf,'DD_CKD')]]<-qq[n+1:n]
  z[[paste0('qBY_',suf,'CKD_DD')]]<-by[1:n];z[[paste0('qBY_',suf,'DD_CKD')]]<-by[n+1:n]
 }
 jointq<-p.adjust(c(z$joint_CKD,z$joint_DD),'BH')
 z$q_joint_CKD_pooled<-jointq[1:n];z$q_joint_DD_pooled<-tail(jointq,n)
 z$q_joint_CKD_separate<-p.adjust(z$joint_CKD,'BH');z$q_joint_DD_separate<-p.adjust(z$joint_DD,'BH')
 z$q_CKD_DD_separate<-p.adjust(z$p_CKD_DD,'BH');z$q_DD_CKD_separate<-p.adjust(z$p_DD_CKD,'BH')
 z
}
tt<-bind_rows(lapply(split(tt,tt$model),adj_models));write_out(tt,'A1_model_tests_country')
summ<-tt%>%group_by(model)%>%summarise(n=n(),p1=sum(p_base==1),p2=sum(p_base==2),
 joint_both_separate=sum(q_joint_CKD_separate<.05&q_joint_DD_separate<.05,na.rm=TRUE),
 joint_CKD_pooled=sum(q_joint_CKD_pooled<.05,na.rm=TRUE),joint_DD_pooled=sum(q_joint_DD_pooled<.05,na.rm=TRUE),joint_both_pooled=sum(q_joint_CKD_pooled<.05&q_joint_DD_pooled<.05,na.rm=TRUE),
 specific_both_separate=sum(q_CKD_DD_separate<.05&q_DD_CKD_separate<.05),
 BH_CKD_to_DD=sum(q_CKD_DD<.05),BH_DD_to_CKD=sum(q_DD_CKD<.05),BH_both=sum(q_CKD_DD<.05&q_DD_CKD<.05),
 BY_both=sum(qBY_CKD_DD<.05&qBY_DD_CKD<.05),HC3_BH_CKD_to_DD=sum(q_HC3_CKD_DD<.05),HC3_BH_DD_to_CKD=sum(q_HC3_DD_CKD<.05),
 HC3_BH_both=sum(q_HC3_CKD_DD<.05&q_HC3_DD_CKD<.05),serial_reject=sum(serial_p<.05,na.rm=TRUE),stable=sum(max_root<1),
 min_df=min(df_residual),median_condition=median(design_condition),.groups='drop')
write_out(summ,'A2_model_summary');print(summ,width=Inf)
# Focused unit-root audit: trend ADF in tseries was the previous all-or-none gate.
adfr<-function(x,type='drift',lags=1){a<-ur.df(x,type=type,lags=lags,selectlags='Fixed');c(stat=a@teststat[1,1],critical5=a@cval[1,'5pct'],reject=a@teststat[1,1]<a@cval[1,'5pct'])}
unit<-list();k<-0
for(d in ds)for(dis in c('CKD','DD'))for(ending in c(2019,2021))for(difference in 0:2)for(typ in c('drift','trend')){
 x<-log(d[[dis]][d$year<=ending]);if(difference>0)x<-diff(x,differences=difference)
 a<-adfr(x,typ,1);k<-k+1
 unit[[k]]<-data.frame(location_id=d$location_id[1],disease=dis,ending=ending,difference=difference,deterministic=typ,stat=a['stat'],critical5=a['critical5'],reject=a['reject'])
}
un<-bind_rows(unit);write_out(un,'A4_unit_root_audit')
write_out(un%>%group_by(disease,ending,difference,deterministic)%>%summarise(n=n(),reject=sum(reject),.groups='drop'),'A5_unit_root_summary')

# Exact original VARX equivalent: log1p, fixed p=2, constant, contemporaneous SDI.
# No data clipping in the main audit forecast; clipped original point results retained separately.
fcvar<-function(d,h,scenario='hold'){
 Y<-cbind(CKD=log1p(d$CKD),DD=log1p(d$DD));X<-matrix(d$SDI,ncol=1,dimnames=list(NULL,'SDI'))
 fit<-VAR(Y,p=2,type='const',exogen=X)
 fit$call$exogen<-X
 future<-if(scenario=='hold')rep(tail(d$SDI,1),h) else pmin(1,pmax(0,tail(d$SDI,1)+seq_len(h)*mean(tail(diff(d$SDI),5))))
 ff<-predict(fit,n.ahead=h,ci=.95,dumvar=matrix(future,ncol=1,dimnames=list(NULL,'SDI')))
 zz<-bind_rows(lapply(c('CKD','DD'),function(dis){x<-ff$fcst[[dis]];data.frame(disease=dis,h=1:h,median=expm1(x[,'fcst']),lower=expm1(x[,'lower']),upper=expm1(x[,'upper']),clipped=pmin(1.5*max(d[[dis]]),pmax(min(d[[dis]]),expm1(x[,'fcst']))),serial_p=serialp(fit,8),max_root=max(roots(fit)))}))
 zz
}
fcar<-function(d,h,scenario='hold'){
 future<-if(scenario=='hold')rep(tail(d$SDI,1),h) else pmin(1,pmax(0,tail(d$SDI,1)+seq_len(h)*mean(tail(diff(d$SDI),5))))
 bind_rows(lapply(c('CKD','DD'),function(dis){
  a<-auto.arima(log1p(d[[dis]]),xreg=matrix(d$SDI,ncol=1,dimnames=list(NULL,'SDI')),seasonal=FALSE,max.p=2,max.q=2,max.d=2,stepwise=TRUE,approximation=FALSE)
  f<-forecast(a,xreg=matrix(future,ncol=1,dimnames=list(NULL,'SDI')),h=h,level=95)
  data.frame(disease=dis,h=1:h,median=expm1(as.numeric(f$mean)),lower=expm1(as.numeric(f$lower)),upper=expm1(as.numeric(f$upper)),clipped=pmin(1.5*max(d[[dis]]),pmax(min(d[[dis]]),expm1(as.numeric(f$mean)))),serial_p=NA,max_root=NA)
 }))
}
ff<-list();cv<-list();k<-0;check<-list()
for(di in seq_along(ds)){
 d<-ds[[di]];id<-d$location_id[1]
 # Numerical equivalence of vars QR fitting and the original MTS normal-equation implementation.
 if(di%in%c(1,41,81,121,161,204)){
  invisible(capture.output(m<-MTS::VARX(log1p(as.matrix(d[,c('CKD','DD')])),p=2,xt=matrix(d$SDI,ncol=1),include.mean=TRUE)))
  invisible(capture.output(mp<-MTS::VARXpred(m,newxt=matrix(rep(tail(d$SDI,1),9),ncol=1),hstep=9)))
  z<-fcvar(d,9);check[[length(check)+1]]<-data.frame(location_id=id,max_abs_difference=max(abs(c(expm1(mp$pred[,1]),expm1(mp$pred[,2]))-z$median)),max_relative_difference=max(abs(c(expm1(mp$pred[,1]),expm1(mp$pred[,2]))-z$median)/pmax(1,abs(z$median))))
 }
 for(model in c('Original_VARX_p2','ARIMAX_log1p')){
  fn<-if(model=='Original_VARX_p2')fcvar else fcar
  for(sc in c('hold','trend')){
   z<-fn(d,9,sc);z$location_id<-id;z$location_name<-d$location_name[1];z$model<-model;z$scenario<-sc;ff[[length(ff)+1]]<-z
  }
  for(origin in c(2009,2012,2015,2018)){
   h<-min(9,2021-origin);train<-d[d$year<=origin,];test<-d[d$year>origin&d$year<=origin+h,]
   z<-fn(train,h);z$actual<-c(test$CKD,test$DD);z$location_id<-id;z$origin<-origin;z$model<-model
   z$APE<-100*abs(z$median-z$actual)/z$actual;z$APE_clipped<-100*abs(z$clipped-z$actual)/z$actual
   z$covered<-z$actual>=z$lower&z$actual<=z$upper;cv[[length(cv)+1]]<-z
  }
 }
 for(origin in c(2009,2012,2015,2018))for(dis in c('CKD','DD')){
  h<-min(9,2021-origin);train<-d[d$year<=origin,];test<-d[d$year>origin&d$year<=origin+h,];pred<-tail(train[[dis]],1)
  z<-data.frame(disease=dis,h=1:h,median=pred,lower=NA,upper=NA,clipped=pred,serial_p=NA,max_root=NA,
   actual=test[[dis]],location_id=id,origin=origin,model='Naive',APE=100*abs(pred-test[[dis]])/test[[dis]],APE_clipped=100*abs(pred-test[[dis]])/test[[dis]],covered=NA)
  cv[[length(cv)+1]]<-z
 }
 if(di%%20==0)cat('Original forecast audit',di,'/204\n')
}
f<-bind_rows(ff);v<-bind_rows(cv)
write_out(f,'A6_original_forecast_2030');write_out(v,'A7_forecast_validation');write_out(bind_rows(check),'A8_MTS_equivalence')
cvsummary<-function(v) v%>%group_by(model,disease)%>%summarise(n=n(),median_APE=median(APE),mean_APE=mean(APE),RMSE=sqrt(mean((median-actual)^2)),median_APE_clipped=median(APE_clipped),mean_APE_clipped=mean(APE_clipped),coverage95=if(all(is.na(covered)))NA_real_ else mean(covered),.groups='drop')
write_out(cvsummary(v[v$origin%in%c(2012,2015,2018)&v$h<=3,]),'A9_short_horizon_same_origins_summary')
write_out(v%>%group_by(model,disease,origin)%>%summarise(n=n(),median_APE=median(APE),mean_APE=mean(APE),RMSE=sqrt(mean((median-actual)^2)),coverage95=if(all(is.na(covered)))NA_real_ else mean(covered),.groups='drop'),'A10_validation_by_origin')
write_out(v%>%group_by(model,disease,h)%>%summarise(n=n(),median_APE=median(APE),mean_APE=mean(APE),RMSE=sqrt(mean((median-actual)^2)),coverage95=if(all(is.na(covered)))NA_real_ else mean(covered),.groups='drop'),'A11_validation_by_horizon')
write_out(cvsummary(v[v$h==9,]),'A12_nine_year_endpoint_summary')
write_out(cvsummary(v[v$origin==2009,]),'A13_pre_pandemic_validation_summary')
fwide<-tidyr::pivot_wider(f[f$h==9,],id_cols=c(location_id,location_name,model,scenario),names_from=disease,values_from=c(median,lower,upper,clipped))
fwide<-fwide%>%arrange(location_id)%>%group_by(model,scenario)%>%mutate(qCKD=ntile(median_CKD,4),qDD=ntile(median_DD,4),qCKD_clip=ntile(clipped_CKD,4),qDD_clip=ntile(clipped_DD,4),HH=qCKD==4&qDD==4,HH_clip=qCKD_clip==4&qDD_clip==4,same=qCKD==qDD,same_clip=qCKD_clip==qDD_clip)%>%ungroup()
write_out(fwide,'A14_2030_point_classification')
write_out(fwide%>%group_by(model,scenario)%>%summarise(n=n(),HH=sum(HH),HH_clip=sum(HH_clip),same=sum(same),same_clip=sum(same_clip),CKD_clipped=sum(abs(median_CKD-clipped_CKD)>1e-8),DD_clipped=sum(abs(median_DD-clipped_DD)>1e-8),.groups='drop'),'A15_2030_point_classification_summary')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('Temporal audit calculations complete\n')
