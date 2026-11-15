options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(vars);library(forecast)})
select<-dplyr::select
source('scripts/config.R')
out<-core_dir; dat<-read.csv(file.path(out,'S1_analysis_panel.csv'))%>%filter(sex_name=='Both')
set.seed(20261105)
savecsv<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
varforecast<-function(d,h=9,B=1000,scenario='hold'){
 Y<-cbind(CKD=log(d$CKD),DD=log(d$DD));dy<-diff(Y);xx<-matrix(diff(d$SDI),ncol=1,dimnames=list(NULL,'dSDI'))
 p<-as.integer(VARselect(dy,lag.max=2,type='const',exogen=xx)$selection['SC(n)']);fit<-VAR(dy,p=p,type='const',exogen=xx)
 future<-if(scenario=='hold')rep(tail(d$SDI,1),h) else pmin(1,pmax(0,tail(d$SDI,1)+seq_len(h)*mean(tail(diff(d$SDI),5))))
 dx<-diff(c(tail(d$SDI,1),future));cf<-Bcoef(fit); n<-nrow(dy);hist<-dy;point<-matrix(NA,h,2)
 for(t in seq_len(h)){v<-c(unlist(lapply(seq_len(p),function(j)hist[nrow(hist)+1-j,])),1,dx[t]);new<-as.numeric(cf%*%v);hist<-rbind(hist,new);point[t,]<-colSums(tail(hist,t))+Y[nrow(Y),]}
 residual<-resid(fit);residual<-sweep(residual,2,colMeans(residual)); paths<-array(NA,c(h,2,B))
 for(b in seq_len(B)){hist<-dy;cur<-Y[nrow(Y),];ee<-residual[sample.int(nrow(residual),h,TRUE),,drop=FALSE];for(t in seq_len(h)){v<-c(unlist(lapply(seq_len(p),function(j)hist[nrow(hist)+1-j,])),1,dx[t]);new<-as.numeric(cf%*%v)+ee[t,];hist<-rbind(hist,new);cur<-cur+new;paths[t,,b]<-cur}}
 stable<-max(roots(fit))<1;serial<-tryCatch(serial.test(fit,lags.pt=4,type='PT.adjusted')$serial$p.value,error=function(e)NA_real_)
 bind_rows(lapply(1:2,function(j)data.frame(disease=c('CKD','DD')[j],h=seq_len(h),median=exp(point[,j]),lower=exp(apply(paths[,j,],1,quantile,.025)),upper=exp(apply(paths[,j,],1,quantile,.975)),p=p,order=paste0('VAR(',p,') on differences'),stable=stable,serial_p=serial,clipped=FALSE)))
}
arforecast<-function(d,h=9,scenario='hold'){
 future<-if(scenario=='hold')rep(tail(d$SDI,1),h) else pmin(1,pmax(0,tail(d$SDI,1)+seq_len(h)*mean(tail(diff(d$SDI),5))))
 bind_rows(lapply(c('CKD','DD'),function(dis){fit<-auto.arima(log(d[[dis]]),xreg=matrix(d$SDI,ncol=1,dimnames=list(NULL,'SDI')),seasonal=FALSE,max.p=2,max.q=2,max.d=2,stepwise=TRUE,approximation=FALSE);f<-forecast(fit,xreg=matrix(future,ncol=1,dimnames=list(NULL,'SDI')),h=h,level=95);oo<-arimaorder(fit); lag<-min(6,length(resid(fit))-1);fp<-oo[1]+oo[3];pv<-Box.test(resid(fit),lag=lag,type='Ljung-Box',fitdf=fp)$p.value;data.frame(disease=dis,h=seq_len(h),median=exp(as.numeric(f$mean)),lower=exp(as.numeric(f$lower)),upper=exp(as.numeric(f$upper)),p=oo[1],order=paste0('ARIMA(',paste(oo,collapse=','),')'),stable=NA,serial_p=pv,clipped=FALSE)}))
}
allf<-list(); val<-list();err<-list();i<-0
for(id in sort(unique(dat$location_id))){
 d<-dat[dat$location_id==id,];d<-d[order(d$year),];i<-i+1
 for(model in c('VARX','ARIMAX')){
  fn<-if(model=='VARX')varforecast else arforecast
  for(sc in c('hold','trend')){
   z<-tryCatch(fn(d,h=9,scenario=sc),error=function(e){err[[length(err)+1]]<<-data.frame(location_id=id,model=model,stage=sc,message=conditionMessage(e));NULL})
   if(!is.null(z)){z$location_id<-id;z$location_name<-d$location_name[1];z$model<-model;z$scenario<-sc;z$year<-2021+z$h;allf[[length(allf)+1]]<-z}
  }
  for(origin in c(2012,2015,2018)){
   train<-filter(d,year<=origin);test<-filter(d,year>origin,year<=origin+3)
   z<-tryCatch(if(model=='VARX')fn(train,h=3,B=300) else fn(train,h=3),error=function(e){err[[length(err)+1]]<<-data.frame(location_id=id,model=model,stage=paste0('CV',origin),message=conditionMessage(e));NULL})
   if(!is.null(z)){z$actual<-c(test$CKD,test$DD);z$location_id<-id;z$origin<-origin;z$model<-model;z$APE<-100*abs(z$actual-z$median)/z$actual;z$covered<-z$actual>=z$lower&z$actual<=z$upper;val[[length(val)+1]]<-z}
  }
 }
 # The last-observation benchmark is computed on the identical origins and outcomes.
 for(origin in c(2012,2015,2018))for(dis in c('CKD','DD')){test<-filter(d,year>origin,year<=origin+3);pred<-d[[dis]][d$year==origin];val[[length(val)+1]]<-data.frame(disease=dis,h=1:3,median=pred,lower=NA,upper=NA,p=NA,order='Last observation',stable=NA,serial_p=NA,clipped=FALSE,actual=test[[dis]],location_id=id,origin=origin,model='Naive',APE=100*abs(test[[dis]]-pred)/test[[dis]],covered=NA)}
 if(i%%20==0)cat('Forecasts',i,'/204\n')
}
f<-bind_rows(allf);v<-bind_rows(val);savecsv(f,'S14_forecasts');savecsv(v,'S15_rolling_validation');savecsv(bind_rows(err),'S16_forecast_failures')
summ<-v%>%group_by(model,disease)%>%summarise(n=n(),median_APE=median(APE),mean_APE=mean(APE),RMSE=sqrt(mean((actual-median)^2)),coverage95=if(all(is.na(covered)))NA_real_ else mean(covered,na.rm=TRUE),.groups='drop');savecsv(summ,'S17_forecast_validation_summary')
diagn<-f%>%filter(scenario=='hold',h==1)%>%group_by(model,disease)%>%summarise(n=n(),residual_autocorrelation=sum(serial_p<.05,na.rm=TRUE),unstable=if(all(is.na(stable)))NA_integer_ else sum(!stable,na.rm=TRUE),.groups='drop');savecsv(diagn,'S18_forecast_diagnostics')
wide<-f%>%filter(year==2030)%>%select(location_id,location_name,model,scenario,disease,median,lower,upper)%>%pivot_wider(names_from=disease,values_from=c(median,lower,upper))%>%arrange(location_id)%>%group_by(model,scenario)%>%mutate(qCKD=ntile(median_CKD,4),qDD=ntile(median_DD,4),high_high=qCKD==4&qDD==4,same_quartile=qCKD==qDD)%>%ungroup();savecsv(wide,'S19_forecast_classification_2030')
savecsv(wide%>%group_by(model,scenario)%>%summarise(high_high=sum(high_high),same_quartile=sum(same_quartile),.groups='drop'),'S20_forecast_classification_summary')
cat('FORECAST COMPLETE\n')
