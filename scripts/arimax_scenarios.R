options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(forecast);library(dplyr);library(tidyr)})
source('scripts/config.R');out<-temporal_dir
args<-commandArgs(trailingOnly=TRUE);part<-if(length(args))as.integer(args[1]) else 1L
rawscale<-length(args)>=2&&args[2]=='raw'
dat<-read.csv(file.path(core_dir,'S1_analysis_panel.csv'));dat<-dat[dat$sex_name=='Both',]
ds<-split(dat,dat$location_id);ds<-ds[seq_along(ds)%%2==part%%2]
ff<-list();vv<-list();errors<-list()
one<-function(d,dis,h,sc='hold'){
 # Archived automatic ARIMAX specification; log1p is the manuscript scenario model.
 y<-if(rawscale)d[[dis]] else log1p(d[[dis]]);xreg<-as.numeric(d$SDI)
 fit<-auto.arima(y,xreg=xreg)
 future<-if(sc=='hold')rep(tail(d$SDI,1),h) else pmin(1,pmax(0,tail(d$SDI,1)+seq_len(h)*mean(tail(diff(d$SDI),5))))
 f<-forecast(fit,xreg=as.numeric(future),h=h,level=95)
 inv<-if(rawscale)identity else expm1
 point<-inv(as.numeric(f$mean));oo<-arimaorder(fit)
 pv<-tryCatch(Box.test(resid(fit),lag=6,fitdf=oo[1]+oo[3],type='Ljung-Box')$p.value,error=function(e)NA_real_)
 data.frame(location_id=d$location_id[1],location_name=d$location_name[1],disease=dis,h=1:h,median=point,lower=inv(as.numeric(f$lower)),upper=inv(as.numeric(f$upper)),clipped=if(rawscale)point else pmin(1.5*max(d[[dis]]),pmax(min(d[[dis]]),point)),order=paste0('ARIMA(',paste(oo,collapse=','),')'),serial_p=as.numeric(pv),model=if(rawscale)'Original_ARIMAX_raw' else 'Original_ARIMAX_default')
}
for(i in seq_along(ds)){
 d<-ds[[i]];d<-d[order(d$year),]
 for(dis in c('CKD','DD')){
  for(sc in c('hold','trend')){
   z<-tryCatch(one(d,dis,9,sc),error=function(e){errors[[length(errors)+1]]<<-data.frame(location_id=d$location_id[1],disease=dis,stage=sc,error=conditionMessage(e));NULL})
   if(!is.null(z)){z$scenario<-sc;ff[[length(ff)+1]]<-z}
  }
  for(origin in c(2009,2012,2015,2018)){
   h<-min(9,2021-origin);train<-d[d$year<=origin,];test<-d[d$year>origin&d$year<=origin+h,]
   z<-tryCatch(one(train,dis,h),error=function(e){errors[[length(errors)+1]]<<-data.frame(location_id=d$location_id[1],disease=dis,stage=paste0('CV',origin),error=conditionMessage(e));NULL})
   if(!is.null(z)){z$origin<-origin;z$actual<-test[[dis]];z$APE<-100*abs(z$median-z$actual)/z$actual;z$APE_clipped<-100*abs(z$clipped-z$actual)/z$actual;z$covered<-z$actual>=z$lower&z$actual<=z$upper;vv[[length(vv)+1]]<-z}
  }
 }
 if(i%%10==0)cat('Original ARIMAX part',part,'country',i,'/',length(ds),'\n')
}
write.csv(bind_rows(ff),file.path(out,paste0(if(rawscale)'A37_raw_ARIMAX_forecasts_part' else 'A27_original_ARIMAX_forecasts_part',part,'.csv')),row.names=FALSE,na='')
write.csv(bind_rows(vv),file.path(out,paste0(if(rawscale)'A38_raw_ARIMAX_validation_part' else 'A28_original_ARIMAX_validation_part',part,'.csv')),row.names=FALSE,na='')
write.csv(bind_rows(errors),file.path(out,paste0(if(rawscale)'A39_raw_ARIMAX_errors_part' else 'A29_original_ARIMAX_errors_part',part,'.csv')),row.names=FALSE,na='')
cat('Original ARIMAX part',part,'complete\n')
