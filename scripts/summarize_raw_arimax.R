options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
source('scripts/config.R');out<-temporal_dir
readpart<-function(prefix)bind_rows(lapply(1:2,function(p)read.csv(file.path(out,paste0(prefix,'_part',p,'.csv')))))
f<-readpart('A37_raw_ARIMAX_forecasts');v<-readpart('A38_raw_ARIMAX_validation')
write.csv(f,file.path(out,'A37_raw_ARIMAX_forecasts.csv'),row.names=FALSE,na='')
write.csv(v,file.path(out,'A38_raw_ARIMAX_validation.csv'),row.names=FALSE,na='')
stopifnot(length(unique(f$location_id))==204,nrow(f)==204*2*2*9,nrow(v)==204*2*27)
others<-read.csv(file.path(out,'A7_forecast_validation.csv'))%>%filter(model%in%c('Original_VARX_p2','Naive'))
logdefault<-read.csv(file.path(out,'A28_original_ARIMAX_validation.csv'))
v<-bind_rows(v,others,logdefault)
sumcv<-function(z)z%>%group_by(model,disease)%>%summarise(n=n(),median_APE=median(APE),mean_APE=mean(APE),RMSE=sqrt(mean((median-actual)^2)),coverage95=if(all(is.na(covered)))NA_real_ else mean(covered),.groups='drop')
write.csv(sumcv(v%>%filter(origin%in%c(2012,2015,2018),h<=3)),file.path(out,'A40_all_fixed_spec_short_validation.csv'),row.names=FALSE,na='')
write.csv(sumcv(v%>%filter(h==9)),file.path(out,'A41_all_fixed_spec_nine_year_endpoint.csv'),row.names=FALSE,na='')
fw<-f%>%filter(h==9)%>%select(location_id,location_name,scenario,disease,median,lower,upper,order,serial_p)%>%pivot_wider(names_from=disease,values_from=c(median,lower,upper,order,serial_p))%>%arrange(location_id)%>%group_by(scenario)%>%mutate(qCKD=ntile(median_CKD,4),qDD=ntile(median_DD,4),HH=qCKD==4&qDD==4,same=qCKD==qDD)%>%ungroup()
write.csv(fw%>%group_by(scenario)%>%summarise(n=n(),HH=sum(HH),same=sum(same),.groups='drop'),file.path(out,'A42_raw_ARIMAX_2030_summary.csv'),row.names=FALSE,na='')
write.csv(f%>%filter(h==1,scenario=='hold')%>%group_by(disease)%>%summarise(n=n(),residual_reject=sum(serial_p<.05,na.rm=TRUE),residual_NA=sum(is.na(serial_p)),invalid_point=sum(!is.finite(median)|median<0),.groups='drop'),file.path(out,'A44_raw_ARIMAX_diagnostics.csv'),row.names=FALSE,na='')

# Validation of the manuscript's actual rank-based target, using the same frozen predictions.
# This adds no fitted candidate and makes no independence or calibrated-probability claim.
wide<-v%>%select(location_id,model,origin,h,disease,median,actual)%>%pivot_wider(names_from=disease,values_from=c(median,actual))%>%arrange(location_id)%>%group_by(model,origin,h)%>%mutate(pCKD=ntile(median_CKD,4),pDD=ntile(median_DD,4),aCKD=ntile(actual_CKD,4),aDD=ntile(actual_DD,4),predHH=pCKD==4&pDD==4,actualHH=aCKD==4&aDD==4,predSame=pCKD==pDD,actualSame=aCKD==aDD,predType=sign(pCKD-pDD),actualType=sign(aCKD-aDD))
rankvalid<-wide%>%summarise(n=n(),exact_two_quartiles_accuracy=mean(pCKD==aCKD&pDD==aDD),three_category_accuracy=mean(predType==actualType),pred_HH=sum(predHH),actual_HH=sum(actualHH),HH_TP=sum(predHH&actualHH),HH_precision=if(sum(predHH)>0)sum(predHH&actualHH)/sum(predHH) else NA_real_,HH_recall=if(sum(actualHH)>0)sum(predHH&actualHH)/sum(actualHH) else NA_real_,HH_Jaccard=sum(predHH&actualHH)/sum(predHH|actualHH),pred_same=sum(predSame),actual_same=sum(actualSame),.groups='drop')
write.csv(rankvalid,file.path(out,'A45_quartile_target_validation.csv'),row.names=FALSE,na='')
