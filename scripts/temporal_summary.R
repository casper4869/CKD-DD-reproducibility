options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
source('scripts/config.R');out<-temporal_dir
d<-read.csv(file.path(out,'A1_model_tests_country.csv'))
d<-bind_rows(lapply(split(d,d$model),function(z){n<-nrow(z);q<-p.adjust(c(z$joint_CKD,z$joint_DD),'BH');z$q_joint_CKD_pooled<-head(q,n);z$q_joint_DD_pooled<-tail(q,n);z}))
write.csv(d,file.path(out,'A1_model_tests_country.csv'),row.names=FALSE,na='')
s<-d%>%group_by(model)%>%summarise(n=n(),p1=sum(p_base==1),p2=sum(p_base==2),
 joint_CKD_pooled=sum(q_joint_CKD_pooled<.05,na.rm=TRUE),joint_DD_pooled=sum(q_joint_DD_pooled<.05,na.rm=TRUE),joint_both_pooled=sum(q_joint_CKD_pooled<.05&q_joint_DD_pooled<.05,na.rm=TRUE),
 joint_both_separate=sum(q_joint_CKD_separate<.05&q_joint_DD_separate<.05,na.rm=TRUE),
 specific_both_separate=sum(q_CKD_DD_separate<.05&q_DD_CKD_separate<.05),
 BH_CKD_to_DD=sum(q_CKD_DD<.05),BH_DD_to_CKD=sum(q_DD_CKD<.05),BH_both=sum(q_CKD_DD<.05&q_DD_CKD<.05),
 BH_CKD_only=sum(q_CKD_DD<.05&q_DD_CKD>=.05),BH_DD_only=sum(q_CKD_DD>=.05&q_DD_CKD<.05),BH_neither=sum(q_CKD_DD>=.05&q_DD_CKD>=.05),
 BY_both=sum(qBY_CKD_DD<.05&qBY_DD_CKD<.05),HC3_BH_CKD_to_DD=sum(q_HC3_CKD_DD<.05),HC3_BH_DD_to_CKD=sum(q_HC3_DD_CKD<.05),
 HC3_BH_both=sum(q_HC3_CKD_DD<.05&q_HC3_DD_CKD<.05),serial_reject=sum(serial_p<.05,na.rm=TRUE),stable=sum(max_root<1),
 min_df=min(df_residual),median_condition=median(design_condition),.groups='drop')
write.csv(s,file.path(out,'A2_model_summary.csv'),row.names=FALSE,na='')
fc<-read.csv(file.path(out,'A6_original_forecast_2030.csv'));cv<-read.csv(file.path(out,'A7_forecast_validation.csv'))
actual<-read.csv(file.path(core_dir,'S1_analysis_panel.csv'))%>%filter(sex_name=='Both',year==2021)%>%arrange(location_id)%>%mutate(qCKD=ntile(CKD,4),qDD=ntile(DD,4),HH=qCKD==4&qDD==4,same=qCKD==qDD)
write.csv(actual,file.path(out,'A23_baseline_2021_classification.csv'),row.names=FALSE,na='')
write.csv(fc%>%filter(scenario=='hold',model=='Original_VARX_p2',h==9)%>%group_by(disease)%>%summarise(n=n(),max_root_max=max(max_root),stable=sum(max_root<1),n_point_over_1e6=sum(median>1e6),n_point_over_1e10=sum(median>1e10),median_point=median(median),max_point=max(median),min_point=min(median),.groups='drop'),file.path(out,'A24_original_2030_extremes.csv'),row.names=FALSE,na='')
write.csv(cv%>%filter(h==9)%>%group_by(model,disease,origin)%>%summarise(n=n(),median_APE=median(APE),mean_APE=mean(APE),RMSE=sqrt(mean((median-actual)^2)),coverage95=if(all(is.na(covered)))NA_real_ else mean(covered),.groups='drop'),file.path(out,'A25_nine_year_endpoints_by_origin.csv'),row.names=FALSE,na='')
write.csv(cv%>%filter(origin%in%c(2012,2015,2018),h<=3)%>%group_by(model,disease,location_id)%>%summarise(mean_APE=mean(APE),RMSE=sqrt(mean((median-actual)^2)),.groups='drop')%>%pivot_wider(names_from=model,values_from=c(mean_APE,RMSE)),file.path(out,'A26_short_horizon_country_metrics.csv'),row.names=FALSE,na='')
cat('Baseline HH=',sum(actual$HH),'same=',sum(actual$same),'\n')
