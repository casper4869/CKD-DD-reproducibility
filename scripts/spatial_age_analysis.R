options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(spgwr);library(spdep);library(sandwich)})
source('scripts/config.R')
out <- spatial_dir
set.seed(2026110403)
savecsv <- function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
perm_two_sided_plus1 <- function(mc){
  observed <- tail(mc$res,1);simulated <- head(mc$res,-1)
  min(1,2*min((1+sum(simulated>=observed))/(length(simulated)+1),(1+sum(simulated<=observed))/(length(simulated)+1)))
}
# No p-values for within-country longitudinal correlations: annual estimates are
# serially dependent and their model-based uncertainty draws are unavailable.
ages <- c('5-9 years','10-14 years',paste0(seq(15,90,5),'-',seq(19,94,5),' years'),'95+ years')
age0 <- as.integer(sub('-.*','',sub('\\+ years','',ages)))
ac <- list()
for(s in c('female','male')){
  a <- read.csv(file.path(private,paste0('CKD_DD_',s,'_LONG_1992_2021.csv')))
  a <- a[a$age_name%in%ages,]
  stopifnot(nrow(a)==204*19*30,!anyDuplicated(a[,c('location_name','age_name','year')]),all(a$ckd_val>0),all(a$dd_val>0))
  # Verify all 2021 values independently against the extracted raw GBD input.
  raw21 <- read.csv(file.path(core_dir,'S5_age_source.csv'))
  v <- inner_join(filter(a,year==2021),filter(raw21,tolower(sex_name)==s),by=c('location_name','age_name'))
  stopifnot(nrow(v)==204*19,max(abs(v$ckd_val-v$CKD))<1e-8,max(abs(v$dd_val-v$DD))<1e-8)
  ac[[s]] <- a%>%arrange(year)%>%group_by(location_name,age_name)%>%summarise(n_year=n(),raw_r=cor(ckd_val,dd_val),dlog_r=cor(diff(log(ckd_val)),diff(log(dd_val))),.groups='drop')%>%mutate(sex=s)
}
ac <- bind_rows(ac);savecsv(ac,'age_withincountry_time_correlations')
asum <- ac%>%group_by(sex,age_name)%>%summarise(n_country=n(),median_raw_r=median(raw_r),q25_raw_r=quantile(raw_r,.25),q75_raw_r=quantile(raw_r,.75),median_dlog_r=median(dlog_r),q25_dlog_r=quantile(dlog_r,.25),q75_dlog_r=quantile(dlog_r,.75),.groups='drop')%>%mutate(age_start=age0[match(age_name,ages)])%>%arrange(sex,age_start)
savecsv(asum,'age_time_summary')
contr <- list()
for(s in c('female','male'))for(estimand in c('raw_r','dlog_r')){
  a <- ac[ac$sex==s,];a$age_start <- age0[match(a$age_name,ages)]
  a$band <- ifelse(a$age_start>=50,'older',ifelse(a$age_start>=15,'younger','excluded'))
  w <- a%>%filter(band!='excluded')%>%group_by(location_name,band)%>%summarise(mean_r=mean(.data[[estimand]]),.groups='drop')%>%pivot_wider(names_from=band,values_from=mean_r)
  w$difference <- w$older-w$younger
  # Conditional, descriptive country-resampling interval. Does not propagate
  # GBD estimation uncertainty, annual-series sampling or spatial dependence.
  boot <- replicate(2000,mean(sample(w$difference,nrow(w),replace=TRUE)))
  ci <- quantile(boot,c(.025,.975))
  contr[[paste(s,estimand)]] <- data.frame(sex=s,estimand=estimand,n_country=nrow(w),mean_difference=mean(w$difference),median_difference=median(w$difference),conditional_boot_lower=ci[1],conditional_boot_upper=ci[2],positive_countries=sum(w$difference>0))
}
savecsv(bind_rows(contr),'age_time_contrast_descriptive')
# Audit the first-round high-high computation and ties.
dat <- read.csv(file.path(core_dir,'S1_analysis_panel.csv'))
qt <- read.csv(file.path(core_dir,'S2_quartile_country_year.csv'))
qa <- dat%>%group_by(sex_name,year)%>%summarise(n=n(),CKD_ties=n()-n_distinct(CKD),DD_ties=n()-n_distinct(DD),.groups='drop')
stopifnot(all(qa$n==204),all(qa$CKD_ties==0),all(qa$DD_ties==0))
savecsv(qa,'quartile_tie_audit')
hh <- qt%>%filter(year==2021)%>%group_by(sex_name)%>%summarise(n=n(),CKD_top=sum(q_CKD==4),DD_top=sum(q_DD==4),observed=sum(q_CKD==4&q_DD==4),.groups='drop')%>%mutate(expected=CKD_top*DD_top/n,p_upper=phyper(observed-1,CKD_top,n-CKD_top,DD_top,lower.tail=FALSE),q_BH=p.adjust(p_upper,'BH'))
savecsv(hh,'quartile_overlap_independent_check')
cat('Age and quartile independent checks complete\n')
# Join only public geometry coordinates to the locally generated disease panel.
geo0 <- read.csv(file.path(core_dir,'S11_spatial_inputs.csv'))
ref <- read.csv(file.path(reference_dir,'geometry_reference_points.csv'))
geo0 <- merge(geo0,ref,by=c('location_id','location_name'),sort=FALSE)
geo0 <- geo0[order(geo0$location_id),]
stopifnot(nrow(geo0)==204)
write.csv(geo0,file.path(out,'spatial_inputs_geometry.csv'),row.names=FALSE)

geo <- read.csv(file.path(out,'spatial_inputs_geometry.csv'))
stopifnot(nrow(geo)==204,!anyDuplicated(geo$location_id))
spresults <- list();diag <- list();olsresults <- list()
for(coordinate in c('original_city_points','land_geometry_points')){
  coords <- if(coordinate=='original_city_points')as.matrix(geo[,c('lng','lat')])else as.matrix(geo[,c('geometry_lng','geometry_lat')])
  kn <- spdep::knearneigh(coords,k=4,longlat=TRUE);nb <- spdep::knn2nb(kn,sym=TRUE);lw <- spdep::nb2listw(nb,style='W')
  for(target in c('CKD','DD'))for(spec in c('SDI_PM25','SDI_only')){
    source <- if(target=='CKD')'DD'else'CKD'
    form <- as.formula(paste0('z',target,' ~ z',source,' + zSDI',if(spec=='SDI_PM25')' + zPM25'else''))
    bw <- spgwr::gwr.sel(form,data=geo,coords=coords,longlat=TRUE,verbose=FALSE)
    fit <- spgwr::gwr(form,data=geo,coords=coords,bandwidth=bw,longlat=TRUE,hatmatrix=TRUE,se.fit=TRUE)
    rr <- as.data.frame(fit$SDF);nm <- paste0('z',source);se <- rr[[paste0(nm,'_se_EDF')]];edf <- fit$results$edf
    p <- 2*pt(-abs(rr[[nm]]/se),df=edf)
    moran <- spdep::moran.mc(rr$gwr.e,lw,nsim=999,alternative='two.sided')
    design <- model.matrix(form,geo)
    cn <- sapply(seq_len(nrow(geo)),function(i){dist <- sp::spDistsN1(coords,coords[i,],longlat=TRUE);w <- spgwr::gwr.Gauss(dist^2,bw);kappa(design*sqrt(w),exact=TRUE)})
    key <- paste(coordinate,target,spec)
    spresults[[key]] <- data.frame(location_id=geo$location_id,location_name=geo$location_name,coordinate=coordinate,target=target,spec=spec,coefficient=rr[[nm]],se_EDF=se,p_conditional=p,q_BH_204=p.adjust(p,'BH'),q_BY_204=p.adjust(p,'BY'),local_R2=rr$localR2,condition_number=cn,residual=rr$gwr.e)
    diag[[key]] <- data.frame(coordinate=coordinate,target=target,spec=spec,n=204,bandwidth_km=bw,edf=edf,R2=1-fit$results$rss/sum((geo[[paste0('z',target)]]-mean(geo[[paste0('z',target)]]))^2),AICc=fit$results$AICc,negative=sum(rr[[nm]]<0),positive=sum(rr[[nm]]>0),BH204=sum(p.adjust(p,'BH')<.05),BY204=sum(p.adjust(p,'BY')<.05),condition_over30=sum(cn>30),Moran_I=unname(moran$statistic),Moran_p_spdep=moran$p.value,Moran_p_plus1_two_sided=perm_two_sided_plus1(moran))
    # OLS coefficient intervals with heteroscedasticity-robust HC3 SE;
    # spatial dependence is reported rather than ignored as a resolved issue.
    ol <- lm(form,data=geo);seol <- sqrt(diag(sandwich::vcovHC(ol,type='HC3')))[nm];b <- coef(ol)[nm];crit <- qt(.975,df.residual(ol));mt <- spdep::moran.mc(resid(ol),lw,nsim=999,alternative='two.sided')
    olsresults[[key]] <- data.frame(coordinate=coordinate,target=target,spec=spec,coefficient=b,HC3_se=seol,HC3_lower=b-crit*seol,HC3_upper=b+crit*seol,p_HC3=2*pt(-abs(b/seol),df.residual(ol)),R2=summary(ol)$r.squared,Moran_I=unname(mt$statistic),Moran_p_spdep=mt$p.value,Moran_p_plus1_two_sided=perm_two_sided_plus1(mt))
    cat(key,'complete\n')
  }
}
spresults <- bind_rows(spresults)%>%group_by(coordinate,spec)%>%mutate(q_BH_408=p.adjust(p_conditional,'BH'),q_BY_408=p.adjust(p_conditional,'BY'))%>%ungroup()
savecsv(spresults,'GWR_sensitivity_coefficients');savecsv(bind_rows(diag),'GWR_sensitivity_diagnostics');savecsv(bind_rows(olsresults),'global_OLS_HC3_comparator')
cons <- list()
for(target in c('CKD','DD'))for(spec in c('SDI_PM25','SDI_only')){
  a <- spresults%>%filter(.data$target==.env$target,.data$spec==.env$spec)%>%select(location_id,coordinate,coefficient)%>%pivot_wider(names_from=coordinate,values_from=coefficient)
  cons[[paste(target,spec)]] <- data.frame(target=target,spec=spec,n=204,same_sign=sum(sign(a$original_city_points)==sign(a$land_geometry_points)),coefficient_Pearson_r=cor(a$original_city_points,a$land_geometry_points),median_absolute_change=median(abs(a$original_city_points-a$land_geometry_points)),max_absolute_change=max(abs(a$original_city_points-a$land_geometry_points)))
}
savecsv(bind_rows(cons),'GWR_coordinate_agreement')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
