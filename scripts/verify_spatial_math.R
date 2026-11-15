options(stringsAsFactors=FALSE)
source('scripts/config.R')
out <- spatial_dir
geo <- read.csv(file.path(out,'spatial_inputs_geometry.csv'))
rr <- read.csv(file.path(out,'GWR_sensitivity_coefficients.csv'))
dd <- read.csv(file.path(out,'GWR_sensitivity_diagnostics.csv'))
checks <- list()
for(k in seq_len(nrow(dd))){
 d <- dd[k,];a <- rr[rr$coordinate==d$coordinate&rr$target==d$target&rr$spec==d$spec,]
 stopifnot(identical(a$location_id,geo$location_id))
 src <- if(d$target=='CKD')'DD'else'CKD'
 f <- as.formula(paste0('z',d$target,' ~ z',src,' + zSDI',if(d$spec=='SDI_PM25')' + zPM25'else''))
 X <- model.matrix(f,geo);y <- geo[[paste0('z',d$target)]]
 coords <- if(d$coordinate=='original_city_points')as.matrix(geo[,c('lng','lat')])else as.matrix(geo[,c('geometry_lng','geometry_lat')])
 b <- sapply(seq_len(nrow(geo)),function(i){dist <- sp::spDistsN1(coords,coords[i,],longlat=TRUE);w <- exp(-.5*(dist/d$bandwidth_km)^2);stats::lm.wfit(X,y,w)$coefficients[paste0('z',src)]})
 err <- max(abs(b-a$coefficient))
 stopifnot(err<1e-8)
 checks[[k]] <- data.frame(coordinate=d$coordinate,target=d$target,spec=d$spec,n_coefficient_checks=204,max_absolute_error=err)
}
old <- read.csv(file.path(core_dir,'S12_GWR_coefficients.csv'))
orig <- subset(rr,coordinate=='original_city_points'&spec=='SDI_PM25')
old <- subset(old,spec=='primary')
m <- merge(orig,old,by=c('target','location_id'),suffixes=c('_audit','_first_revision'))
stopifnot(nrow(m)==408,max(abs(m$coefficient_audit-m$coefficient_first_revision))<1e-8)
write.csv(do.call(rbind,checks),file.path(out,'independent_weighted_regression_check.csv'),row.names=FALSE)
writeLines(c('8 models x 204 local target coefficients independently reproduced using lm.wfit and explicit Gaussian exp(-0.5*(distance/bandwidth)^2) weights; max abs error below 1e-8.','408 original-city SDI+PM local coefficients agree with first-revision S12 within 1e-8 (CSV round-trip inputs).','No primary result was selected using these audit results.'),file.path(out,'numerical_validation.txt'))
cat('Spatial mathematical checks passed\n')
