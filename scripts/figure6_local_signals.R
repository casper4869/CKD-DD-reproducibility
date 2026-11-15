options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(maps)})
source('scripts/config.R');out<-file.path(result,'figure_sources');dir.create(out,showWarnings=FALSE)
g <- read.csv(file.path(core_dir,'S12_GWR_coefficients.csv'))
gg <- read.csv(file.path(spatial_dir,'GWR_sensitivity_coefficients.csv'))%>%filter(spec=='SDI_PM25')
world <- map_data('world')%>%filter(region!='Antarctica')
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),legend.title=element_text(size=8),legend.text=element_text(size=7),plot.title=element_text(size=9,face='bold'),strip.text=element_text(size=8,face='bold'),legend.position='bottom',plot.tag=element_text(size=11,face='bold')))
pointmap <- function(d,var,label){ggplot()+geom_polygon(data=world,aes(long,lat,group=group),fill='#F0F1F2',colour='white',linewidth=.12)+geom_point(data=d,aes(x=lng,y=lat,colour=.data[[var]]),size=1.6,alpha=.9)+coord_quickmap(xlim=c(-180,180),ylim=c(-60,85),expand=FALSE)+theme_void(base_size=8)+theme(legend.position='bottom',legend.key.width=grid::unit(11,'mm'),legend.title=element_text(size=7),legend.text=element_text(size=7),plot.title=element_text(size=9,face='bold'))}
lim <- max(abs(filter(g,spec=='primary')$coefficient))
p <- list();source <- list();ring_counts <- list()
for(tar in c('CKD','DD')){
 z <- filter(g,target==tar,spec=='primary')
 stopifnot(nrow(z)==204,!anyNA(z$q_BH))
 sig <- filter(z,q_BH<.05)
 p[[length(p)+1]] <- pointmap(z,'coefficient','Local coefficient')+
  scale_colour_gradient2(low='#3B6E9A',mid='white',high='#B45C58',midpoint=0,limits=c(-lim,lim),name='Local coefficient')+
  geom_point(data=sig,aes(x=lng,y=lat),shape=1,colour='black',size=1.9,stroke=.28,inherit.aes=FALSE)+ggtitle(paste(tar,'as outcome'))
 source[[tar]] <- transform(z,black_ring=q_BH<.05)
 ring_counts[[tar]] <- data.frame(target=tar,n_shown=nrow(z),ring_n=nrow(sig),ring_positive=sum(sig$coefficient>0),ring_negative=sum(sig$coefficient<0),color_limit=lim)
}
for(tar in c('CKD','DD')){
 z <- filter(gg,target==tar)%>%select(location_id,coordinate,coefficient)%>%pivot_wider(names_from=coordinate,values_from=coefficient)
 stopifnot(nrow(z)==204)
 p[[length(p)+1]] <- ggplot(z,aes(original_city_points,land_geometry_points))+geom_abline(slope=1,intercept=0,linetype=2,linewidth=.3,colour='grey50')+geom_point(size=1.3,colour='#617D8B',alpha=.7)+labs(x='Archived city / reference points',y='Geometry reference points',title=paste(tar,'coordinate sensitivity'))
}
plot <- wrap_plots(p,ncol=2)+plot_annotation(tag_levels='A')
stem <- file.path(figure_dir,'Figure_6')
svglite::svglite(paste0(stem,'.svg'),width=183/25.4,height=150/25.4);print(plot);dev.off()
cairo_pdf(paste0(stem,'.pdf'),width=183/25.4,height=150/25.4,family='Arial');print(plot);dev.off()
ragg::agg_tiff(paste0(stem,'.tif'),width=183,height=150,units='mm',res=600,compression='lzw');print(plot);dev.off()
ragg::agg_png(paste0(stem,'.png'),width=183,height=150,units='mm',res=300);print(plot);dev.off()
write.csv(bind_rows(source),file.path(out,'Figure_6_map_source_with_all_ring_flags.csv'),row.names=FALSE)
write.csv(bind_rows(ring_counts),file.path(out,'Figure_6_ring_audit.csv'),row.names=FALSE)
write.csv(data.frame(figure='Figure 6',change='Presentation-only: overlay black open rings on every original-coordinate local coefficient whose existing BH204 q<0.05, positive or negative. All 204 locations remain in each map. Coordinate-sensitivity panels unchanged.',backend='R',width_mm=183,height_mm=150,tiff_dpi=600,png_dpi=300,colour_limit_unchanged=lim),file.path(out,'Figure_6_contract.csv'),row.names=FALSE)
writeLines(c('Figure 6. Local ecological coefficients and coordinate sensitivity.','A-B, local other-disease coefficients in separate 2021 CKD-outcome and DD-outcome GWR models, adjusted for SDI and PM2.5, using archived reference points. All 204 countries and territories are displayed in each panel. Variables were standardized across countries, coefficients varied locally, and fixed Gaussian kernels used cross-validated great-circle bandwidths. Color denotes coefficient sign and magnitude on the unchanged common scale. Black open rings mark all local coefficients with approximate conditional BH-adjusted q<0.05 within the 204-location outcome model, regardless of coefficient sign; they do not establish causal effects or account for bandwidth selection and all residual dependence. C-D, coefficients from the same adjustment set using archived city/reference points versus documented geometry-based points; dashed lines indicate equality. These are cross-sectional ecological contrasts. Full coefficients, approximate uncertainty, local fit, condition numbers, population and SDI-only sensitivities, and Moran diagnostics are supplied as derived data.'),file.path(out,'Figure_6_caption.txt'))
cat('Figure 6 presentation-only update complete\n')
