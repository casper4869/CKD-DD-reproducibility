options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
source('scripts/config.R')
out <- spatial_dir
dir.create(file.path(out,'figures'),showWarnings=FALSE)
contract <- data.frame(figure='Figure S3',claim='Age patterns of national time-series correlation depend on preprocessing; they do not imply individual-level age effects.',evidence='Female and male panels compare median and interquartile range across 204 within-country time-series correlations by age.',archetype='quantitative grid',backend='R',width_mm=183,height_mm=92,formats='PNG/PDF/SVG/600dpi LZW TIFF',uncertainty='Bands show between-country IQR, not confidence intervals. No inferential p-values.',source='age_time_summary.csv and age_withincountry_time_correlations.csv')
write.csv(contract,file.path(out,'figure_S3_contract.csv'),row.names=FALSE)
a <- read.csv(file.path(out,'age_time_summary.csv'))
z <- bind_rows(transmute(a,sex,age_start,method='Rates in levels',median=median_raw_r,lower=q25_raw_r,upper=q75_raw_r),transmute(a,sex,age_start,method='Changes in log rates',median=median_dlog_r,lower=q25_dlog_r,upper=q75_dlog_r))
z$method <- factor(z$method,levels=c('Rates in levels','Changes in log rates'))
write.csv(z,file.path(out,'figure_S3_source_data.csv'),row.names=FALSE)
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),panel.grid=element_blank(),plot.title=element_text(size=9,face='bold'),legend.title=element_blank(),legend.text=element_text(size=8),plot.tag=element_text(size=10,face='bold')))
pal <- c('Rates in levels'='#285F83','Changes in log rates'='#C77F3F')
one <- function(s,label){
 ggplot(filter(z,sex==s),aes(age_start,median,colour=method,fill=method,group=method))+
 geom_hline(yintercept=0,linewidth=.3,colour='grey70')+
 geom_vline(xintercept=50,linewidth=.3,linetype='dotted',colour='grey60')+
 geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.13,colour=NA)+
 geom_line(linewidth=.65)+geom_point(size=1.2)+
 scale_colour_manual(values=pal)+scale_fill_manual(values=pal)+
 scale_x_continuous(breaks=c(5,20,35,50,65,80,95),labels=c('5','20','35','50','65','80','95+'),limits=c(5,95),expand=expansion(mult=.02))+
 scale_y_continuous(breaks=seq(-1,1,.5),limits=c(-1,1),expand=expansion(mult=.02))+
 labs(title=label,x='Age band (starting age, years)',y='Within-country Pearson correlation')
}
p <- (one('female','Female')+one('male','Male'))+plot_layout(guides='collect')+plot_annotation(tag_levels='A',caption='Median and IQR across 204 countries; 30 annual levels or 29 annual log changes, 1992-2021.',theme=theme(plot.caption=element_text(size=7,hjust=0),plot.margin=margin(4,4,4,4))) & theme(legend.position='bottom')
w <- 183/25.4;h <- 92/25.4;stem <- file.path(figure_dir,'Figure_S2')
ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='in',res=180);print(p);dev.off()
ragg::agg_tiff(paste0(stem,'.tif'),width=w,height=h,units='in',res=600,compression='lzw');print(p);dev.off()
cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,family='Arial');print(p);dev.off()
svglite::svglite(paste0(stem,'.svg'),width=w,height=h);print(p);dev.off()
writeLines(c('Figure S2. Age-specific distributions of within-country CKD-DD time-series correlations, 1992-2021.','For each country, sex and age group, Pearson correlations were computed across 30 annual incidence-rate estimates and across 29 first differences of log rates. Lines show the median and shaded bands the interquartile range across 204 countries and territories. The dotted vertical line marks age 50. These are descriptive distributions, not confidence intervals or tests of an individual-level age effect. Estimates do not propagate GBD estimation uncertainty; serial dependence and shared model structure limit inferential interpretation.'),file.path(out,'Figure_S2_caption.txt'))
