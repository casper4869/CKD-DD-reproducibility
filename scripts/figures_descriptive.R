suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(maps)})
source('scripts/config.R');a<-core_dir;o<-figure_dir
read<-function(n)read.csv(file.path(a,paste0(n,'.csv')))
qdat<-read('S2_quartile_country_year');geo<-read('S11_spatial_inputs'); mp<-read.csv(file.path(reference_dir,'map_names.csv')); world<-map_data('world')%>%filter(region!='Antarctica')
pal<-c('Same quartile'='#7A6CA6','Higher CKD quartile'='#3B6E9A','Higher DD quartile'='#C9854C'); sigpal<-c('Both predictive tests'='#8C6BB1','CKD to DD only'='#3B6E9A','DD to CKD only'='#C9854C','Neither test'='#C9CDD1')
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),legend.title=element_text(size=8),legend.text=element_text(size=7),plot.title=element_text(size=9,face='bold'),strip.text=element_text(size=8,face='bold'),legend.position='bottom',plot.tag=element_text(size=11,face='bold')))
saveplot<-function(p,name,w=183,h=150){
 svglite::svglite(file.path(o,paste0(name,'.svg')),width=w/25.4,height=h/25.4);print(p);dev.off()
 cairo_pdf(file.path(o,paste0(name,'.pdf')),width=w/25.4,height=h/25.4,family='Arial');print(p);dev.off()
 ragg::agg_tiff(file.path(o,paste0(name,'.tif')),width=w,height=h,units='mm',res=600,compression='lzw');print(p);dev.off()
 ragg::agg_png(file.path(o,paste0(name,'.png')),width=w,height=h,units='mm',res=300);print(p);dev.off()
}
pointmap<-function(d,var,label,limits=NULL){ggplot()+geom_polygon(data=world,aes(long,lat,group=group),fill='#F0F1F2',colour='white',linewidth=.12)+geom_point(data=d,aes(x=lng,y=lat,colour=.data[[var]]),size=1.6,alpha=.9)+scale_colour_gradient(low='#DBE9EF',high='#225A7A',name=label,limits=limits)+coord_quickmap(xlim=c(-180,180),ylim=c(-60,85),expand=FALSE)+theme_void(base_size=8)+theme(legend.position='bottom',legend.key.width=grid::unit(11,'mm'),legend.title=element_text(size=7),legend.text=element_text(size=7),plot.title=element_text(size=9,face='bold'))}
d<-filter(qdat,year==2021,sex_name!='Both')%>%left_join(geo%>%select(location_id,lng,lat),by='location_id')
p<-list(); for(dis in c('CKD','DD')) for(s in c('Female','Male'))p[[length(p)+1]]<-pointmap(filter(d,sex_name==s),dis,paste(dis,'ASIR'))+ggtitle(paste(s,dis))
for(s in c('Female','Male')){z<-filter(d,sex_name==s); p[[length(p)+1]]<-ggplot(z,aes(CKD,DD,colour=type))+geom_point(size=1.5,alpha=.8)+geom_vline(xintercept=quantile(z$CKD,.75),linetype=2,linewidth=.3,colour='grey50')+geom_hline(yintercept=quantile(z$DD,.75),linetype=2,linewidth=.3,colour='grey50')+scale_colour_manual(values=pal,name=NULL)+theme(legend.position='none')+scale_x_log10()+scale_y_log10()+labs(x='CKD ASIR per 100,000',y='DD ASIR per 100,000',title=paste(s,'continuous distribution'))}
saveplot(wrap_plots(p,ncol=2)+plot_annotation(tag_levels='A'),'Figure_1',183,220)
ag<-read('S7_age_correlations');ag$age<-as.numeric(sub('^([0-9]+).*','\\1',ag$age_name));pr<-read('S6_age_profile_descriptive')%>%left_join(geo%>%select(location_id,lng,lat),by='location_id');p<-list()
for(s in c('Female','Male')){z<-filter(pr,sex_name==s);p[[length(p)+1]]<-pointmap(z,'profile_r','Age-profile Pearson r',c(-1,1))+ggtitle(paste(s,'age profiles'))+scale_colour_gradient2(low='#3B6E9A',mid='white',high='#B45C58',midpoint=0,limits=c(-1,1),name='Age-profile r')}
for(s in c('Female','Male'))p[[length(p)+1]]<-ggplot(filter(ag,sex_name==s),aes(age,rho))+geom_hline(yintercept=0,colour='grey75',linewidth=.3)+geom_ribbon(aes(ymin=lower,ymax=upper),fill='#D9E7EF')+geom_line(colour='#3B6E9A',linewidth=.6)+geom_point(aes(shape=q_BH<.05),size=1.7)+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),labels=c('No','Yes'),name='BH q < 0.05')+scale_x_continuous(breaks=c(5,20,35,50,65,80,95),labels=c('5-9','20-24','35-39','50-54','65-69','80-84','95+'))+coord_cartesian(ylim=c(-1,1))+labs(x='Age group (years)',y='Across-country Spearman rho',title=s)
saveplot(wrap_plots(p,ncol=2)+plot_annotation(tag_levels='A'),'Figure_2',183,150)
z<-read('S3_quartile_counts')%>%filter(sex_name!='Both');p<-ggplot(z,aes(year,n,colour=type))+geom_line(linewidth=.7)+scale_colour_manual(values=pal,name=NULL)+facet_wrap(~sex_name,ncol=2)+scale_x_continuous(breaks=c(1992,2000,2010,2021))+labs(x='Year',y='Countries and territories (n)')
saveplot(p,'Figure_3',183,85)
