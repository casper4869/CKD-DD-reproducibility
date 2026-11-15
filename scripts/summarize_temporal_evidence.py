"""Read-only re-summary of already fitted CKD-DD results. No refits or new thresholds.

All comparisons retain the complete 204-country panel. No resampling, new hypothesis
tests, tuning, alpha changes, or favorable-origin selection are performed.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd

import os
BASE=Path(os.environ.get('CKDDD_RESULTS','results'))
AUD=BASE/'temporal'
OUT=BASE/'evidence_temporal'
OUT.mkdir(parents=True,exist_ok=True)
ALL_ORIGINS=[2012,2015,2018]
KEY=['disease','location_id','origin','h']
sources=[AUD/'A7_forecast_validation.csv',AUD/'A28_original_ARIMAX_validation.csv',
         AUD/'A38_raw_ARIMAX_validation.csv',BASE/'core'/'S15_rolling_validation.csv',
         AUD/'A1_model_tests_country.csv',AUD/'A19_matched_model_pre_pandemic.csv',
         BASE/'core'/'S9_predictive_tests_diagnostics.csv']

a7=pd.read_csv(sources[0])
f=pd.concat([a7,pd.read_csv(sources[1]),pd.read_csv(sources[2])],ignore_index=True)
v1=pd.read_csv(sources[3])
v1=v1[v1.model!='Naive'].copy()
v1['model']=v1.model.map({'VARX':'Differenced_VARX_BIC','ARIMAX':'ARIMAX_log_max2'})
f=pd.concat([f,v1],ignore_index=True)
assert not f.duplicated(KEY+['model']).any()
naive=f[f.model=='Naive'][KEY+['APE','median','actual']].rename(columns={'APE':'APE_Naive','median':'point_Naive','actual':'actual_Naive'})
paired=f[f.model!='Naive'].merge(naive,on=KEY,how='left',validate='many_to_one')
assert paired.APE_Naive.notna().all()
assert np.allclose(paired.actual,paired.actual_Naive)
paired['delta_APE_pp']=paired.APE-paired.APE_Naive
paired['square_error']=(paired['median']-paired.actual)**2
paired['square_error_Naive']=(paired.point_Naive-paired.actual)**2
location_names=pd.read_csv(AUD/'A1_model_tests_country.csv').query("model=='original_raw_endo'").set_index('location_id').location_name
paired['location_name']=paired.location_id.map(location_names)
short=paired[paired.origin.isin(ALL_ORIGINS)&(paired.h<=3)].copy()
assert all(short.groupby(['model','disease','origin','h']).size()==204)
short.to_csv(OUT/'P01_short_horizon_paired_predictions.csv',index=False)

def summarize(z):
    d=z.delta_APE_pp
    return pd.Series({'n_predictions':len(z),'n_countries':z.location_id.nunique(),
        'mean_APE_model':z.APE.mean(),'mean_APE_Naive':z.APE_Naive.mean(),
        'mean_APE_relative_reduction_pct':100*(1-z.APE.mean()/z.APE_Naive.mean()),
        'median_APE_model':z.APE.median(),'median_APE_Naive':z.APE_Naive.median(),
        'paired_mean_delta_APE_pp':d.mean(),'paired_median_delta_APE_pp':d.median(),
        'paired_delta_Q1_pp':d.quantile(.25),'paired_delta_Q3_pp':d.quantile(.75),
        'predictions_lower_APE':int((d<0).sum()),'predictions_higher_APE':int((d>0).sum()),'predictions_tied_APE':int((d==0).sum()),
        'RMSE_model':np.sqrt(z.square_error.mean()),'RMSE_Naive':np.sqrt(z.square_error_Naive.mean()),
        'coverage95_model':z.covered.mean()})

short.groupby(['model','disease'],sort=True).apply(summarize,include_groups=False).reset_index().to_csv(OUT/'P02_short_horizon_all_models_summary.csv',index=False)
short.groupby(['model','disease','origin','h'],sort=True).apply(summarize,include_groups=False).reset_index().to_csv(OUT/'P03_short_horizon_all_origin_horizon_cells.csv',index=False)
short.groupby(['model','disease','origin'],sort=True).apply(summarize,include_groups=False).reset_index().to_csv(OUT/'P04_short_horizon_each_origin.csv',index=False)
short.groupby(['model','disease','h'],sort=True).apply(summarize,include_groups=False).reset_index().to_csv(OUT/'P05_short_horizon_each_step.csv',index=False)

def country_metrics(z):
    return pd.Series({'n_evaluations':len(z),'mean_APE_model':z.APE.mean(),'mean_APE_Naive':z.APE_Naive.mean(),
        'delta_mean_APE_pp':z.APE.mean()-z.APE_Naive.mean(),
        'median_APE_model':z.APE.median(),'median_APE_Naive':z.APE_Naive.median(),
        'delta_median_APE_pp':z.APE.median()-z.APE_Naive.median(),
        'RMSE_model':np.sqrt(z.square_error.mean()),'RMSE_Naive':np.sqrt(z.square_error_Naive.mean())})

c=short.groupby(['model','disease','location_id'],sort=True).apply(country_metrics,include_groups=False).reset_index()
c['location_name']=c.location_id.map(location_names)
assert all(c.n_evaluations==9)
c.to_csv(OUT/'P06_country_paired_short_metrics.csv',index=False)
co=short.groupby(['model','disease','location_id','origin'],sort=True).apply(country_metrics,include_groups=False).reset_index()
co['location_name']=co.location_id.map(location_names)
co.to_csv(OUT/'P07_country_origin_paired_metrics.csv',index=False)

cs=[]
for (model,dis),d in c.groupby(['model','disease']):
    z=co[(co.model==model)&(co.disease==dis)]
    wins=z.assign(win=z.delta_mean_APE_pp<0).groupby('location_id').win.sum()
    cells=short[(short.model==model)&(short.disease==dis)].groupby(['origin','h']).delta_APE_pp.mean()
    cs.append({'model':model,'disease':dis,'n_countries':len(d),
        'countries_lower_9evaluation_meanAPE':int((d.delta_mean_APE_pp<0).sum()),
        'countries_higher_9evaluation_meanAPE':int((d.delta_mean_APE_pp>0).sum()),
        'countries_lower_9evaluation_medianAPE':int((d.delta_median_APE_pp<0).sum()),
        'countries_lower_RMSE':int((d.RMSE_model<d.RMSE_Naive).sum()),
        'median_country_delta_meanAPE_pp':d.delta_mean_APE_pp.median(),
        'Q1_country_delta_meanAPE_pp':d.delta_mean_APE_pp.quantile(.25),
        'Q3_country_delta_meanAPE_pp':d.delta_mean_APE_pp.quantile(.75),
        'countries_meanAPE_better_all_3_origins':int((wins==3).sum()),
        'countries_meanAPE_better_2plus_origins':int((wins>=2).sum()),
        'countries_meanAPE_better_no_origin':int((wins==0).sum()),
        'origin_horizon_cells_lower_meanAPE':int((cells<0).sum()),'origin_horizon_cells_total':len(cells)})
pd.DataFrame(cs).to_csv(OUT/'P08_country_consistency_summary.csv',index=False)

long=paired[(paired.h==9)&paired.origin.isin([2009,2012])].copy()
assert all(long.groupby(['model','disease','origin']).size()==204)
long.to_csv(OUT/'P09_nine_year_endpoint_paired_predictions.csv',index=False)
long.groupby(['model','disease'],sort=True).apply(summarize,include_groups=False).reset_index().to_csv(OUT/'P10_nine_year_endpoint_summary.csv',index=False)
long.groupby(['model','disease','origin'],sort=True).apply(summarize,include_groups=False).reset_index().to_csv(OUT/'P11_nine_year_endpoint_each_origin.csv',index=False)

# Frozen disease-equation rejections; every family is 204x2, alpha=.05 as before.
tests=pd.read_csv(sources[4]);pre=pd.read_csv(sources[5]);pre=pre[pre.ending==2019].copy()
diff=pd.read_csv(sources[6])
status=[]
def flags_frame(d,label,method,c1,c2):
    return pd.DataFrame({'location_id':d.location_id,'location_name':d.location_name,'model':label,'correction':method,
        'q_CKD_DD':d[c1],'q_DD_CKD':d[c2],
        'CKD_DD_reject':d[c1]<.05,'DD_CKD_reject':d[c2]<.05,
        'both_reject':(d[c1]<.05)&(d[c2]<.05)})
for model,d in tests.groupby('model'):
    for label,cols in [('BH408',['q_CKD_DD','q_DD_CKD']),('BY408',['qBY_CKD_DD','qBY_DD_CKD']),
                       ('HC3_BH408',['q_HC3_CKD_DD','q_HC3_DD_CKD']),('HC3_BY408',['qBY_HC3_CKD_DD','qBY_HC3_DD_CKD'])]:
        status.append(flags_frame(d,model,label,*cols))
status.append(flags_frame(pre,'original_raw_endo_ending2019','BH408','q_CKD_DD','q_DD_CKD'))
for spec,d in diff.groupby('spec'):
    for label,cols in [('BH408',['q_CKD_DD','q_DD_CKD']),('BY408',['qBY_CKD_DD','qBY_DD_CKD'])]:
        status.append(flags_frame(d,'revision_v1_'+spec,label,*cols))
status=pd.concat(status,ignore_index=True)
assert all(status.groupby(['model','correction']).size()==204)
status.to_csv(OUT/'P12_all_existing_test_statuses.csv',index=False)
ss=[]
for (model,correction),d in status.groupby(['model','correction']):
    ss.append({'model':model,'correction':correction,'n_countries':len(d),'CKD_DD_reject':int(d.CKD_DD_reject.sum()),
        'DD_CKD_reject':int(d.DD_CKD_reject.sum()),'both_reject':int(d.both_reject.sum()),
        'CKD_only':int((d.CKD_DD_reject&~d.DD_CKD_reject).sum()),
        'DD_only':int((~d.CKD_DD_reject&d.DD_CKD_reject).sum()),
        'neither':int((~d.CKD_DD_reject&~d.DD_CKD_reject).sum())})
pd.DataFrame(ss).to_csv(OUT/'P13_all_existing_test_summary.csv',index=False)

primary=status[(status.model=='original_raw_endo')&(status.correction=='BH408')].set_index('location_id')
cross=primary[['location_name','CKD_DD_reject','DD_CKD_reject','both_reject']].add_prefix('primary_')
ret=[]
for (model,correction),d in status.groupby(['model','correction']):
    z=d.set_index('location_id').loc[primary.index]
    for direct in ['CKD_DD_reject','DD_CKD_reject','both_reject']:
        a=primary[direct];b=z[direct]
        ret.append({'model':model,'correction':correction,'direction':direct,'n_countries':204,
                    'primary_reject':int(a.sum()),'sensitivity_reject':int(b.sum()),
                    'rejected_in_both':int((a&b).sum()),'primary_only':int((a&~b).sum()),
                    'sensitivity_only':int((~a&b).sum()),'neither':int((~a&~b).sum())})
    if (model=='original_raw_endo' and correction!='BH408') or model=='original_raw_endo_ending2019':
        for direct in ['CKD_DD_reject','DD_CKD_reject','both_reject']:
            cross[f'{model}__{correction}__{direct}']=z[direct]
pd.DataFrame(ret).to_csv(OUT/'P14_primary_signal_retention_all_specs.csv',index=False)
cross.reset_index().to_csv(OUT/'P15_primary_BH_BY_HC3_period_country_status.csv',index=False)

intersection=[]
def named_status(model,correction,direct):
    return status[(status.model==model)&(status.correction==correction)].set_index('location_id').loc[primary.index,direct]
sets=[('primary_BH_and_BY',[('original_raw_endo','BH408'),('original_raw_endo','BY408')]),
      ('primary_BH_and_HC3BH',[('original_raw_endo','BH408'),('original_raw_endo','HC3_BH408')]),
      ('primary_BH_and_ending2019BH',[('original_raw_endo','BH408'),('original_raw_endo_ending2019','BH408')]),
      ('primary_BY_and_ending2019BH',[('original_raw_endo','BY408'),('original_raw_endo_ending2019','BH408')]),
      ('primary_BH_BY_HC3BH_ending2019BH',[('original_raw_endo','BH408'),('original_raw_endo','BY408'),('original_raw_endo','HC3_BH408'),('original_raw_endo_ending2019','BH408')]),
      ('primary_BH_BY_HC3BY_ending2019BH',[('original_raw_endo','BH408'),('original_raw_endo','BY408'),('original_raw_endo','HC3_BY408'),('original_raw_endo_ending2019','BH408')])]
for label,conditions in sets:
    for direct in ['CKD_DD_reject','DD_CKD_reject','both_reject']:
        keep=np.logical_and.reduce([named_status(model,corr,direct).to_numpy() for model,corr in conditions])
        intersection.append({'intersection':label,'direction':direct,'n_all_countries':204,'n_retained':int(keep.sum()),
                             'interpretation':'Descriptive overlap of frozen rejections; not a new corrected test or validation certificate'})
pd.DataFrame(intersection).to_csv(OUT/'P16_frozen_signal_intersections.csv',index=False)

manifest={'n_countries':204,'short_origins':ALL_ORIGINS,'short_horizons':[1,2,3],
          'notes':['No model fitting or coefficient updates','No alpha or multiplicity changes','No new hypothesis tests or statistical significance claims for paired errors','Country/origin/horizon comparisons are descriptive; annual and spatial dependence remain'],
          'inputs':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]}
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(pd.DataFrame(cs).query("model=='Original_VARX_p2'").to_string(index=False))
print(pd.DataFrame(intersection).to_string(index=False))
