from pathlib import Path
from itertools import combinations
import hashlib
import json
import numpy as np
import pandas as pd

import os
base=Path(os.environ.get('CKDDD_RESULTS','results'))
out=base/'evidence_spatial_age'
out.mkdir(parents=True,exist_ok=True)
sources={
    'age_profiles':base/'core/S6_age_profile_descriptive.csv',
    'quartile_panel':base/'core/S2_quartile_country_year.csv',
    'local_coefficients':base/'spatial/GWR_sensitivity_coefficients.csv',
}
age=pd.read_csv(sources['age_profiles'])
assert len(age)==408 and not age.duplicated(['location_id','sex_name']).any()
assert set(age.n_age)=={19}
rows=[]
for sex,a in age.groupby('sex_name'):
    r=a.profile_r
    rows.append(dict(sex=sex,n_country=len(a),n_age=19,positive=int((r>0).sum()),negative=int((r<0).sum()),zero=int((r==0).sum()),positive_percent=100*(r>0).mean(),negative_percent=100*(r<0).mean(),median_r=r.median(),q25_r=r.quantile(.25),q75_r=r.quantile(.75),min_r=r.min(),max_r=r.max()))
pd.DataFrame(rows).to_csv(out/'P1_age_profile_summary.csv',index=False)
ap=age.pivot(index=['location_id','location_name'],columns='sex_name',values='profile_r').reset_index()
ap['positive_female']=ap.Female>0;ap['positive_male']=ap.Male>0
ap.to_csv(out/'P2_country_age_profiles_both_sexes.csv',index=False)

qd=pd.read_csv(sources['quartile_panel'])
hh=qd[(qd.year==2021)&(qd.sex_name=='Both')&(qd.q_CKD==4)&(qd.q_DD==4)].copy().sort_values('location_id')
assert len(hh)==10
hist=qd[(qd.sex_name=='Both')&(qd.q_CKD==4)&(qd.q_DD==4)].groupby('location_id').year.agg(['size','min','max']).rename(columns={'size':'number_high_high_years_1992_2021','min':'first_high_high_year','max':'last_high_high_year'})
hh=hh.merge(hist,on='location_id',how='left')
for sex in ['Female','Male']:
    a=qd[(qd.sex_name==sex)&(qd.year==2021)][['location_id','q_CKD','q_DD']].rename(columns={'q_CKD':f'{sex}_q_CKD','q_DD':f'{sex}_q_DD'})
    hh=hh.merge(a,on='location_id',how='left')
hh.to_csv(out/'P3_high_high_2021_all_ten_locations.csv',index=False)

sp=pd.read_csv(sources['local_coefficients'])
assert len(sp)==204*2*2*2
assert not sp.duplicated(['location_id','target','coordinate','spec']).any()
sp['condition']=sp.coordinate+'__'+sp.spec
family_cols=['q_BH_204','q_BY_204','q_BH_408','q_BY_408']
rows=[]
for (target,coord,spec),a in sp.groupby(['target','coordinate','spec']):
    assert len(a)==204
    d=dict(target=target,coordinate=coord,spec=spec,n_country=len(a),positive=int((a.coefficient>0).sum()),negative=int((a.coefficient<0).sum()),zero=int((a.coefficient==0).sum()),median_coefficient=a.coefficient.median(),q25_coefficient=a.coefficient.quantile(.25),q75_coefficient=a.coefficient.quantile(.75))
    for q in family_cols:
        d[q+'_positive']=int(((a[q]<.05)&(a.coefficient>0)).sum())
        d[q+'_negative']=int(((a[q]<.05)&(a.coefficient<0)).sum())
        d[q+'_total']=int((a[q]<.05).sum())
    rows.append(d)
pd.DataFrame(rows).to_csv(out/'P4_GWR_all_conditions_symmetric_signal_counts.csv',index=False)

stable=[];country=[];pairwise=[]
for target,a in sp.groupby('target'):
    co=a.pivot(index=['location_id','location_name'],columns='condition',values='coefficient').sort_index()
    assert co.shape==(204,4) and not co.isna().any().any()
    cp=co.gt(0).all(axis=1);cn=co.lt(0).all(axis=1)
    d=dict(target=target,n_country=len(co),n_fixed_conditions=4,stable_positive=int(cp.sum()),stable_negative=int(cn.sum()),mixed_sign=int((~(cp|cn)).sum()),stable_positive_percent=100*cp.mean(),stable_negative_percent=100*cn.mean())
    oc=co.reset_index()
    oc['stable_positive_all4']=cp.to_numpy();oc['stable_negative_all4']=cn.to_numpy();oc['mixed_sign_all4']=~(cp|cn).to_numpy()
    for q in family_cols:
        qp=a.pivot(index=['location_id','location_name'],columns='condition',values=q).reindex(index=co.index,columns=co.columns)
        reject_all=qp.lt(.05).all(axis=1)
        reject_any=qp.lt(.05).any(axis=1)
        d[q+'_stable_positive_reject_all4']=int((cp&reject_all).sum())
        d[q+'_stable_negative_reject_all4']=int((cn&reject_all).sum())
        d[q+'_reject_any4_positive_any']=int((qp.lt(.05)&co.gt(0)).any(axis=1).sum())
        d[q+'_reject_any4_negative_any']=int((qp.lt(.05)&co.lt(0)).any(axis=1).sum())
        oc[q+'_reject_all4']=reject_all.to_numpy()
        oc[q+'_reject_any4']=reject_any.to_numpy()
    oc.insert(0,'target',target);country.append(oc);stable.append(d)
    for x,y in combinations(co.columns,2):
        pairwise.append(dict(target=target,condition_1=x,condition_2=y,n_country=204,same_sign=int((np.sign(co[x])==np.sign(co[y])).sum()),both_positive=int(((co[x]>0)&(co[y]>0)).sum()),both_negative=int(((co[x]<0)&(co[y]<0)).sum()),coefficient_Pearson_r=co[x].corr(co[y])))
pd.DataFrame(stable).to_csv(out/'P5_GWR_stable_signals_across_all_four_conditions.csv',index=False)
pd.concat(country,ignore_index=True).to_csv(out/'P6_GWR_country_all_condition_stability.csv',index=False)
pd.DataFrame(pairwise).to_csv(out/'P7_GWR_pairwise_condition_agreement.csv',index=False)

manifest={k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for k,p in sources.items()}
meta=dict(mode='Read-only summary of existing fixed results; no new fitted model, selection threshold, geography, covariate, or significance family.',source_manifest=manifest,fixed_grid=['original city/reference points vs documented geometry points','SDI+PM2.5 vs SDI-only','CKD-outcome and DD-outcome separately'],significance_threshold='Previously fixed q<0.05; all four saved BH/BY adjustment columns are reported symmetrically.',age_note='Descriptive Pearson correlation across 19 age groups within each country in 2021; no patient-level or longitudinal age effect.',spatial_note='Conditional local coefficient tests remain approximate; sign persistence is a descriptive robustness summary, not independent replication or causal validation. The all-four intersection is not a new conjunction test, and no separately calibrated FDR is claimed for this intersection.',high_high_note='Descriptive simultaneous upper-quartile membership; not excess overlap or individual comorbidity.',number_age_profile_positive_both_sexes=int((ap.positive_female&ap.positive_male).sum()),number_age_profile_negative_both_sexes=int(((~ap.positive_female)&(~ap.positive_male)).sum()))
(out/'audit_definition_and_input_hashes.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(stable).to_string(index=False))
print(pd.read_csv(out/'P1_age_profile_summary.csv').to_string(index=False))
print(hh[['location_name','CKD','DD','number_high_high_years_1992_2021']].to_string(index=False))

# Compact views of the same symmetric summaries; no additional model/test.
signals=pd.read_csv(out/'P4_GWR_all_conditions_symmetric_signal_counts.csv')
signals[['target','coordinate','spec','n_country','positive','negative','q_BH_204_positive','q_BH_204_negative','q_BY_408_positive','q_BY_408_negative']].to_csv(out/'P8_compact_paper_table_all_fixed_settings.csv',index=False)
stable=pd.read_csv(out/'P5_GWR_stable_signals_across_all_four_conditions.csv')
stable[['target','n_country','n_fixed_conditions','stable_positive','stable_negative','mixed_sign','q_BH_204_stable_positive_reject_all4','q_BH_204_stable_negative_reject_all4','q_BY_204_stable_positive_reject_all4','q_BY_204_stable_negative_reject_all4','q_BY_408_stable_positive_reject_all4','q_BY_408_stable_negative_reject_all4']].to_csv(out/'P9_compact_stability_summary_all_four_settings.csv',index=False)
