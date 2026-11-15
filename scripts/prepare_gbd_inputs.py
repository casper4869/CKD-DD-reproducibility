"""Filter authorized official IHME CSV exports locally; never downloads data."""
import argparse
from pathlib import Path
import pandas as pd

AGES=[f'{a}-{a+4} years' for a in range(5,95,5)]+['95+ years']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ckd',type=Path,required=True,help='Official CKD export containing all requested years and ages')
    p.add_argument('--dd',type=Path,required=True,help='Official DD export containing all requested years and ages')
    p.add_argument('--output',type=Path,default=Path('data/private'))
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1]
    ids=set(pd.read_csv(root/'data/reference/geometry_reference_points.csv').location_id)
    extracts={}
    for dis,cause,source in [('CKD','Chronic kidney disease',a.ckd),('DD','Depressive disorders',a.dd)]:
        parts=[]
        for d in pd.read_csv(source,chunksize=300000):
            keep=(d.measure_name=='Incidence')&(d.metric_name=='Rate')&(d.cause_name==cause)&d.location_id.isin(ids)&d.year.between(1992,2021)
            parts.append(d.loc[keep])
        d=pd.concat(parts,ignore_index=True)
        assert not d.duplicated(['location_id','sex_name','year','age_name']).any()
        asir=d[(d.age_name=='Age-standardized')&d.sex_name.isin(['Both','Female','Male'])]
        assert len(asir)==204*3*30
        age=d[(d.year==2021)&d.sex_name.isin(['Female','Male'])]
        assert len(age[age.age_name.isin(AGES)])==204*2*19
        asir.to_csv(a.output/f'{dis}_ASIR_1992_2021.csv',index=False)
        age.to_csv(a.output/f'{dis}_age_2021.csv',index=False)
        extracts[dis]=d
    for sex in ['Female','Male']:
        left=extracts['CKD'].loc[lambda d:(d.sex_name==sex)&d.age_name.isin(AGES),['location_name','age_name','year','val']].rename(columns={'val':'ckd_val'})
        right=extracts['DD'].loc[lambda d:(d.sex_name==sex)&d.age_name.isin(AGES),['location_name','age_name','year','val']].rename(columns={'val':'dd_val'})
        z=left.merge(right,on=['location_name','age_name','year'],validate='one_to_one')
        assert len(z)==204*19*30
        z.to_csv(a.output/f'CKD_DD_{sex.lower()}_LONG_1992_2021.csv',index=False)
    print('Local extraction complete. SDI, population and PM2.5 matrices must be supplied separately; run validate_inputs.py.')


if __name__=='__main__':main()
