"""Validate authorized local input CSVs without network access (standard library)."""
import csv
import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = Path(os.environ.get('CKDDD_PRIVATE_DATA', ROOT/'data/private'))
AGES = [f'{a}-{a+4} years' for a in range(5,95,5)] + ['95+ years']


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def check_inputs(private=PRIVATE):
    refs = read(ROOT/'data/reference/geometry_reference_points.csv')
    expected_ids = {int(r['location_id']) for r in refs}
    expected_names = {r['location_name'] for r in refs}
    for disease, cause in [('CKD','Chronic kidney disease'),('DD','Depressive disorders')]:
        rows=read(private/f'{disease}_ASIR_1992_2021.csv')
        keys={(int(r['location_id']),r['sex_name'],int(r['year'])) for r in rows}
        assert len(rows)==len(keys)==204*3*30
        assert {k[0] for k in keys}==expected_ids
        assert {k[1] for k in keys}=={'Both','Female','Male'}
        assert {k[2] for k in keys}==set(range(1992,2022))
        assert all(r['cause_name']==cause and r['age_name']=='Age-standardized' and r['measure_name']=='Incidence' and r['metric_name']=='Rate' and math.isfinite(float(r['val'])) and float(r['val'])>0 for r in rows)
        rows=[r for r in read(private/f'{disease}_age_2021.csv') if r['age_name'] in AGES]
        keys={(int(r['location_id']),r['sex_name'],r['age_name']) for r in rows}
        assert len(rows)==len(keys)==204*2*19
        assert {k[0] for k in keys}==expected_ids and {k[1] for k in keys}=={'Female','Male'}
        assert all(int(r['year'])==2021 and r['cause_name']==cause and r['metric_name']=='Rate' and r['measure_name']=='Incidence' and math.isfinite(float(r['val'])) and float(r['val'])>=0 for r in rows)
    for name in ['SDI.csv','Population.csv','PM25.csv']:
        rows=read(private/name)
        assert len(rows)==204 and {r['location_name'] for r in rows}==expected_names
        assert all(math.isfinite(float(r[f'val_{y}'])) for r in rows for y in range(1992,2022))
    for sex in ['female','male']:
        rows=[r for r in read(private/f'CKD_DD_{sex}_LONG_1992_2021.csv') if r['age_name'] in AGES]
        keys={(r['location_name'],r['age_name'],int(r['year'])) for r in rows}
        assert len(rows)==len(keys)==204*19*30
        assert {k[0] for k in keys}==expected_names
        assert {k[2] for k in keys}==set(range(1992,2022))
        assert all(math.isfinite(float(r[k])) and float(r[k])>0 for r in rows for k in ['ckd_val','dd_val'])
    return {'status':'pass','n_locations':204,'n_years':30,'n_age_groups':19,'input_directory':str(private)}


if __name__=='__main__':
    print(check_inputs())
