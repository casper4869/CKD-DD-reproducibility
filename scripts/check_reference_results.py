"""Compare locally recomputed aggregate outputs with the locked reference results.

No data access or network calls. Run after all documented analysis steps.
Missing results fail, rather than silently passing a partial reproduction.
"""
import csv
import json
import math
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULT=Path(os.environ.get('CKDDD_RESULTS', ROOT/'results'))
MAP={
    'S4_threshold_sensitivity.csv':'core/S4_threshold_sensitivity.csv',
    'S8_age_contrast.csv':'core/S8_age_contrast.csv',
    'A2_model_summary.csv':'temporal/A2_model_summary.csv',
    'A18_residual_diagnostics_summary.csv':'temporal/A18_residual_diagnostics_summary.csv',
    'A20_matched_model_pre_pandemic_summary.csv':'temporal/A20_matched_model_pre_pandemic_summary.csv',
    'A30_original_model_short_validation.csv':'temporal/A30_original_model_short_validation.csv',
    'A31_original_model_nine_year_endpoint.csv':'temporal/A31_original_model_nine_year_endpoint.csv',
    'GWR_sensitivity_diagnostics.csv':'spatial/GWR_sensitivity_diagnostics.csv',
    'GWR_coordinate_agreement.csv':'spatial/GWR_coordinate_agreement.csv',
    'S3_age_reference_summary.csv':'age_reference/S3_age_reference_summary.csv',
}


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    return rows


def same_cell(a,b):
    if a==b:return True
    if a in ('','NA','NaN') and b in ('','NA','NaN'):return True
    try:
        aa,bb=float(a),float(b)
        if math.isnan(aa) and math.isnan(bb):return True
        return math.isclose(aa,bb,rel_tol=1e-7,abs_tol=1e-9)
    except ValueError:return False


def compare(expected,actual):
    a,b=read(expected),read(actual)
    if len(a)!=len(b):return {'pass':False,'reason':'Row count differs','expected_rows':len(a),'actual_rows':len(b)}
    if not a or list(a[0])!=list(b[0]):return {'pass':False,'reason':'Column order/schema differs'}
    # All files use deterministic order from the executed scripts; compare every cell.
    different=[{'row':i+2,'column':k} for i,(x,y) in enumerate(zip(a,b)) for k in x if not same_cell(x[k],y[k])]
    return {'pass':not different,'rows':len(a),'columns':len(a[0]),'different_cells':different[:20],'n_different':len(different)}


def main():
    checks={}
    for name,target in MAP.items():
        p=RESULT/target
        checks[name]=compare(ROOT/'reference_results'/name,p) if p.exists() else {'pass':False,'reason':'Recomputed file missing'}
    result={'status':'pass' if all(x['pass'] for x in checks.values()) else 'fail','relative_tolerance':1e-7,'absolute_tolerance':1e-9,'checks':checks}
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['status']=='pass' else 1)


if __name__=='__main__':main()
