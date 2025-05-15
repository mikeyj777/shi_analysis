'''

This will generate inverse (1/x) values as well as log values for each numeric data element.

In addition, the percentile of each value is stored as well as the normalized ([x-mu] / sig) score.

'''

from math import log10, pi
import pandas as pd
from scipy.stats import zscore

df = pd.read_csv('raw_data_inhalation_non_catastrophic.csv')
sparse_range_columns = ['elev_m', 'release_orientation_rad']

def get_data_from_pandas(elem):
    if isinstance(elem, pd.Series):
        elem = elem.values[0]
    return elem

def is_outlier(k_stats, v):
    k_mean = k_stats['mean']
    k_std = k_stats['std']
    return abs(v - k_mean) > 3*k_std

def processing(k, v, k_stats, operation):
    if k not in sparse_range_columns:  # columns where there are only a smaller range of inputs.  outlier analysis not applicable
        if is_outlier(k_stats=k_stats, v=v):
            return None
    ans = None
    try:
        if operation == 'inv':
            ans = 2 # for elevations of zero, inverse will be stored as "2".  larger than inverse elevation of 1, but not nearing infinity
            ans = 1/v
        if operation == 'log':
            ans = -10
            ans = log10(v)
        if operation == 'pct':
            k_min = k_stats['min']
            k_max = k_stats['max']
            ans = 1
            ans = (v - k_min) / (k_max - k_min)
        if operation == 'norm':
            k_mean = k_stats['mean']
            k_std = k_stats['std']
            ans = 1
            ans = (v - k_mean) / k_std
    except:
        pass
    
    return ans

output = []
df_stats = df.describe()
for _, row in df.iterrows():
    row_for_iter = row.to_dict()
    row_out = row.to_dict()
    for k,v in row_for_iter.items():
        if k not in df_stats:
            continue
        if k == 'study_id':
            continue
        if not (isinstance(v, float) or isinstance(v, int)):
            continue
        if pd.isna(v):
            continue
        k_stats_dict = df_stats[k].to_dict()
        for operation in ['inv', 'log']:
            k_inner = f'{operation}_{k}'
            row_out[k_inner] = processing(k=k, v=v, k_stats=k_stats_dict, operation=operation)
            for oper_ranker in ['pct', 'norm']:
                val = row_out[k_inner]
                row_out[f'{oper_ranker}_{k_inner}'] = processing(k = k_inner, v = val, k_stats=k_stats_dict, operation=oper_ranker)

    output.append(row_out)

df_out = pd.DataFrame(output)
df_out.to_csv('results_complete_take_6.csv', index=False)

apple = 1