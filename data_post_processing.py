'''

This will generate inverse (1/x) values as well as log values for each numeric data element.

In addition, the percentile of each value is stored as well as the normalized ([x-mu] / sig) score.

'''

from math import log10, pi
import pandas as pd

df = pd.read_csv('raw_data_inhalation_non_catastrophic.csv')

def get_data_from_pandas(elem):
    if isinstance(elem, pd.Series):
        elem = elem.values[0]
    return elem

def processing(k, v, k_stats, operation):
    ans = None
    try:
        if operation == 'inv':
            ans = 1/v
        if operation == 'log':
            ans = -10
            ans = log10(v)
        if operation == 'pct':
            k_min = k_stats['min']
            k_max = k_stats['max']
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
        k_25_pct = k_stats_dict['25%']
        k_75_pct = k_stats_dict['75%']
        for operation in ['inv', 'log', 'pct', 'norm']:
            row_out[f'{operation}_{k}'] = None
            if v < k_25_pct:    
                continue
            if v > k_75_pct:
                continue
            row_out[f'{operation}_{k}'] = processing(k=k, v=v, k_stats=k_stats_dict, operation=operation)
        

    output.append(row_out)

df_out = pd.DataFrame(output)
df_out.to_csv('results_complete_take_6.csv', index=False)

apple = 1