from math import log10
import pandas as pd

df = pd.read_csv('results_for_log.csv')

def isnumeric(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)

output = []
for _, row in df.iterrows():
    row_for_iter = row.to_dict()
    row_out = row.to_dict()
    for k, v in row_for_iter.items():
        if k.lower() == 'study_id':
            continue
        if k[:3].lower() == 'log':
            continue
        if not isnumeric(v):
            continue
        ans = 0
        k_new = f'log_{k}'
        if k[-5:].lower() == 'deg_c':
            k_new = k_new[:-1] + 'k'
            if v is not None:
                v += 273.15
        if v is not None and v != 0:
            ans = log10(v)
        row_out[k_new] = ans
    output.append(row_out)

df_out = pd.DataFrame(output)
df_out.to_csv('results_complete.csv')
