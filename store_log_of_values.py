from math import log10, pi
import pandas as pd

df = pd.read_csv('results_complete_take_3.csv')

def get_data_from_pandas(elem):
    if isinstance(elem, pd.Series):
        elem = elem.values[0]
    return elem

output = []
for _, row in df.iterrows():
    row_for_iter = row.to_dict()
    row_out = row.to_dict()
    for k,v in row_for_iter.items():
        if k == 'study_id':
            continue
        if not (isinstance(v, float) or isinstance(v, int)):
            continue
        log_val = -10
        if v > 0:
            log_val = log10(v)
        row_out[f'log_{k}'] = log_val

    output.append(row_out)

df_out = pd.DataFrame(output)
df_out.to_csv('results_complete_take_4.csv', index=False)

apple = 1