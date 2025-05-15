from math import log10, pi
import pandas as pd

df = pd.read_csv('results_complete_take_5.csv')

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
        if pd.isna(v):
            continue
        ans = None
        try:
            ans = 1/v
        except:
            pass
        row_out[f'inv_{k}'] = ans

    output.append(row_out)

df_out = pd.DataFrame(output)
df_out.to_csv('results_complete_take_6.csv', index=False)

apple = 1