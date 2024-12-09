import pandas as pd
from apicalls.api_oop import GOClient, UniProtClient

df_core = pd.read_csv("./prlists/all_col_core_df.csv")
go_client = GOClient()
uniprot_client = UniProtClient()

uniprot_client.get_pr_fnc(df_core.iloc[0]['UniProtAC'])

fncs = []
for idx, row in df_core.iterrows():
    print(row.UniProtAC)
    resp = uniprot_client.get_pr_fnc(row.UniProtAC)
    fncs.append(resp)
df_core['pr_funcs'] = fncs


molecular_function = go_client.get_relative_terms('GO:0003674')
mf_childs = {}
for mf in molecular_function['childrens']:
    print(str(mf))
    resp = go_client.get_relative_terms(str(mf))
    mf_childs[resp['name']] = resp['childrens']

fncs = df_core['pr_funcs'].tolist()
mol_functions = []
fncs[0]
for pr in fncs:
    for k, v in mf_childs.items():
        intersection = list(set(pr) & set(v))
        print(intersection)

mf_childs.values()
df_core.head()
df_core['molecular_function'] = mol_functions
df_core.to_csv("./prlists/all_col_core_df.csv", index=False)
