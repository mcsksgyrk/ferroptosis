import pandas as pd
from apicalls.api_oop import GOClient, UniProtClient
from config import OUTPUTS_DIR

df_core = pd.read_csv(OUTPUTS_DIR / "core_geneProducts.csv")
go_client = GOClient()
uniprot_client = UniProtClient()
df_core

fncs = []
for idx, row in df_core.iterrows():
    resp = uniprot_client.get_pr_fnc(row.uniprot_id)
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

go_to_function = {}
for function_name, go_ids in mf_childs.items():
    for go_id in go_ids:
        go_to_function[go_id] = function_name

results = []
for _, row in df_core.iterrows():
    protein_functions = row.pr_funcs  # List of GO IDs for this protein
    matched_categories = {}

    for go_id in protein_functions:
        if go_id in go_to_function:
            category = go_to_function[go_id]
            if category in matched_categories:
                matched_categories[category].append(go_id)
            else:
                matched_categories[category] = [go_id]
    results.append(matched_categories)

results
mf_childs
df_core.head()
df_core['molecular_function'] = mol_functions
