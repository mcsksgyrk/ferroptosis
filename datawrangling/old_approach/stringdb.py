import pandas as pd
import apicalls.api_oop as api
import json

df_cons = pd.read_csv("./sources/string/9606.protein.links.full.v12.0.txt",
                      delim_whitespace=True)
df_core = pd.read_csv("./prlists/all_col_core_df.csv")
uniprot_client = api.UniProtClient()
df_cons.columns

string_ids = []
for idx, row in df_core.iterrows():
    find = uniprot_client.convert_from_uniprot_id('STRING', [row.UniProtAC], False)
    if find:
        string_ids.append(find)

with open("./prlists/stringIDs.txt", "w") as f:
    json.dump(string_ids, f, indent=4)

dict_list = [tup[0] for tup in string_ids]
string_values = [list(d.values()) for d in dict_list]

hits = df_cons[
    (df_cons.protein1.isin(string_values)) & (df_cons.protein2.isin(string_values))]

clean_dict_list = [d for d in dict_list if d]
id_mapping = {list(d.values())[0]: list(d.keys())[0] for d in clean_dict_list}
hits['protein1'] = hits['protein1'].map(id_mapping).fillna(hits['protein1'])
hits['protein2'] = hits['protein2'].map(id_mapping).fillna(hits['protein2'])
# score minnél nagyobb annál nagyobb a konfidencia 700-tól jó kb
hits.to_csv('./datawrangling/string_cons.csv')
