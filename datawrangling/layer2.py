import apicalls.api_oop as api
import pandas as pd
from database.sqlite_db_api import PsimiSQL


op_interact = pd.read_csv("./sources/omnipath/omnipath_pr_interacts.txt", delim_whitespace=True)
core_df = pd.read_csv('./prlists/all_col_core_df.csv')

uniprot_client = api.UniProtClient()
go_client = api.GOClient()

SQL_SEED = "./database/network_db_seed.sql"
DB_DESTINATION = "./output/omnipath2"
db_api = PsimiSQL(SQL_SEED)

#layer 0, connection between core prs
df = op_interact[
    (op_interact.source.isin(core_df.UniProtAC)) & (op_interact.target.isin(core_df.UniProtAC))]


core_df['molecular_function'] = core_df['molecular_function'].apply(
            lambda x: '|'.join(item.strip().strip("'") for item in x.strip('[]').split(','))
    )
for idx, row in core_df.iterrows():
    aux_dict = dict()
    aux_dict['name'] = row.UniProtAC
    aux_dict['gene_name'] = row.Gene_Name
    aux_dict['tax_id'] = row.tax_id
    aux_dict['pathways'] = row.pathways
    aux_dict['source'] = 'todo'
    cleaned = '|'.join(item.strip().strip("'") for item in row.molecular_function.strip('[]').split(','))
    aux_dict['function'] = cleaned
    db_api.insert_node(aux_dict)

for idx, row in df.iterrows():
    directed = 'false'
    direct = 'flase'
    print(row.to_dict())
    if row.is_directed == 1:
        directed = 'true'
        direct = 'true'
    elif row.is_directed == 0:
        directed = 'false'
    else:
        print("WARNING: unknown direction flag in line: " + idx)
    interaction_types = "is_directed:%s|is_direct:%s" % (directed, direct)
    edge_dict = {
        'source_db': 'OmniPath',
        'interaction_types': interaction_types,
        'layer': 0,
    }
    source_dict = db_api.get_node(row.source, 9606)
    target_dict = db_api.get_node(row.target, 9606)
    db_api.insert_edge(source_dict, target_dict, edge_dict)

#layer 1, direct regulators of core prs
df1 = op_interact[
    (op_interact.target.isin(core_df.UniProtAC)) &
    (~op_interact.source.isin(core_df.UniProtAC))]

layer1_nodes = pd.DataFrame(columns=core_df.columns)
layer1_nodes.columns
layer1_nodes['UniProtAC'] = df1.source
#fncs = []
#for idx, row in layer1_nodes.iterrows():
#    if "_" in row.UniProtAC:
#        continue
#    print(row.UniProtAC)
#    resp = uniprot_client.get_pr_fnc(row.UniProtAC)
#    fncs.append(resp)
#layer1_nodes['pr_funcs'] = fncs

for idx, row in layer1_nodes.iterrows():
    aux_dict = dict()
    aux_dict['name'] = row.UniProtAC
    aux_dict['gene_name'] = "-"
    aux_dict['tax_id'] = 9606
    aux_dict['pathways'] = "-"
    aux_dict['source'] = 'OmniPath'
#    cleaned = '|'.join(item.strip().strip("'") for item in row.molecular_function.strip('[]').split(','))
    aux_dict['function'] = "-"
    db_api.insert_node(aux_dict)

for idx, row in df1.iterrows():
    directed = 'false'
    direct = 'flase'
    print(row.to_dict())
    if row.is_directed == 1:
        directed = 'true'
        direct = 'true'
    elif row.is_directed == 0:
        directed = 'false'
    else:
        print("WARNING: unknown direction flag in line: " + idx)
    interaction_types = "is_directed:%s|is_direct:%s" % (directed, direct)
    edge_dict = {
        'source_db': 'OmniPath',
        'interaction_types': interaction_types,
        'layer': 1,
    }
    source_dict = db_api.get_node(row.source, 9606)
    target_dict = db_api.get_node(row.target, 9606)
    db_api.insert_edge(source_dict, target_dict, edge_dict)
db_api.save_db_to_file(DB_DESTINATION)
