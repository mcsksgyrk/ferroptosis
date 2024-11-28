import apicalls.api_oop as api
import pandas as pd
from database.sqlite_db_api import PsimiSQL

op_interact = pd.read_csv("./sources/omnipath/omnipath_pr_interacts.txt", delim_whitespace=True)
core_df = pd.read_csv('./prlists/all_col_core_df.csv')

reactome_client = api.ReactomeClient()

SQL_SEED = "./database/network_db_seed.sql"
DB_DESTINATION = "./output/omnipath"
db_api = PsimiSQL(SQL_SEED)

df = op_interact[
    (op_interact.source.isin(core_df.UniProtAC)) & (op_interact.target.isin(core_df.UniProtAC))]

for idx, row in core_df.iterrows():
    aux_dict = dict()
    aux_dict['name'] = row.UniProtAC
    aux_dict['gene_name'] = row.Gene_Name
    aux_dict['tax_id'] = row.tax_id
    aux_dict['pathways'] = row.pathways
    db_api.insert_node(aux_dict)
df.columns
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
    }
    source_dict = db_api.get_node(row.source, 9606)
    target_dict = db_api.get_node(row.target, 9606)
    db_api.insert_edge(source_dict, target_dict, edge_dict)
db_api.save_db_to_file(DB_DESTINATION)
