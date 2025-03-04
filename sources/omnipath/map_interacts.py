import apicalls.api_oop as api
import pandas as pd
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from database.sqlite_db_api import PsimiSQL

op_interact = pd.read_csv(SOURCES_DIR / "omnipath" /"omnipath_interactions.txt", delim_whitespace=True)
core_df = pd.read_csv(OUTPUTS_DIR / "nodes_NO_op.csv")


core_df['tax_id'] = 9606
reactome_client = api.ReactomeClient()
pathways = []
for idx, row in core_df.iterrows():
    print(row.uniprot_id)
    pathways.append(reactome_client.map_protein_to_pathways(row.uniprot_id))

core_df['pathways'] = pathways
core_df['pathways'] = core_df['pathways'].apply(lambda x: ';'.join(x))
core_df.to_csv(OUTPUTS_DIR / "nodes_w_pw.csv")

SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed.sql"
DB_DESTINATION = OUTPUTS_DIR / "test.db"
db_api = PsimiSQL(SQL_SEED)

df = op_interact[
    (op_interact.source.isin(core_df.uniprot_id)) & (op_interact.target.isin(core_df.uniprot_id))]
for idx, row in core_df.iterrows():
    aux_dict = dict()
    aux_dict['name'] = row.uniprot_id
    aux_dict['gene_name'] = ""
    aux_dict['tax_id'] = row.tax_id
    aux_dict['pathways'] = row.pathways
    aux_dict['source'] = row.source
    aux_dict['function'] = "kekking"
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
        'layer': 1,
    }
    source_dict = db_api.get_node(row.source, 9606)
    target_dict = db_api.get_node(row.target, 9606)
    db_api.insert_edge(source_dict, target_dict, edge_dict)
db_api.save_db_to_file(str(DB_DESTINATION))
