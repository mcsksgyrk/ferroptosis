import pandas as pd
from database.external_db import DBconnector
from config import OUTPUTS_DIR, SOURCES_DIR


def get_uniprot_id(row, db):
    target_id_query = f"SELECT uniprot_id FROM general_target WHERE target_id = '{row.target_id}'"
    regulator_id_query = f"SELECT External_id FROM general_regulator WHERE regulator_id = '{row.regulator_id}'"
    target_id = db.custom_query(target_id_query)
    regulator_id = db.custom_query(regulator_id_query)
    return [flatten_if_single_list(target_id), flatten_if_single_list(regulator_id)]


def flatten_if_single_list(value):
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


ferreg = DBconnector(OUTPUTS_DIR / "ferreg.db")
omnipath = pd.read_csv(SOURCES_DIR / "omnipath" / "omnipath_interactions.txt", delim_whitespace=True)

get_edges = "SELECT target_id, regulator_id FROM target_regulator_drug_disease_pair"
ferreg_edges = ferreg.query_to_dataframe(get_edges)
ferreg_edges = ferreg_edges[~((ferreg_edges.target_id=='TAR99999') | (ferreg_edges.regulator_id=="."))]
res = ferreg_edges.apply(lambda row: get_uniprot_id(row, ferreg), axis=1)
results_df = pd.DataFrame(res.tolist(), index=ferreg_edges.index, columns=['target', 'source'])
ferreg_edges = pd.concat([ferreg_edges, results_df], axis=1)

core_prs = pd.read_csv(OUTPUTS_DIR / "core_geneProducts.csv", delimiter=",").uniprot_id.tolist()

wtf = []
layer0 = []
layer1 = []
for id, row in ferreg_edges.iterrows():
    if row.target and row.source in core_prs:
        layer0.append(row)
    elif row.target in core_prs and row.source not in core_prs:
        layer1.append(row)
    elif row.target not in core_prs:
        if row.source in core_prs:
            print('?')
        wtf.append(row)
len(wtf)
len(layer0)
len(layer1)
wtf
