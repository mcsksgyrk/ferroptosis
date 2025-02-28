import pandas as pd
from database.external_db import DBconnector
from config import OUTPUTS_DIR


def get_uniprot_id(row, db):
    target_id_query = f"SELECT uniprot_id FROM general_target WHERE target_id = '{row.target_id}'"
    regulator_id_query = f"SELECT External_id FROM general_regulator WHERE regulator_id = '{row.regulator_id}'"
    target_id = db.custom_query(target_id_query)
    regulator_id = db.custom_query(regulator_id_query)
    return [target_id, regulator_id]


ferreg = DBconnector(OUTPUTS_DIR / "ferreg.db")

get_edges = "SELECT target_id, regulator_id FROM target_regulator_drug_disease_pair"
ferreg_edges = ferreg.query_to_dataframe(get_edges)
ferreg_edges = ferreg_edges[~((ferreg_edges.target_id=='TAR99999') | (ferreg_edges.regulator_id=="."))]
res = ferreg_edges.apply(lambda row: get_uniprot_id(row, ferreg), axis=1)
results_df = pd.DataFrame(res.tolist(), index=ferreg_edges.index, columns=['target', 'source'])
ferreg_edges = pd.concat([ferreg_edges, results_df], axis=1)
ferreg_edges
