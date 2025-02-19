import pandas as pd
import sqlite3


def make_sql_query(what, table, where=None, values=None, db="ferreg.db"):
    with sqlite3.connect(db) as conn:
        if values is not None and where is not None:
            placeholders = ",".join("?"*len(values))
            query = f'SELECT {what} FROM {table} WHERE {where} IN ({placeholders})'
            cursor = conn.cursor()
            res = cursor.execute(query, values.tolist()).fetchall()
            res = [pr[0] for pr in res]
        else:
            query = f'SELECT {what} FROM {table}'
            cursor = conn.cursor()
            res = cursor.execute(query).fetchall()
            res = [pr[0] for pr in res]
    return res


def custom_query(query, db):
    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        res = cursor.execute(query)
        res = [s[0] for s in res]
    return res


core_prs = pd.read_csv("./prlists/prs_from_kegg_go_wp.csv")
all_core = pd.read_csv("./prlists/core_proteins.csv")
osszes = pd.read_csv("./prlists/all_col_core_df.csv")

res_coree = make_sql_query("uniprot_id", '"target_data|general_target"', "uniprot_id", core_prs.uniprotID, "ferreg.db")
res_layer = make_sql_query("uniprot_id", '"target_data|general_target"', "uniprot_id", all_core.UniProtAC, "ferreg.db")
res_all = make_sql_query("uniprot_id", '"target_data|general_target"')

# TFR2 and GPX1 only in ferreg
only_in_ferreg = list(set(res_all)-(set(res_all) & set(res_layer)))
only_in_core = list((set(all_core.UniProtAC.tolist())-set(res_all)))

core_prs[~core_prs.uniprotID.isin(only_in_core)]

ferreg_regulators = custom_query('SELECT External_id FROM "disease_drug|general_regulator" WHERE type = "Protein coding"', "ferreg.db")
core_regulators = all_core.UniProtAC.tolist()
osszes_regulators = osszes.UniProtAC.tolist()
intersction = list(set(ferreg_regulators) & set(core_regulators))
sajat_reg = list(set(core_prs.uniprotID.tolist()) & set(osszes_regulators))
set(core_prs.uniprotID.tolist())-set(sajat_reg)
core_prs
