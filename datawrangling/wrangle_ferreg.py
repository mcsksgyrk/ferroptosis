import pandas as pd
from config import SOURCES_DIR, OUTPUTS_DIR
import sqlite3


def make_sql_query(what, table, where=None, values=None, db="ferreg.db"):
    with sqlite3.connect(OUTPUTS_DIR / db) as conn:
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


def compile_query_string(columns, table, filters=None) -> str:
    query = f"SELECT {columns} FROM {table}"
    if filters is not None:
        query += " WHERE "
        for filter in filters:
            if filter == filters[-1]:
                query += filter
            else:
                query += filter + " AND "
    return query


def lists_difference(list1, list2):
    list1_only = list(set(list1)-set(list2))
    list2_only = list(set(list2)-set(list1))
    return [list1_only, list2_only]


def get_caonincal_prs(pr_list):
    return [x for x in pr_list if x[0] in ['P', 'Q', 'O']]


core_prs = pd.read_csv(OUTPUTS_DIR / "core_geneProducts.csv")

ferreg_target_in_cores = make_sql_query("uniprot_id",
                                        "general_target",
                                        "uniprot_id",
                                        core_prs['uniprot_id'])
ferreg_targets = make_sql_query("uniprot_id",
                                "general_target")

filters = [
           "Exp_Organism='Human'",
           "Gene_type_hgnc_locus_type_or_other='gene with protein product'",
           "Confidence='Validated'"
           ]


query = compile_query_string("DISTINCT UniProtAC", "driver", filters)
res = custom_query(query, OUTPUTS_DIR / "ferrdb.db")
len(res)
a, b = lists_difference(core_prs.uniprot_id.tolist(), res)
hc = core_prs.uniprot_id.tolist()
len(set(hc)&set(res))
