import pandas as pd
from config import SOURCES_DIR, OUTPUTS_DIR
from database.external_db import DBconnector
from typing import List, Union, Optional, Any


def lists_difference(list1, list2):
    list1_only = list(set(list1)-set(list2))
    list2_only = list(set(list2)-set(list1))
    return [list1_only, list2_only]


def lists_intersection(list1, list2):
    return list(set(list1) & set(list2))


def get_caonincal_prs(pr_list):
    return [x for x in pr_list if x[0] in ['P', 'Q', 'O']]


def extend_source_column(val, source):
    ext_val = val+";"+source
    return ext_val


def extractor(db: DBconnector,
              columns: str,
              table: str,
              filters: Optional[List[str]] = None) -> List[str]:
    query = db.compile_query_string(columns, table, filters)
    res = db.custom_query(query)
    return res


core_prs = pd.read_csv(OUTPUTS_DIR / "core_geneProducts.csv")

ferreg = DBconnector(OUTPUTS_DIR / "ferreg.db")
ferreg_target_in_core = ferreg.query("uniprot_id",
                                     "general_target",
                                     "uniprot_id",
                                     core_prs.uniprot_id.tolist())
ferreg_target = ferreg.query("uniprot_id",
                             "general_target")
ferrdb = DBconnector(OUTPUTS_DIR / "ferrdb.db")
ferrdb_driver_in_core = ferrdb.query("UniProtAC",
                                     "marker",
                                     "UniProtAC",
                                     core_prs.uniprot_id.tolist())
ferrdb_suppressor_in_core = ferrdb.query("UniProtAC",
                                         "suppressor",
                                         "UniProtAC",
                                         core_prs.uniprot_id.tolist())
in_ferrdb = list(set(ferrdb_suppressor_in_core).union(ferrdb_driver_in_core))
ext_df = core_prs.copy()
ext_df.loc[ext_df['uniprot_id'].isin(ferreg_target_in_core), 'source'] += ";ferreg"
ext_df.loc[ext_df['uniprot_id'].isin(in_ferrdb), 'source'] += ";ferrdb"

"""
ferreg regulator:
    Type = Protein coding
    External_id az uniprot_id...
"""
params = [
          "Type='Protein coding'",
]

query = ferreg.compile_query_string("External_id", "general_regulator", params)
ferreg_regulator = ferreg.custom_query(query)
