import pandas as pd
from database.external_db import DBconnector
from config import OUTPUTS_DIR, SOURCES_DIR
from typing import List, Union, Optional, Any


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


def combine_sources(dataframes_list: list[pd.DataFrame]) -> pd.DataFrame:
    combined_df = pd.concat(dataframes_list)
    combined_df = combined_df.drop_duplicates(subset=['uniprot_id', 'is_core', 'source'])
    grouped = combined_df.groupby(['uniprot_id', 'is_core'])
    result_df = grouped.agg({
        'source': lambda x: ';'.join(sorted(set(x)))
    }).reset_index()
    return result_df


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
query_uniques = [
    {"db": ferreg,
     "table": "general_regulator",
     "id_column": "External_id",
     "params": ["Type='Protein coding'"],
     "source": "ferreg",
     "is_core": 0
     },
    {"db": ferrdb,
     "table": "suppressor",
     "id_column": "UniProtAC",
     "params": ["Gene_type_hgnc_locus_type_or_other='gene with protein product'"],
     "source": "ferrdb",
     "is_core": 0
     },
    {"db": ferrdb,
     "table": "marker",
     "id_column": "UniProtAC",
     "params": ["Gene_type_hgnc_locus_type_or_other='gene with protein product'"],
     "source": "ferrdb",
     "is_core": 0
     },
    {"db": ferrdb,
     "table": "driver",
     "id_column": "UniProtAC",
     "params": ["Gene_type_hgnc_locus_type_or_other='gene with protein product'"],
     "source": "ferrdb",
     "is_core": 0
     }
]
core_ids = core_prs.uniprot_id.tolist()
df_list = []
for q in query_uniques:
    res = extractor(q['db'], q['id_column'], q['table'], q['params'])
    res_unique = lists_difference(res, core_ids)[0]
    df = pd.DataFrame(
        {"uniprot_id": res_unique,
         'source': q['source'],
         'is_core': q['is_core']})
    df_list.append(df)

not_core_df = combine_sources(df_list)
nodes_df = pd.concat([core_prs, not_core_df])
nodes_df = nodes_df[~nodes_df.uniprot_id.str.contains("_")]
nodes_df.to_csv(OUTPUTS_DIR / "nodes_NO_op.csv", index=False)

ferreg = DBconnector(OUTPUTS_DIR / "ferreg.db")

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
