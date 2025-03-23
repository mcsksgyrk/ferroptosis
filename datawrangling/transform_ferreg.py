from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from typing import List, Union, Optional, Any
import pandas as pd
from database.external_db import DBconnector
from database.sqlite_db_api2 import PsimiSQL
from apicalls.api_oop import PubChemClient


def extractor(db: DBconnector,
              columns: str,
              table: str,
              filters: Optional[List[str]] = None) -> List[str]:
    query = db.compile_query_string(columns, table, filters)
    res = db.custom_query(query)
    return res


ferreg = DBconnector(OUTPUTS_DIR / "ferreg.db")
query_uniques = [
    {"db": ferreg,
     "table": "general_drug",
     "id_column": "drug_name,drug_type,ttd_id,drugmap_id",
     "params": ["Type='Protein coding'"],
     "source": "ferreg",
     "is_core": 0
     }
]
ligand_cols_name = query_uniques[0]['id_column'].split(',')
res = extractor(ferreg, query_uniques[0]["id_column"], query_uniques[0]['table'])
ligand_df = pd.DataFrame(data=res, columns=ligand_cols_name)

pubchem_client = PubChemClient()
for _, row in ligand_df.iterrows():
    for column, val in row.items():
        res = pubchem_client.get_primary_cid_or_sid(val)
        if isinstance(res, int):
            print(res)
            break

SQL_SEED = PROJECT_ROOT/"database"/"network_db_seed2.sql"
DB_DESTINATION = OUTPUTS_DIR/"ferreg_network.db"
ferreg = DBconnector(OUTPUTS_DIR/"ferrdb.db")
db_api = PsimiSQL(SQL_SEED)
