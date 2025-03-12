
import pandas as pd
from typing import Dict, List
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from apicalls.api_oop import KEGGClient, UniProtClient
from database.sqlite_db_api2 import PsimiSQL
from database.external_db import DBconnector


ferrdb = DBconnector(OUTPUTS_DIR / "ferrdb.db")

cols = ["UniProtAC","Symbol_or_reported_abbr","Exp_organism"]
driver_prs_q = f"SELECT {', '.join(cols)} FROM driver WHERE Gene_type_hgnc_locus_type_or_other='gene with protein product'"
driver_prs = ferrdb.custom_query(driver_prs_q)
pr_df_cols = ['unbiprot_id', 'gene_name', 'topology']
df_drivers = pd.DataFrame(data=driver_prs, columns=pr_df_cols)

ferrdb.get_table_names()
