from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from typing import List, Union, Optional, Any
import pandas as pd
from database.external_db import DBconnector
from database.sqlite_db_api2 import PsimiSQL
from apicalls.api_oop import PubChemClient


SQL_SEED = PROJECT_ROOT/"database"/"network_db_seed2.sql"
DB_DESTINATION = OUTPUTS_DIR/"ferreg_network.db"
ferreg = DBconnector(OUTPUTS_DIR/"ferrdb.db")
db_api = PsimiSQL(SQL_SEED)
pubchem_client = PubChemClient()
pubchem_client.get_primary_cid_or_sid("Tongxinluo")
