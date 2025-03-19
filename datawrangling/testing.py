import pandas as pd
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from database.sqlite_db_api2 import PsimiSQL


SQL_SEED = PROJECT_ROOT/"database"/"network_db_seed2.sql"
db_api = PsimiSQL(SQL_SEED)


db_api.import_from_db_file(str(OUTPUTS_DIR/"ferrdb_network.db"))

node = db_api.get_node_by_id(1)
test = {'name': 'P62917', 'primary_id_type': 'uniprot_id', 'display_name': 'RPL8',
        'tax_id': 9606, 'type': 'protein', 'pathways': '',
        'source': 'FerrDB_driver_prs',
        'function': '', 'uniprot_id': 'P62917', 'gene_name': 'RPL8',
        'hgnc_id': 'HGNC:10368', 'ensembl_id': 'ENSG00000161016'}

test1 = {'name': 'P00100', 'primary_id_type': 'uniprot_id', 'display_name': 'RPL8',
        'tax_id': 9606, 'type': 'protein', 'pathways': '',
        'source': 'FerrDB_driver_prs',
        'function': '', 'uniprot_id': 'P62917', 'gene_name': 'RPL8',
        'hgnc_id': 'HGNC:10368', 'ensembl_id': 'ENSG00000161016'}

test2 = {'name': 'ENSG00000161016', 'primary_id_type': 'uniprot_id', 'display_name': 'RPL8',
        'tax_id': 9606, 'type': 'protein', 'pathways': '',
        'source': 'FerrDB_driver_prs',
        'function': '', 'uniprot_id': 'P62917', 'gene_name': 'RPL8',
        'hgnc_id': 'HGNC:10368', 'ensembl_id': 'ENSG00000161016'}
db_api.check_if_node_exists(test)
db_api.check_if_node_exists(test1)
db_api.check_if_node_exists(test2)
db_api.insert_node(test1)
db_api.get_node_by_name('P00100')
