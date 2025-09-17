import pandas as pd
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from apicalls.pubchem import PubChemClient
from database.external_db import DBconnector


def make_compound_set(f_path):
    compounds = set()
    with open(f_path, 'r') as f:
        for line in f:
            kegg_entry = line.strip().split('\t', 1)[1]
            names = [name.strip().lower() for name in kegg_entry.split(';')]
            compounds.update(names)
        return compounds


db_path = OUTPUTS_DIR / "final.db"
f_path = SOURCES_DIR / "kegg/kegg_compounds.txt"
db = DBconnector(db_path)

query = """
    SELECT n.*, ni.*
    FROM node n
    INNER JOIN node_identifier as ni ON n.id = ni.node_id
    WHERE n.type IN ('compound', 'nd', 'small_molecule')
"""
compound_df = db.query_to_dataframe(query)
kegg_compounds = make_compound_set(f_path)
pubchem = PubChemClient()
unique_compounds = set()
duplicated_names = set()
for idx, row in compound_df.iterrows():
    inchi = pubchem.name_to_inchikey(row['display_name'])
    if inchi not in unique_compounds:
        unique_compounds.add(inchi)
    elif inchi in unique_compounds:
        duplicated_names.add(row['display_name'])
len(duplicated_names)
