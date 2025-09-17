import pandas as pd
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from apicalls.pubchem import PubChemClient
from database.external_db import DBconnector


def make_kegg_dict(f):
    kegg_dict = dict()
    rev_dict = dict()
    for line in f:
        kegg_line = line.strip().split('\t', 1)
        if len(kegg_line) < 2:
            continue
        kegg_entry = kegg_line[0]
        names = [name.strip().lower() for name in kegg_line[1].split(';')]
        for name in names:
            rev_dict[name] = kegg_entry
        kegg_dict[kegg_entry] = names
    return kegg_dict, rev_dict


def make_compound_dict(f_path):
    with open(f_path, 'r') as f:
        kegg_dict, rev_dict = make_kegg_dict(f)
    return kegg_dict, rev_dict


db_path = OUTPUTS_DIR / "final.db"
cpd_path = SOURCES_DIR / "kegg/kegg_compounds.txt"
drug_path = SOURCES_DIR / "kegg/kegg_drugs.txt"
db = DBconnector(db_path)

query = """
    SELECT n.*, ni.*
    FROM node n
    INNER JOIN node_identifier as ni ON n.id = ni.node_id
    WHERE n.type IN ('compound', 'nd', 'small_molecule')
"""
compound_df = db.query_to_dataframe(query)
compounds, revc = make_compound_dict(cpd_path)
drugs, revd = make_compound_dict(drug_path)
shiet = set()
for idx, row in compound_df.iterrows():
    bby = row['display_name'].lower()
    if revc.get(bby).source_db:
        continue
    elif revd.get(bby):
        continue
    else:
        shiet.add(bby)
