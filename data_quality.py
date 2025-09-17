from config import OUTPUTS_DIR, SOURCES_DIR
from database.external_db import DBconnector
import re
from apicalls.uniprot import UniProtClient

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


def is_uniprot_id(identifier):
    if identifier is None:
        return False
    identifier_str = str(identifier).upper().strip()
    swiss_prot = r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])$'
    trembl = r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]$'
    return bool(re.match(swiss_prot, identifier_str) or re.match(trembl, identifier_str))


cpd_path = SOURCES_DIR / "kegg/kegg_compounds.txt"
drug_path = SOURCES_DIR / "kegg/kegg_drugs.txt"
compounds, revc = make_compound_dict(cpd_path)
drugs, revd = make_compound_dict(drug_path)
db_path = OUTPUTS_DIR / "final.db"
db = DBconnector(db_path)

query = """
    SELECT n.*
    FROM node n
    WHERE n.source_db LIKE "%ferrdb%"
    AND n.primary_id_type = "uniprot_id"
"""

false_uniprots = []
res = db.query_to_dataframe(query)
for idx, row in res.iterrows():
    if is_uniprot_id(row['name']):
        continue
    else:
        false_uniprots.append(row['name'])

uniprot = UniProtClient()
for id in false_uniprots:
    uniprot.convert_to_uniprot_id("Gene_Name", [id])
