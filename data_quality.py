from config import OUTPUTS_DIR, SOURCES_DIR, Path
from database.external_db import DBconnector
import re
from apicalls.uniprot import UniProtClient
import sqlite3
from typing import List
import pickle


class TestInterface:
    def __init__(self, db_path: str):
        self.db_path = db_path

    def _get_connection(self):
        return sqlite3.connect(self.db_path)

    def update_entry(self, query: str):
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query)

    def custom_query(self, query: str) -> List:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            res = cursor.execute(query).fetchall()
            if not res:
                return []
            if len(res[0]) == 1:
                return [row[0] for row in res]
            else:
                return res


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
db_path = OUTPUTS_DIR / "ferr_test.db"
db = DBconnector(db_path)

if Path("filename.pickle").exists():
    uniprot = UniProtClient()
    with open('filename.pickle', 'rb') as handle:
        existing_gn_up_pairs = pickle.load(handle)

query = """
    SELECT n.*
    FROM node n
    WHERE n.source_db LIKE "%ferrdb%"
    AND n.primary_id_type = "uniprot_id"
"""

uniprot_but_not_pr = []
pr_w_invalid_uniprot = []
res = db.query_to_dataframe(query)
for idx, row in res.iterrows():
    if is_uniprot_id(row['name']):
        continue
    else:
        pr_w_invalid_uniprot.append(row['name'])
pr_w_invalid_uniprot.remove('_NA_')

if len(pr_w_invalid_uniprot) != 0:
    r, f = uniprot.batch_convert_to_uniprot_id("Gene_Name", pr_w_invalid_uniprot, human=True)
    if not existing_gn_up_pairs:
        with open('filename.pickle', 'wb') as handle:
            pickle.dump(r, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        for k, v in r.items():
            existing_gn_up_pairs[k] = v

test_db = TestInterface(db_path)
for geneName, uniprotID in existing_gn_up_pairs.items():
    existing_uniprot_entry = test_db.custom_query(f"SELECT id, name FROM node WHERE name = '{uniprotID}'")
    if existing_uniprot_entry:
        print(f"{geneName} id {uniprotID} exists")
        continue
    print(f"Updating {geneName} to {uniprotID}")
    query = f"""
        UPDATE node
        SET name = '{uniprotID}',
            primary_id_type = 'uniprot_id'
        WHERE name = '{geneName}'
    """
    test_db.update_entry(query)
    query = f"""
        UPDATE edge
        SET interactor_a_node_name = '{uniprotID}'
        WHERE interactor_a_node_name = '{geneName}'
    """
    test_db.update_entry(query)

    query = f"""
        UPDATE edge
        SET interactor_b_node_name = '{uniprotID}'
        WHERE interactor_b_node_name = '{geneName}'
    """
    test_db.update_entry(query)

for geneName, uniprotID in r.items():
    query = f"""
        SELECT n.id, n.name, ni.id_type, ni.id_value, ni.is_primary
        FROM node n
        LEFT JOIN node_identifier ni ON n.id = ni.node_id
        WHERE n.name = '{uniprotID}'
    """
    results = test_db.custom_query(query)
    print(results)

    if not results:
        print(f"Node with name {uniprotID} not found")
        continue

nd_query = """
    SELECT name FROM node
    WHERE type = 'nd'
"""

checking_bro = test_db.custom_query(nd_query)
uniprot = UniProtClient()
r, f = uniprot.batch_convert_to_uniprot_id("GeneCard", checking_bro, human=True)
uniprot.batch_convert_to_uniprot_id("Gene_Name", ['HERC2'], human=True)

res[res['name']=='PRAP1']
is_uniprot_id('PRAP1')
uniprot.convert_to_uniprot_id('Gene_Name', ['PRAP1'], True)
