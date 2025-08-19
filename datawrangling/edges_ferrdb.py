
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
from database.external_db import DBconnector
from apicalls.uniprot import UniProtClient
import re


def parse_ferrdb_pw(row):
    if row.Pathway == '_NA_':
        return []
    else:
        edge_dict = []
        edges = row['Pathway'].split(',')
        for edge in edges:
            nodes = re.split(r':-:|:\+:', edge)
            d = {
                'interaction_type': check_edge_sign(edge),
                'interactor_a_node_name': nodes[0],
                'interactor_b_node_name': nodes[1]
            }
            edge_dict.append(d)
        return edge_dict


def check_edge_sign(edge):
    if ':-:' in edge:
        return 'down'
    if ':+:' in edge:
        return 'up'


ferrdb = DBconnector(OUTPUTS_DIR / "ferrdb.db")
uniprot = UniProtClient()

nodes_query = """
SELECT Symbol_or_reported_abbr, UniProtAC, Pathway, PMID FROM driver
WHERE Confidence = 'Validated'
AND LOWER(Exp_organism) LIKE '%human%'
AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
UNION
SELECT Symbol, UniProtAC, Pathway, PMID FROM suppressor
WHERE Confidence = 'Validated'
AND LOWER(Exp_organism) LIKE '%human%'
AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
UNION
SELECT Symbol, UniProtAC, Pathway, PMID FROM marker
WHERE Confidence = 'Validated'
AND LOWER(Exp_organism) LIKE '%human%'
AND Gene_type_hgnc_locus_type_or_other = 'gene with protein product'
"""
nodes_df = ferrdb.query_to_dataframe(nodes_query)

nodes_df['edges'] = nodes_df.apply(parse_ferrdb_pw, axis=1)
single_step_edges = nodes_df[nodes_df['edges'].apply(len)<2]
multiple_step_edges = nodes_df[nodes_df['edges'].apply(len)>1]
single_step_edges
multiple_step_edges
